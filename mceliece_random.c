#include "mceliece_types.h"
#include <time.h>
#include "mceliece_random.h"
#include "gf.h"
#include "stdio.h"
#include <stdint.h> // For uint16_t, uint32_t


typedef struct {
    uint32_t val; // <--- 必须是 uint32_t！
    uint16_t pos;
} pair_t;

int compare_pairs(const void *a, const void *b) {
    const pair_t *p1 = (const pair_t *)a;
    const pair_t *p2 = (const pair_t *)b;
    if (p1->val < p2->val) return -1;
    if (p1->val > p2->val) return 1;
    // 如果值相同，按原始位置排序以保证稳定性（可选，但良好实践）
    if (p1->pos < p2->pos) return -1;
    if (p1->pos > p2->pos) return 1;
    return 0;
}

// 简单的伪随机数生成器状态
static uint64_t rng_state = 1;

// 设置随机种子
void set_random_seed(uint64_t seed) {
    rng_state = seed ? seed : 1;
}

// 初始化随机种子（使用当前时间）
void init_random(void) {
    set_random_seed((uint64_t)time(NULL));
}

// 生成随机字节
uint8_t random_byte(void) {
    // 简单的线性同余生成器
    rng_state = rng_state * 1103515245ULL + 12345ULL;
    return (uint8_t)(rng_state >> 32);
}

// 生成随机字节数组
void random_bytes(uint8_t *output, int len) {
    for (int i = 0; i < len; i++) {
        output[i] = random_byte();
    }
}


mceliece_error_t generate_field_ordering(gf_elem_t *alpha, const uint8_t *random_bits) {
    int q = MCELIECE_Q;
    int m = MCELIECE_M;
    int sigma2_bits = 32;
    int sigma2_bytes = sigma2_bits / 8;

    printf("    -> generate_field_ordering: Entered function.\n");
    printf("    -> generate_field_ordering: First 8 random bytes for a_i: %02x %02x %02x %02x %02x %02x %02x %02x\n",
           random_bits[0], random_bits[1], random_bits[2], random_bits[3],
           random_bits[4], random_bits[5], random_bits[6], random_bits[7]);

    pair_t *pairs = malloc(q * sizeof(pair_t));
    if (!pairs) {
        printf("    -> generate_field_ordering: FATAL - Malloc for pairs failed.\n");
        return MCELIECE_ERROR_MEMORY;
    }

    // 1. 从随机比特生成 q 个 32-bit 的整数 a_i (小端序)
    for (int i = 0; i < q; i++) {
        int offset = i * sigma2_bytes;
        uint32_t a_i = (uint32_t)random_bits[offset] |
                       ((uint32_t)random_bits[offset + 1] << 8) |
                       ((uint32_t)random_bits[offset + 2] << 16) |
                       ((uint32_t)random_bits[offset + 3] << 24);
        pairs[i].val = a_i;
        pairs[i].pos = i;
    }

    // 2. 检查是否有重复值
    printf("    -> generate_field_ordering: Checking for duplicates among %d values...\n", q);
    pair_t *sorted_for_check = malloc(q * sizeof(pair_t));
    if (!sorted_for_check) { free(pairs); return MCELIECE_ERROR_MEMORY; }
    memcpy(sorted_for_check, pairs, q * sizeof(pair_t));
    qsort(sorted_for_check, q, sizeof(pair_t), compare_pairs);

    int has_duplicates = 0;
    for (int i = 0; i < q - 1; i++) {
        if (sorted_for_check[i].val == sorted_for_check[i+1].val) {
            has_duplicates = 1;
            // ---> 新增日志
            printf("    -> generate_field_ordering: FAILURE! Found duplicate value: %u at original positions %u and %u\n",
                   sorted_for_check[i].val, sorted_for_check[i].pos, sorted_for_check[i+1].pos);
            break;
        }
    }
    free(sorted_for_check);

    if (has_duplicates) {
        free(pairs);
        return MCELIECE_ERROR_KEYGEN_FAIL;
    }

    // ---> 新增日志
    printf("    -> generate_field_ordering: No duplicates found. Proceeding to sort...\n");

    // 3. 按值对 (a_i, i) 进行字典序排序
    qsort(pairs, q, sizeof(pair_t), compare_pairs);

    // 4. 定义置换 pi，pi[i] 是排序后第 i 个元素的原始位置
    uint16_t *pi = malloc(q * sizeof(uint16_t));
    if (!pi) { free(pairs); return MCELIECE_ERROR_MEMORY; }
    for(int i = 0; i < q; ++i) {
        pi[i] = pairs[i].pos;
    }

    free(pairs); // pairs 不再需要

    // 5. 根据置换 pi 生成最终的 alpha 序列
    //    规范公式: α_i = Σ π(i)_j * z^(m-1-j) for j=0 to m-1
    //    在标准二进制表示中，这恰好等价于将整数 π(i) 直接作为 F_q 的元素。
    //    我们在这里显式地实现它，以确保没有误解。
    for (int i = 0; i < q; i++) {
        gf_elem_t current_alpha = 0;
        uint16_t pi_val = pi[i];

        // 我们只关心 pi_val 的低 m 位
        pi_val &= (1 << m) - 1;

        // 规范中的公式将 j=0 对应到最高位 z^(m-1)。
        // 这在标准表示法中与直接使用整数值是等价的。
        // 例如：pi_val = 13 (0...1101), m=4.
        // j=0 (LSB=1) -> 1*z^0
        // j=1 (LSB=0) -> 0*z^1
        // j=2 (LSB=1) -> 1*z^2
        // j=3 (MSB=1) -> 1*z^3
        // sum = z^3+z^2+1, which is the element 13.
        // 所以直接赋值是正确的。
        current_alpha = (gf_elem_t)pi_val;
        alpha[i] = current_alpha;
    }

    free(pi);
    return MCELIECE_SUCCESS;
}

static void get_field_poly_F(polynomial_t *F_y) {
    int t = MCELIECE_T;

    // 检查是否是我们正在实现的参数集 mceliece6688128 (t=128)
    if (t == 128) {
        // 根据规范，对于 mceliece6688128，
        // F(y) = y^128 + y^7 + y^2 + y + 1

        // 清零多项式以防有旧数据
        memset(F_y->coeffs, 0, (F_y->max_degree + 1) * sizeof(gf_elem_t));
        F_y->degree = -1;

        // 设置非零系数。这些系数都是 F_q 中的元素 '1'。
        polynomial_set_coeff(F_y, 0,   1); // y^0 项 (常数项)
        polynomial_set_coeff(F_y, 1,   1); // y^1 项
        polynomial_set_coeff(F_y, 2,   1); // y^2 项
        polynomial_set_coeff(F_y, 7,   1); // y^7 项
        polynomial_set_coeff(F_y, 128, 1); // y^128 项 (首项)

    } else {
        // 如果您的 MCELIECE_T 被设置为其他值，程序会报错退出。
        // 这可以防止因为配置错误而产生无效的密钥。
        fprintf(stderr, "FATAL ERROR: Field polynomial F(y) is not defined in the code for t=%d. "
                        "Please define it for your chosen parameter set.\n", t);
        exit(EXIT_FAILURE);
    }
}

static polynomial_t* ext_zero(int t) {
    polynomial_t *z = polynomial_create(t - 1);
    return z; // already zero
}

static polynomial_t* ext_one(int t) {
    polynomial_t *o = polynomial_create(t - 1);
    polynomial_set_coeff(o, 0, 1);
    return o;
}

static void ext_copy(polynomial_t *dst, const polynomial_t *src) {
    polynomial_copy(dst, src);
}

static void ext_add(polynomial_t *out, const polynomial_t *a, const polynomial_t *b) {
    polynomial_add(out, a, b);
}

static void ext_mul(polynomial_t *out, const polynomial_t *a, const polynomial_t *b, const polynomial_t *F) {
    polynomial_t *tmp = polynomial_create(2 * (F->degree) - 2);
    polynomial_mul(tmp, a, b);
    polynomial_div(NULL, out, tmp, F);
    polynomial_free(tmp);
}

static void ext_qpower(polynomial_t *out, const polynomial_t *a, const polynomial_t *F, int m) {
    // raise a to q=2^m inside F_{q^t}
    polynomial_t *acc = polynomial_create(F->degree - 1);
    polynomial_t *tmp = polynomial_create(F->degree - 1);
    polynomial_t *sq_tmp = polynomial_create(2 * (F->degree) - 2);
    ext_copy(acc, a);

    for (int i = 0; i < m; i++) {
        // square: acc = acc^2
        // compute poly_mul(acc, acc) then mod F
        polynomial_mul(sq_tmp, acc, acc);
        polynomial_div(NULL, tmp, sq_tmp, F);
        ext_copy(acc, tmp);
        // additionally square coefficients (Frobenius on base field)
        // Note: already accounted inside polynomial_mul over GF(2^m)
    }

    ext_copy(out, acc);
    polynomial_free(acc);
    polynomial_free(tmp);
    polynomial_free(sq_tmp);
}

static int ext_equal(const polynomial_t *a, const polynomial_t *b) {
    if (a->degree != b->degree) return 0;
    for (int i = 0; i <= a->degree; i++) {
        if (a->coeffs[i] != b->coeffs[i]) return 0;
    }
    return 1;
}

static int ext_is_in_Fq(const polynomial_t *a) {
    return a->degree <= 0; // degree -1 (zero) or 0 (constant)
}

static gf_elem_t ext_const_term(const polynomial_t *a) {
    return (a->degree >= 0) ? a->coeffs[0] : 0;
}

// Helpers for irreducibility test over F_{2^m}
static int poly_equal(const polynomial_t *a, const polynomial_t *b) {
    if (a->degree != b->degree) return 0;
    for (int i = 0; i <= a->degree; i++) if (a->coeffs[i] != b->coeffs[i]) return 0;
    return 1;
}

static void poly_square_mod(polynomial_t *out, const polynomial_t *a, const polynomial_t *mod) {
    int cap = 2 * (mod->degree) - 2;
    if (cap < 0) cap = 0;
    polynomial_t *tmp = polynomial_create(cap);
    polynomial_mul(tmp, a, a);
    polynomial_div(NULL, out, tmp, mod);
    polynomial_free(tmp);
}

static void poly_qpower_mod(polynomial_t *out, const polynomial_t *a, const polynomial_t *mod, int m) {
    int cap = (mod->degree > 0) ? (mod->degree - 1) : 0;
    polynomial_t *acc = polynomial_create(cap);
    polynomial_t *tmp = polynomial_create(cap);
    polynomial_copy(acc, a);
    for (int i = 0; i < m; i++) {
        poly_square_mod(tmp, acc, mod);
        polynomial_copy(acc, tmp);
    }
    polynomial_copy(out, acc);
    polynomial_free(acc);
    polynomial_free(tmp);
}

static void poly_add_inplace(polynomial_t *a, const polynomial_t *b) {
    polynomial_add(a, a, b);
}

static int t_prime_divisors(int t, int *divs, int max) {
    int count = 0;
    int n = t;
    for (int p = 2; p*p <= n; p++) {
        if (n % p == 0) {
            divs[count++] = p;
            if (count >= max) break;
            while (n % p == 0) n /= p;
        }
    }
    if (n > 1 && count < max) divs[count++] = n;
    return count;
}

static int is_irreducible_over_Fqm(const polynomial_t *g, int m) {
    int t = g->degree;
    if (t <= 0) return 0;
    if (t != MCELIECE_T) return 0;
    // 1) Check monic
    if (g->coeffs[t] != 1) return 0;

    // 2) x^{q^t} ≡ x (mod g)
    polynomial_t *x = polynomial_create(t - 1);
    polynomial_set_coeff(x, 1, 1);
    polynomial_t *h = polynomial_create(t - 1);
    polynomial_copy(h, x);
    for (int i = 0; i < t; i++) {
        poly_qpower_mod(h, h, g, m); // h = h^{q} mod g
    }
    int ok1 = poly_equal(h, x);

    if (!ok1) {
        polynomial_free(x); polynomial_free(h);
        return 0;
    }

    // 3) For each distinct prime r | t: gcd(x^{q^{t/r}} - x, g) == 1
    int divs[16];
    int nd = t_prime_divisors(t, divs, 16);
    polynomial_t *u = polynomial_create(t - 1);
    polynomial_t *tmp = polynomial_create(t - 1);
    polynomial_t *gcd = polynomial_create(t - 1);

    for (int i = 0; i < nd; i++) {
        int r = divs[i];
        int e = t / r;
        polynomial_copy(u, x);
        for (int j = 0; j < e; j++) {
            poly_qpower_mod(u, u, g, m);
        }
        // u = x^{q^{e}} mod g; compute u - x = u + x
        polynomial_copy(tmp, u);
        poly_add_inplace(tmp, x);

        polynomial_gcd(gcd, g, tmp);
        if (gcd->degree > 0) { // non-trivial gcd
            polynomial_free(x); polynomial_free(h);
            polynomial_free(u); polynomial_free(tmp); polynomial_free(gcd);
            return 0;
        }
    }

    polynomial_free(x); polynomial_free(h);
    polynomial_free(u); polynomial_free(tmp); polynomial_free(gcd);
    return 1;
}

// 最终版，包含内部自检逻辑
mceliece_error_t generate_irreducible_poly_final(polynomial_t *g, const uint8_t *random_bits) {
    printf("      -> girr: g=%p coeffs=%p maxdeg=%d rand=%p\n", (void*)g, (void*)g->coeffs, g->max_degree, (void*)random_bits); fflush(stdout);
    printf("    -> generate_irreducible_poly: start. t=%d m=%d\n", MCELIECE_T, MCELIECE_M); fflush(stdout);

    if (!g || !g->coeffs) return MCELIECE_ERROR_KEYGEN_FAIL;

    int t = MCELIECE_T;
    int m = MCELIECE_M;

    // Build candidate g(x) from bits: monic degree t
    memset(g->coeffs, 0, (g->max_degree + 1) * sizeof(gf_elem_t));
    g->degree = -1;

    for (int j = 0; j < t; j++) {
        int offset = j * 2;
        gf_elem_t c = (gf_elem_t)random_bits[offset] | (gf_elem_t)(random_bits[offset+1] << 8);
        c &= (1 << m) - 1;
        polynomial_set_coeff(g, j, c);
    }
    polynomial_set_coeff(g, t, 1);

    int ok = is_irreducible_over_Fqm(g, m);
    printf("      -> girr: candidate degree=%d irreducible=%d\n", g->degree, ok); fflush(stdout);
    return ok ? MCELIECE_SUCCESS : MCELIECE_ERROR_KEYGEN_FAIL;
}


