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

// 最终版，包含内部自检逻辑
mceliece_error_t generate_irreducible_poly_final(polynomial_t *g, const uint8_t *random_bits) {
    int t = MCELIECE_T;
    int m = MCELIECE_M;
    mceliece_error_t ret = MCELIECE_SUCCESS;

    // --- 准备工作 ---
    polynomial_t *F_y = polynomial_create(t);
    polynomial_t *beta_poly = polynomial_create(t - 1);
    matrix_fq_t *A = matrix_fq_create(t, t);
    gf_elem_t *b = malloc(t * sizeof(gf_elem_t));
    polynomial_t *last_power = polynomial_create(t - 1);
    polynomial_t *temp_power = polynomial_create(2 * t - 2);
    polynomial_t *temp_mod = polynomial_create(t - 1);
    gf_elem_t *solution = NULL;

    if (!F_y || !beta_poly || !A || !b || !last_power || !temp_power || !temp_mod) {
        ret = MCELIECE_ERROR_MEMORY;
        goto cleanup;
    }

    get_field_poly_F(F_y);

    // --- 1. 构造 beta 多项式 ---
    for (int j = 0; j < t; j++) {
        int offset = j * 2;
        gf_elem_t beta_j = (gf_elem_t)random_bits[offset] | (gf_elem_t)(random_bits[offset+1] << 8);
        beta_j &= (1 << m) - 1;
        polynomial_set_coeff(beta_poly, j, beta_j);
    }

    // --- 2. 构造矩阵 A 和向量 b ---
    polynomial_set_coeff(last_power, 0, 1); // beta^0 = 1
    for (int i = 0; i < t; i++) {
         A->data[i * t + 0] = (i == 0) ? 1 : 0;
    }

    for (int k = 1; k <= t; k++) {
        polynomial_mul(temp_power, last_power, beta_poly);
        polynomial_div(NULL, temp_mod, temp_power, F_y);

        if (k < t) {
            for (int i = 0; i < t; i++) {
                A->data[i * t + k] = (i <= temp_mod->degree) ? temp_mod->coeffs[i] : 0;
            }
        } else {
            for (int i = 0; i < t; i++) {
                b[i] = (i <= temp_mod->degree) ? temp_mod->coeffs[i] : 0;
            }
        }
        polynomial_copy(last_power, temp_mod);
    }

    // --- 3. 求解线性方程组 ---
    solution = solve_linear_system(A, b);

    if (!solution) {
        ret = MCELIECE_ERROR_KEYGEN_FAIL;
        goto cleanup;
    }

    // --- 4. 构造多项式 g ---
    for(int i = 0; i < t; i++) {
        polynomial_set_coeff(g, i, solution[i]);
    }
    polynomial_set_coeff(g, t, 1);

    if (g->degree != t) {
        ret = MCELIECE_ERROR_KEYGEN_FAIL;
        goto cleanup;
    }

    // =============================================================
    //   5. 核心自检：验证我们计算出的 g 是否满足 g(beta) = 0
    // =============================================================
    printf("    -> generate_irreducible_poly: Verifying g(beta) == 0...\n");
    polynomial_t *check_result = polynomial_create(t - 1);
    if (!check_result) {
        ret = MCELIECE_ERROR_MEMORY;
        goto cleanup;
    }

    polynomial_eval_at_poly(check_result, g, beta_poly, F_y);

    if (!polynomial_is_zero(check_result)) {
        printf("    -> FATAL INTERNAL ERROR! Calculated g does not satisfy g(beta) = 0.\n");
        printf("    -> The math in polynomial or linear algebra functions is likely wrong.\n");
        ret = MCELIECE_ERROR_KEYGEN_FAIL;
    } else {
        printf("    -> Verification successful: g(beta) is indeed 0.\n");
    }
    polynomial_free(check_result);


cleanup:
    // --- 统一清理所有分配的内存 ---
    polynomial_free(F_y);
    polynomial_free(beta_poly);
    matrix_fq_free(A);
    free(b);
    polynomial_free(last_power);
    polynomial_free(temp_power);
    polynomial_free(temp_mod);
    free(solution);

    return ret;
}


