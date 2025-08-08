#include "mceliece_types.h"
#include "gf.h"
#include <stdio.h>
#include <assert.h>

gf_elem_t *gf_log = NULL;
gf_elem_t *gf_antilog = NULL;


/*
 * 这是一个已知正确的、用于初始化的GF(2^m)乘法实现。
 * 它使用标准的“俄罗斯农夫乘法”(位移和异或)算法。
 */
static gf_elem_t gf_mul_for_init(gf_elem_t a, gf_elem_t b) {
    // Classic McEliece m=13 的不可约多项式是 x^13 + x^4 + x^3 + x + 1
    // 这对应于二进制的 0b1 0000 0000 11011 = 0x201B
    // 我们需要的是去掉最高位的多项式，用于约简：0x001B
    const gf_elem_t reducing_poly = 0x001B;
    gf_elem_t r = 0;

    // 我们将循环 m 次 (m=13)
    for (int i = 0; i < MCELIECE_M; i++) {
        // 如果 b 的当前最低位是 1
        if (b & 1) {
            r ^= a;
        }

        // b 右移一位，准备处理下一位
        b >>= 1;

        // a 左移一位 (相当于 a = a * x)
        // 检查最高位 (x^12) 是否为 1
        if (a & (1 << (MCELIECE_M - 1))) {
            // 如果是，左移后会溢出，需要进行模约简
            // 1. 先左移
            a <<= 1;
            // 2. 然后与约简多项式异或
            a ^= reducing_poly;
        } else {
            // 如果不是，直接左移即可
            a <<= 1;
        }
    }

    return r;
}

void gf_init(void) {
    // ----> 新增：在函数开头进行内存分配 <----
    if (gf_log == NULL) { // 确保只分配一次
        gf_log = malloc(MCELIECE_Q * sizeof(gf_elem_t));
        if (gf_log == NULL) {
            fprintf(stderr, "Failed to allocate memory for gf_log\n");
            exit(1);
        }
    }
    if (gf_antilog == NULL) {
        gf_antilog = malloc(MCELIECE_Q * sizeof(gf_elem_t));
        if (gf_antilog == NULL) {
            fprintf(stderr, "Failed to allocate memory for gf_antilog\n");
            exit(1);
        }
    }

    // ----> 后续的初始化逻辑保持完全不变 <----
    const gf_elem_t generator = 3;
    gf_elem_t p = 1;
    int i;

    assert(MCELIECE_M == 13);
    assert(MCELIECE_Q == 8192);

    for (i = 0; i < MCELIECE_Q - 1; i++) {
        gf_antilog[i] = p;
        gf_log[p] = (gf_elem_t)i;
        p = gf_mul_for_init(p, generator);
    }

    gf_log[0] = 0;

    printf("GF(2^13) tables initialized successfully on the heap.\n");
    fflush(stdout);
}
gf_elem_t gf_add(gf_elem_t a, gf_elem_t b) {
    return a ^ b;
}

// ----> 您为外部使用而定义的、高效的查表版本保持不变 <----

// 新的、快速且正确的 GF(2^13) 乘法
gf_elem_t gf_mul(gf_elem_t a, gf_elem_t b) {
    if (a == 0 || b == 0) {
        return 0;
    }
    int log_a = gf_log[a];
    int log_b = gf_log[b];
    int sum_log = log_a + log_b;
    if (sum_log >= MCELIECE_Q - 1) {
        sum_log -= (MCELIECE_Q - 1);
    }
    return gf_antilog[sum_log];
}

// 新的、快速且正确的 GF(2^13) 求逆
gf_elem_t gf_inv(gf_elem_t a) {
    if (a == 0) {
        return 0;
    }
    // 确保 log[a] 在范围内，避免负数索引
    if (gf_log[a] == 0 && a != 1) return 0; // 处理未初始化或无效情况
    return gf_antilog[(MCELIECE_Q - 1) - gf_log[a]];
}

// GF(2^13)除法
gf_elem_t gf_div(gf_elem_t a, gf_elem_t b) {
    if (b == 0) return 0;  // 除零错误
    return gf_mul(a, gf_inv(b));
}

// GF(2^13)幂运算
gf_elem_t gf_pow(gf_elem_t base, int exp) {
    if (exp == 0) return 1;
    if (base == 0) return 0;

    gf_elem_t result = 1;
    base = base & ((1 << MCELIECE_M) - 1);

    while (exp > 0) {
        if (exp & 1) {
            result = gf_mul(result, base);
        }
        base = gf_mul(base, base);
        exp >>= 1;
    }

    return result;
}

// 从比特向量表示转换为GF元素
gf_elem_t bits_to_gf(const uint8_t *bits, int start_bit) {
    gf_elem_t result = 0;

    for (int i = 0; i < MCELIECE_M; i++) {
        int byte_idx = (start_bit + i) / 8;
        int bit_idx = (start_bit + i) % 8;

        if (bits[byte_idx] & (1 << bit_idx)) {
            result |= (1 << i);
        }
    }

    return result;
}

// 从GF元素转换为比特向量表示
void gf_to_bits(gf_elem_t elem, uint8_t *bits, int start_bit) {
    for (int i = 0; i < MCELIECE_M; i++) {
        int byte_idx = (start_bit + i) / 8;
        int bit_idx = (start_bit + i) % 8;

        if (elem & (1 << i)) {
            bits[byte_idx] |= (1 << bit_idx);
        } else {
            bits[byte_idx] &= ~(1 << bit_idx);
        }
    }
}

// 多项式在GF中的求值
gf_elem_t polynomial_eval(const polynomial_t *poly, gf_elem_t x) {
    if (poly->degree < 0) {
        return 0; // 零多项式
    }

    // 从最高次项系数开始
    gf_elem_t result = poly->coeffs[poly->degree];

    // 向下迭代到常数项
    for (int i = poly->degree - 1; i >= 0; i--) {
        result = gf_mul(result, x);
        result = gf_add(result, poly->coeffs[i]);
    }

    return result;
}


// 设置多项式系数并更新次数
void polynomial_set_coeff(polynomial_t *poly, int degree, gf_elem_t coeff) {
    if (degree > poly->max_degree) return;

    poly->coeffs[degree] = coeff;

    // 更新多项式次数
    if (coeff != 0 && degree > poly->degree) {
        poly->degree = degree;
    } else if (coeff == 0 && degree == poly->degree) {
        // 如果清零了最高次项，需要重新计算次数
        int new_degree = -1;
        for (int i = poly->max_degree; i >= 0; i--) {
            if (poly->coeffs[i] != 0) {
                new_degree = i;
                break;
            }
        }
        poly->degree = new_degree;
    }
}

void polynomial_copy(polynomial_t *dst, const polynomial_t *src) {
    if (!dst || !src) return;

    // 清零目标多项式
    memset(dst->coeffs, 0, (dst->max_degree + 1) * sizeof(gf_elem_t));

    // 确定要复制的项数，不能超过目标的最大容量
    int terms_to_copy = (src->degree < dst->max_degree) ? src->degree : dst->max_degree;

    // 如果源多项式是零多项式
    if (src->degree < 0) {
        dst->degree = -1;
        return;
    }

    // 复制系数
    for (int i = 0; i <= terms_to_copy; i++) {
        dst->coeffs[i] = src->coeffs[i];
    }

    // 设置次数
    dst->degree = terms_to_copy;

    // 如果发生了截断，需要重新检查最高次项是否为0
    while(dst->degree >= 0 && dst->coeffs[dst->degree] == 0) {
        dst->degree--;
    }
}

// 检查多项式是否为零
int polynomial_is_zero(const polynomial_t *poly) {
    return poly->degree < 0;
}

// 多项式加法（GF上就是异或）
void polynomial_add(polynomial_t *result, const polynomial_t *a, const polynomial_t *b) {
    int max_deg = (a->degree > b->degree) ? a->degree : b->degree;

    if (max_deg > result->max_degree) return;

    // 清零结果
    memset(result->coeffs, 0, (result->max_degree + 1) * sizeof(gf_elem_t));

    for (int i = 0; i <= max_deg; i++) {
        gf_elem_t coeff_a = (i <= a->degree) ? a->coeffs[i] : 0;
        gf_elem_t coeff_b = (i <= b->degree) ? b->coeffs[i] : 0;
        result->coeffs[i] = gf_add(coeff_a, coeff_b);
    }

    // 重新计算次数
    result->degree = -1;
    for (int i = max_deg; i >= 0; i--) {
        if (result->coeffs[i] != 0) {
            result->degree = i;
            break;
        }
    }
}


// --- 实现新增的高级多项式运算 ---

void polynomial_mul(polynomial_t *result, const polynomial_t *a, const polynomial_t *b) {
    int deg_res = a->degree + b->degree;
    assert(deg_res <= result->max_degree);

    // 清零结果
    memset(result->coeffs, 0, (result->max_degree + 1) * sizeof(gf_elem_t));
    result->degree = -1;

    for (int i = 0; i <= a->degree; i++) {
        for (int j = 0; j <= b->degree; j++) {
            gf_elem_t term = gf_mul(a->coeffs[i], b->coeffs[j]);
            result->coeffs[i + j] = gf_add(result->coeffs[i + j], term);
        }
    }

    // 更新次数
    for (int i = deg_res; i >= 0; i--) {
        if (result->coeffs[i] != 0) {
            result->degree = i;
            return;
        }
    }
}

void polynomial_div(polynomial_t *q, polynomial_t *r, const polynomial_t *a, const polynomial_t *b) {
    assert(b->degree >= 0); // 不能除以零多项式

    polynomial_copy(r, a); // 余数 r 初始化为 a

    if (q) {
        memset(q->coeffs, 0, (q->max_degree + 1) * sizeof(gf_elem_t));
        q->degree = -1;
    }

    // 如果被除数次数小于除数，商为0，余数为被除数本身
    if (r->degree < b->degree) {
        if (q) q->degree = -1;
        return;
    }

    int deg_b = b->degree;
    gf_elem_t lead_b_inv = gf_inv(b->coeffs[deg_b]);

    // 长除法核心循环
    for (int i = r->degree; i >= deg_b; i--) {
        gf_elem_t coeff = gf_mul(r->coeffs[i], lead_b_inv);

        if (coeff != 0) {
            if (q) {
                polynomial_set_coeff(q, i - deg_b, coeff);
            }
            // 从 r 中减去 (实际上是加上) coeff * x^(i-deg_b) * b
            for (int j = 0; j <= deg_b; j++) {
                gf_elem_t term = gf_mul(coeff, b->coeffs[j]);
                int r_idx = i - deg_b + j;
                if (r_idx <= r->max_degree) {
                    r->coeffs[r_idx] = gf_add(r->coeffs[r_idx], term);
                }
            }
        }
    }

    // 更新余数 r 的真实次数
    int new_r_degree = -1;
    for (int i = deg_b - 1; i >= 0; i--) {
        if (i <= r->max_degree && r->coeffs[i] != 0) {
            new_r_degree = i;
            break;
        }
    }
    r->degree = new_r_degree;
}

void polynomial_gcd(polynomial_t *result, const polynomial_t *a, const polynomial_t *b) {
    polynomial_t *tmp_a = polynomial_create(a->max_degree);
    polynomial_t *tmp_b = polynomial_create(b->max_degree);
    polynomial_t *r = polynomial_create(b->max_degree);

    polynomial_copy(tmp_a, a);
    polynomial_copy(tmp_b, b);

    while (!polynomial_is_zero(tmp_b)) {
        polynomial_div(NULL, r, tmp_a, tmp_b);
        polynomial_copy(tmp_a, tmp_b);
        polynomial_copy(tmp_b, r);
    }

    polynomial_copy(result, tmp_a);

    polynomial_free(tmp_a);
    polynomial_free(tmp_b);
    polynomial_free(r);
}

void polynomial_pow_mod(polynomial_t *result, unsigned long exp, const polynomial_t *f) {
    polynomial_t *s = polynomial_create(f->degree * 2);
    polynomial_set_coeff(s, 1, 1); // s = x

    polynomial_t *tmp = polynomial_create(f->degree * 2);
    polynomial_set_coeff(result, 0, 1); // result = 1

    while (exp > 0) {
        if (exp & 1) {
            polynomial_mul(tmp, result, s);
            polynomial_div(NULL, result, tmp, f);
        }
        polynomial_mul(tmp, s, s);
        polynomial_div(NULL, s, tmp, f);
        exp >>= 1;
    }

    polynomial_free(s);
    polynomial_free(tmp);
}


// 创建 F_q 矩阵
matrix_fq_t* matrix_fq_create(int rows, int cols) {
    matrix_fq_t *mat = malloc(sizeof(matrix_fq_t));
    if (!mat) return NULL;
    mat->rows = rows;
    mat->cols = cols;
    mat->data = calloc((size_t)rows * cols, sizeof(gf_elem_t));
    if (!mat->data) {
        free(mat);
        return NULL;
    }
    return mat;
}

// 释放 F_q 矩阵
void matrix_fq_free(matrix_fq_t *mat) {
    if (mat) {
        free(mat->data);
        free(mat);
    }
}


// =======================================================
//   核心：求解线性方程组 A * x = b
//   返回一个包含解的数组，如果无解或解不唯一则返回 NULL。
//   调用者负责释放返回的数组。
// =======================================================
// =======================================================
//   最终修正版：solve_linear_system
// =======================================================
gf_elem_t* solve_linear_system(const matrix_fq_t *A, const gf_elem_t *b) {
    int n = A->rows;
    if (n != A->cols) return NULL; // 只支持方阵

    // 1. 创建增广矩阵 [A | b]
    matrix_fq_t *aug = matrix_fq_create(n, n + 1);
    if (!aug) return NULL;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug->data[i * (n + 1) + j] = A->data[i * n + j];
        }
        aug->data[i * (n + 1) + n] = b[i];
    }

    // 2. 高斯-若尔当消元
    for (int i = 0; i < n; i++) {
        // --- 寻找主元 (在第 i 列，从第 i 行开始) ---
        int pivot_row = i;
        while (pivot_row < n && aug->data[pivot_row * (n + 1) + i] == 0) {
            pivot_row++;
        }

        if (pivot_row == n) { // 在当前列找不到主元，矩阵奇异
            matrix_fq_free(aug);
            return NULL; // 无唯一解
        }

        // --- 交换行，把主元行换到第 i 行 ---
        if (pivot_row != i) {
            for (int k = 0; k < n + 1; k++) {
                gf_elem_t temp = aug->data[i * (n + 1) + k];
                aug->data[i * (n + 1) + k] = aug->data[pivot_row * (n + 1) + k];
                aug->data[pivot_row * (n + 1) + k] = temp;
            }
        }

        // --- 将主元归一化为 1 (通过乘以其逆元) ---
        // 对整行 (从第 i 列开始即可，因为左边都是0) 进行操作
        gf_elem_t pivot_val = aug->data[i * (n + 1) + i];
        gf_elem_t pivot_inv = gf_inv(pivot_val);
        for (int k = i; k < n + 1; k++) {
            aug->data[i * (n + 1) + k] = gf_mul(aug->data[i * (n + 1) + k], pivot_inv);
        }

        // --- 将其他所有行的第 i 列都消为 0 ---
        for (int j = 0; j < n; j++) {
            if (i != j) {
                gf_elem_t factor = aug->data[j * (n + 1) + i];
                if (factor != 0) {
                    // 对整行 (从第 i 列开始即可) 进行操作
                    for (int k = i; k < n + 1; k++) {
                        gf_elem_t term = gf_mul(factor, aug->data[i * (n + 1) + k]);
                        aug->data[j * (n + 1) + k] = gf_add(aug->data[j * (n + 1) + k], term);
                    }
                }
            }
        }
    }

    // 3. 提取解
    gf_elem_t *solution = malloc(n * sizeof(gf_elem_t));
    if (!solution) {
        matrix_fq_free(aug);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        solution[i] = aug->data[i * (n + 1) + n];
    }

    matrix_fq_free(aug);
    return solution;
}


void polynomial_eval_at_poly(polynomial_t *result, const polynomial_t *g, const polynomial_t *beta, const polynomial_t *F) {
    if (g->degree < 0) {
        // 如果 g 是零多项式，结果也是零
        polynomial_set_coeff(result, 0, 0);
        return;
    }

    int t = F->degree;
    polynomial_t *temp_mul = polynomial_create(2 * t);
    polynomial_t *temp_mod = polynomial_create(t - 1);

    // 使用霍纳法则：从最高次项开始
    // result 初始化为 g 的最高次项系数
    polynomial_set_coeff(result, 0, g->coeffs[g->degree]);

    for (int i = g->degree - 1; i >= 0; i--) {
        // result = (result * beta) + g_i

        // 1. 计算 result * beta
        polynomial_mul(temp_mul, result, beta);

        // 2. 取模 F(y)
        polynomial_div(NULL, temp_mod, temp_mul, F);

        // 3. 加上 g 的下一个系数 g_i
        temp_mod->coeffs[0] = gf_add(temp_mod->coeffs[0], g->coeffs[i]);
        polynomial_set_coeff(temp_mod, 0, temp_mod->coeffs[0]); // 更新次数

        // 4. 将结果复制回 result，用于下一次迭代
        polynomial_copy(result, temp_mod);
    }

    // 清理
    polynomial_free(temp_mul);
    polynomial_free(temp_mod);
}