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
    // According to Classic McEliece specification, for m=13 the irreducible polynomial is:
    // x^13 + x^4 + x^3 + x + 1
    // In binary: 10000000011011 = 0x201B
    // For reduction, we use the polynomial without the leading bit: 0x001B
    const gf_elem_t reducing_poly = 0x001B;
    gf_elem_t r = 0;
    gf_elem_t temp_a = a;  // 使用临时变量，不修改原始参数

    // 我们将循环 m 次 (m=13)
    for (int i = 0; i < MCELIECE_M; i++) {
        // 如果 b 的当前最低位是 1
        if (b & 1) {
            r ^= temp_a;
        }

        // b 右移一位，准备处理下一位
        b >>= 1;

        // temp_a 左移一位 (相当于 temp_a = temp_a * x)
        // 检查最高位 (x^12) 是否为 1
        if (temp_a & (1 << (MCELIECE_M - 1))) {
            // 如果是，左移后会溢出，需要进行模约简
            // 1. 先左移
            temp_a <<= 1;
            // 2. 然后与约简多项式异或
            temp_a ^= reducing_poly;
        } else {
            // 如果不是，直接左移即可
            temp_a <<= 1;
        }
        // Ensure temp_a stays within the field bounds (13 bits)
        temp_a &= ((1 << MCELIECE_M) - 1);
    }

    // Ensure result stays within the field bounds (13 bits)
    return r & ((1 << MCELIECE_M) - 1);
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

    // Initialize GF(2^13) lookup tables using generator element 3
    const gf_elem_t generator = 3;
    gf_elem_t p = 1;
    int i;

    assert(MCELIECE_M == 13);
    assert(MCELIECE_Q == 8192);

    for (i = 0; i < MCELIECE_Q - 1; i++) {
        gf_antilog[i] = p;
        gf_log[p] = (gf_elem_t)i;
        gf_elem_t old_p = p;
        p = gf_mul_for_init(p, generator);
        
        // Check for infinite loop condition
        if (i > 0 && p == 1) {
            break;
        }
        
        // Prevent potential infinite loops
        if (i > 0 && p == old_p) {
            fprintf(stderr, "ERROR: GF table generation stuck at iteration %d\n", i);
            exit(1);
        }
    }

    gf_log[0] = 0;
}
gf_elem_t gf_add(gf_elem_t a, gf_elem_t b) {
    return a ^ b;
}

// ----> 您为外部使用而定义的、高效的查表版本保持不变 <----

// Optimized GF(2^13) multiplication using lookup tables
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

// Efficient polynomial evaluation using Horner's method
gf_elem_t polynomial_eval(const polynomial_t *poly, gf_elem_t x) {
    if (poly->degree < 0) {
        return 0; // Zero polynomial
    }

    // Use Horner's method: start with highest degree coefficient
    gf_elem_t result = poly->coeffs[poly->degree];

    // Iterate down to constant term
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














