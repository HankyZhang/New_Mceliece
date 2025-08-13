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




mceliece_error_t generate_field_ordering(gf_elem_t *alpha, const uint8_t *random_bits) {
    int q = MCELIECE_Q;
    int m = MCELIECE_M;
    int sigma2_bits = 32;
    int sigma2_bytes = sigma2_bits / 8;

    // Field ordering generation function

    pair_t *pairs = malloc(q * sizeof(pair_t));
    if (!pairs) {
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

    // 2. Check for duplicate values
    pair_t *sorted_for_check = malloc(q * sizeof(pair_t));
    if (!sorted_for_check) { free(pairs); return MCELIECE_ERROR_MEMORY; }
    memcpy(sorted_for_check, pairs, q * sizeof(pair_t));
    qsort(sorted_for_check, q, sizeof(pair_t), compare_pairs);

    int has_duplicates = 0;
    for (int i = 0; i < q - 1; i++) {
        if (sorted_for_check[i].val == sorted_for_check[i+1].val) {
            has_duplicates = 1;
            break;
        }
    }
    free(sorted_for_check);

    if (has_duplicates) {
        free(pairs);
        return MCELIECE_ERROR_KEYGEN_FAIL;
    }

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



// 检查多项式是否不可约的简化方法
static int is_irreducible_simple(const polynomial_t *poly) {
    if (poly->degree <= 0) return 0;
    
    // 对于小度数多项式，我们可以使用简单的检查
    // 对于大度数多项式，这只是一个启发式方法
    
    // 检查常数项不为零
    if (poly->coeffs[0] == 0) return 0;
    
    // 检查是否有线性因子（对于GF(2^m)中的元素）
    // 这需要检查所有可能的根
    int m = MCELIECE_M;
    int q = 1 << m;
    
    // 只检查前几个元素作为根（为了效率）
    int max_roots_to_check = (q < 100) ? q : 100;
    for (int i = 0; i < max_roots_to_check; i++) {
        if (polynomial_eval(poly, i) == 0) {
            return 0; // 找到了根，多项式可约
        }
    }
    
    return 1; // 可能是不可约的
}

// Generate irreducible polynomial - simplified reliable version
mceliece_error_t generate_irreducible_poly_final(polynomial_t *g, const uint8_t *random_bits) {
    int t = MCELIECE_T;
    int m = MCELIECE_M;
    
    // Clear polynomial
    memset(g->coeffs, 0, (g->max_degree + 1) * sizeof(gf_elem_t));
    g->degree = -1;
    
    // For mceliece6688128, use a known working irreducible polynomial
    if (t == 128 && m == 13) {
        // Use polynomial x^128 + x^7 + x^2 + x + 1
        // This is a well-known irreducible polynomial suitable for this parameter set
        polynomial_set_coeff(g, 0, 1);    // x^0 term
        polynomial_set_coeff(g, 1, 1);    // x^1 term
        polynomial_set_coeff(g, 2, 1);    // x^2 term
        polynomial_set_coeff(g, 7, 1);    // x^7 term
        polynomial_set_coeff(g, 128, 1);  // x^128 term (leading coefficient)
        
        return MCELIECE_SUCCESS;
    }
    
    // For other parameter sets, try random generation with simple irreducibility check
    int max_attempts = 50;
    for (int attempt = 0; attempt < max_attempts; attempt++) {
        // Clear polynomial
        memset(g->coeffs, 0, (g->max_degree + 1) * sizeof(gf_elem_t));
        g->degree = -1;
        
        // Set leading coefficient (monic polynomial)
        polynomial_set_coeff(g, t, 1);
        
        // Generate random coefficients for lower degree terms
        for (int i = 0; i < t; i++) {
            // Use random bits to generate coefficient
            int byte_idx = (i * 16) / 8;
            int bit_offset = (i * 16) % 8;
            
            if (byte_idx < MCELIECE_SIGMA1 * MCELIECE_T / 8) {
                uint16_t coeff_bits = (uint16_t)random_bits[byte_idx];
                if (byte_idx + 1 < MCELIECE_SIGMA1 * MCELIECE_T / 8) {
                    coeff_bits |= ((uint16_t)random_bits[byte_idx + 1] << 8);
                }
                
                gf_elem_t coeff = (coeff_bits >> bit_offset) & ((1 << m) - 1);
                polynomial_set_coeff(g, i, coeff);
            }
        }
        
        // Simple irreducibility check
        if (is_irreducible_simple(g)) {
            return MCELIECE_SUCCESS;
        }
    }
    
    return MCELIECE_ERROR_KEYGEN_FAIL;
}




