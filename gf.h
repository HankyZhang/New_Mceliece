#ifndef GF_H
#define GF_H

#include <stdint.h>
#include "mceliece_types.h"  // 假设里面定义了 gf_elem_t 和 polynomial_t 等类型

#ifdef __cplusplus
extern "C" {
#endif

    // GF(2^13)加法（异或）
    gf_elem_t gf_add(gf_elem_t a, gf_elem_t b);

    // GF(2^13)乘法
    gf_elem_t gf_mul(gf_elem_t a, gf_elem_t b);

    // GF(2^13)求逆
    gf_elem_t gf_inv(gf_elem_t a);

    // GF(2^13)除法
    gf_elem_t gf_div(gf_elem_t a, gf_elem_t b);

    // GF(2^13)幂运算
    gf_elem_t gf_pow(gf_elem_t base, int exp);

    // 从比特向量表示转换为GF元素
    gf_elem_t bits_to_gf(const uint8_t *bits, int start_bit);
    extern gf_elem_t *gf_log;
    extern gf_elem_t *gf_antilog;

    void gf_init(void);
    // 从GF元素转换为比特向量表示
    void gf_to_bits(gf_elem_t elem, uint8_t *bits, int start_bit);

    // 多项式在GF中的求值
    gf_elem_t polynomial_eval(const polynomial_t *poly, gf_elem_t x);

    // 设置多项式系数并更新次数
    void polynomial_set_coeff(polynomial_t *poly, int degree, gf_elem_t coeff);

    // 多项式复制
    void polynomial_copy(polynomial_t *dst, const polynomial_t *src);

    // 检查多项式是否为零
    int polynomial_is_zero(const polynomial_t *poly);

    // 多项式加法（GF上就是异或）
    void polynomial_add(polynomial_t *result, const polynomial_t *a, const polynomial_t *b);
    // 多项式乘法: result(x) = a(x) * b(x)
    void polynomial_mul(polynomial_t *result, const polynomial_t *a, const polynomial_t *b);

    // 多项式除法: q(x) = a(x) / b(x), r(x) = a(x) mod b(x)
    void polynomial_div(polynomial_t *q, polynomial_t *r, const polynomial_t *a, const polynomial_t *b);

    // 多项式最大公约数: result(x) = gcd(a(x), b(x))
    void polynomial_gcd(polynomial_t *result, const polynomial_t *a, const polynomial_t *b);

    // 多项式 x^exp mod f(x)
    void polynomial_pow_mod(polynomial_t *result, unsigned long exp, const polynomial_t *f);


    // 函数原型
    matrix_fq_t* matrix_fq_create(int rows, int cols);
    void matrix_fq_free(matrix_fq_t *mat);
    gf_elem_t* solve_linear_system(const matrix_fq_t *A, const gf_elem_t *b);

    void polynomial_eval_at_poly(polynomial_t *result, const polynomial_t *g, const polynomial_t *beta, const polynomial_t *F);
#ifdef __cplusplus
}
#endif

#endif // GF_H
