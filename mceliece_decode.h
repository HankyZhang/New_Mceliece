#ifndef MCELIECE_DECODE_H
#define MCELIECE_DECODE_H

#include <stdint.h>
#include "mceliece_types.h"  // 你的类型定义文件

#ifdef __cplusplus
extern "C" {
#endif

    // Berlekamp-Massey算法：求解线性反馈移位寄存器
    mceliece_error_t berlekamp_massey(const gf_elem_t *syndrome,
                                     polynomial_t *sigma, polynomial_t *omega);

    // Chien搜索：寻找错误定位多项式的根
    mceliece_error_t chien_search(const polynomial_t *sigma, const gf_elem_t *alpha,
                                 int *error_positions, int *num_errors);

    // Forney算法：计算错误值（对于二元码，错误值总是1）
    mceliece_error_t forney_algorithm(const polynomial_t *sigma, const polynomial_t *omega,
                                     const gf_elem_t *alpha, const int *error_positions,
                                     int num_errors, gf_elem_t *error_values);

    // 完整的解码算法
    mceliece_error_t decode_goppa(const uint8_t *received, const polynomial_t *g,
                                 const gf_elem_t *alpha, uint8_t *error_vector,
                                 int *decode_success);

#ifdef __cplusplus
}
#endif

#endif // MCELIECE_DECODE_H
