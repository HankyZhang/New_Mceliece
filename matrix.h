#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>
#include "mceliece_types.h"
#include "gf.h"


// 矩阵创建与释放
matrix_t* matrix_create(int rows, int cols);
void matrix_free(matrix_t *mat);

// 设置与读取矩阵元素（按bit操作）
void matrix_set_bit(matrix_t *mat, int row, int col, int value);
int matrix_get_bit(const matrix_t *mat, int row, int col);

// 行列基本操作
void matrix_swap_rows(matrix_t *mat, int row1, int row2);
void matrix_swap_cols(matrix_t *mat, int col1, int col2);
void matrix_xor_rows(matrix_t *mat, int row_dst, int row_src);

// 高斯消元与系统形式检查
mceliece_error_t matrix_rref(matrix_t *mat, int *pivot_cols, int *rank);
int matrix_is_systematic(const matrix_t *mat);
int reduce_to_systematic_form(matrix_t *H);
mceliece_error_t matrix_to_systematic(matrix_t *mat, int *p, matrix_t *T_out);
// 打印与拷贝
void print_matrix(const char *label, const matrix_t *mat);
void matrix_copy(matrix_t *dst, const matrix_t *src);

// 向量操作
void matrix_vector_multiply(const matrix_t *mat, const uint8_t *vec, uint8_t *result);

// 生成校验矩阵 H = [I | T]
matrix_t* construct_parity_check_matrix(const matrix_t *T);

// MatGen 与编码
mceliece_error_t mat_gen(const polynomial_t *g, const gf_elem_t *alpha,
                        matrix_t *T_out);
void encode_vector(const uint8_t *error_vector, const matrix_t *T, uint8_t *ciphertext);

// syndrome 计算
void compute_syndrome(const uint8_t *received, const polynomial_t *g,
                      const gf_elem_t *alpha, gf_elem_t *syndrome);

#endif // MATRIX_H
