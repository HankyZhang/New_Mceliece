#include "mceliece_types.h"
#include "matrix.h"
#include "gf.h"
#include "vector.h"
#include "matrix.h"
#include "mceliece.h"
#include <stdio.h>
#include <stdlib.h> // For malloc and free
#include <string.h> // For memset

// 矩阵行交换
void matrix_swap_rows(matrix_t *mat, int row1, int row2) {
    if (row1 == row2 || row1 >= mat->rows || row2 >= mat->rows) return;
    
    for (int col = 0; col < mat->cols_bytes; col++) {
        uint8_t temp = mat->data[row1 * mat->cols_bytes + col];
        mat->data[row1 * mat->cols_bytes + col] = mat->data[row2 * mat->cols_bytes + col];
        mat->data[row2 * mat->cols_bytes + col] = temp;
    }
}
int reduce_to_systematic_form(matrix_t *H) {
    int mt = H->rows;
    int n = H->cols;
    int i, j;

    // Create column permutation array
    int *col_perm = malloc(n * sizeof(int));
    if (!col_perm) return -1;
    
    // Initialize column permutation
    for (i = 0; i < n; i++) {
        col_perm[i] = i;
    }

    // --- Forward elimination to form upper triangular matrix ---
    for (i = 0; i < mt; i++) {
        // 1. Find pivot (starting from column i, starting from row i)
        int pivot_row = -1;
        int pivot_col = -1;
        
        for (int col = i; col < n; col++) {
            for (int row = i; row < mt; row++) {
                if (matrix_get_bit(H, row, col) == 1) {
                    pivot_row = row;
                    pivot_col = col;
                    break;
                }
            }
            if (pivot_row != -1) break;
        }

        if (pivot_row == -1) {
            // No pivot found, matrix is singular
            free(col_perm);
            return -1;
        }

        // 2. Swap pivot row to current row i
        if (pivot_row != i) {
            matrix_swap_rows(H, i, pivot_row);
        }

        // 3. Swap pivot column to current column i
        if (pivot_col != i) {
            matrix_swap_cols(H, i, pivot_col);
            // Update column permutation
            int temp = col_perm[i];
            col_perm[i] = col_perm[pivot_col];
            col_perm[pivot_col] = temp;
        }

        // 4. Eliminate all elements below the pivot in column i
        for (j = i + 1; j < mt; j++) {
            if (matrix_get_bit(H, j, i) == 1) {
                matrix_xor_rows(H, j, i);
            }
        }
    }

    // --- Back elimination to form identity matrix ---
    for (i = mt - 1; i >= 0; i--) {
        for (j = 0; j < i; j++) {
            if (matrix_get_bit(H, j, i) == 1) {
                matrix_xor_rows(H, j, i);
            }
        }
    }

    free(col_perm);
    return 0; // Success
}

// 矩阵列交换
void matrix_swap_cols(matrix_t *mat, int col1, int col2) {
    if (col1 == col2 || col1 >= mat->cols || col2 >= mat->cols) return;
    
    for (int row = 0; row < mat->rows; row++) {
        int bit1 = matrix_get_bit(mat, row, col1);
        int bit2 = matrix_get_bit(mat, row, col2);
        matrix_set_bit(mat, row, col1, bit2);
        matrix_set_bit(mat, row, col2, bit1);
    }
}

// 矩阵行异或（row_dst = row_dst XOR row_src）
void matrix_xor_rows(matrix_t *mat, int row_dst, int row_src) {
    if (row_dst >= mat->rows || row_src >= mat->rows) return;
    
    for (int col = 0; col < mat->cols_bytes; col++) {
        mat->data[row_dst * mat->cols_bytes + col] ^= 
            mat->data[row_src * mat->cols_bytes + col];
    }
}

// 寻找指定列的第一个非零行（从start_row开始）
int matrix_find_pivot(const matrix_t *mat, int col, int start_row) {
    for (int row = start_row; row < mat->rows; row++) {
        if (matrix_get_bit(mat, row, col)) {
            return row;
        }
    }
    return -1;  // 未找到
}


// 检查矩阵是否为系统形式
int matrix_is_systematic(const matrix_t *mat) {
    if (mat->rows > mat->cols) return 0;
    
    // 检查前mat->rows列是否构成单位矩阵
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->rows; j++) {
            int expected = (i == j) ? 1 : 0;
            if (matrix_get_bit(mat, i, j) != expected) {
                return 0;
            }
        }
    }
    
    return 1;
}



mceliece_error_t mat_gen(const polynomial_t *g, const gf_elem_t *alpha,
                        matrix_t *T_out) { // 注意：p_out 参数被移除了
    if (!g || !alpha || !T_out) {
        return MCELIECE_ERROR_INVALID_PARAM;
    }

    int n = MCELIECE_N;
    int t = MCELIECE_T;
    int m = MCELIECE_M;
    int mt = m * t;
    int k = n - mt;

    // Create Goppa parity-check matrix H
    matrix_t *H = matrix_create(mt, n);
    if (!H) return MCELIECE_ERROR_MEMORY;

    // According to specification 1.2.7, construct Goppa code parity-check matrix
    // M[i,j] = alpha[j]^i / g(alpha[j]) for i=0,...,t-1 and j=0,...,n-1
    
    for (int i = 0; i < t; i++) {
        for (int j = 0; j < n; j++) {
            // Calculate alpha[j]^i using efficient exponentiation
            gf_elem_t alpha_power = gf_pow(alpha[j], i);
            
            // Calculate g(alpha[j])
            gf_elem_t g_alpha = polynomial_eval(g, alpha[j]);
            
            // Check if g(alpha[j]) is zero (would cause division by zero)
            if (g_alpha == 0) {
                matrix_free(H);
                return MCELIECE_ERROR_KEYGEN_FAIL;
            }
            
            // Calculate M[i,j] = alpha[j]^i / g(alpha[j])
            gf_elem_t M_ij = gf_div(alpha_power, g_alpha);
            
            // Expand GF(2^m) element into m binary bits
            for (int bit = 0; bit < m; bit++) {
                int bit_value = (M_ij >> bit) & 1;
                matrix_set_bit(H, i * m + bit, j, bit_value);
            }
        }
    }
    
    // Convert H to systematic form [I_mt | T]
    if (reduce_to_systematic_form(H) != 0) {
        matrix_free(H);
        return MCELIECE_ERROR_KEYGEN_FAIL; // Matrix is singular
    }

    // At this point H is in the form [I_mt | T]
    // Extract public key T from the right side of H
    for (int i = 0; i < mt; i++) {
        for (int j = 0; j < k; j++) {
            int bit = matrix_get_bit(H, i, mt + j);
            matrix_set_bit(T_out, i, j, bit);
        }
    }

    matrix_free(H);
    return MCELIECE_SUCCESS;
}

// Encode算法：C = He，其中H = (I_mt | T)
void encode_vector(const uint8_t *error_vector, const matrix_t *T, uint8_t *ciphertext) {
    if (!error_vector || !T || !ciphertext) return;
    
    int mt = MCELIECE_M * MCELIECE_T;
    int mt_bytes = (mt + 7) / 8;
    
    // 清零密文
    memset(ciphertext, 0, mt_bytes);
    
    // C = H * e，其中H = (I_mt | T)
    // 由于H的前mt列是单位矩阵，所以前mt位的C直接等于e的前mt位
    
    // 复制e的前mt位到C
    for (int i = 0; i < mt; i++) {
        // 1. 从 error_vector 获取比特值 (0 或 1)
        int bit = vector_get_bit(error_vector, i);

        // 2. 将获取到的比特值直接设置到 ciphertext 中
        vector_set_bit(ciphertext, i, bit);
    }
    
    // 计算T矩阵与e的后k位的乘积，并异或到C中
    for (int row = 0; row < mt; row++) {
        int sum = 0;
        for (int col = 0; col < T->cols; col++) {
            int e_bit = vector_get_bit(error_vector, mt + col);
            int T_bit = matrix_get_bit(T, row, col);
            sum ^= (e_bit & T_bit);
        }
        
        if (sum) {
            // 直接异或当前位
            int current_bit = vector_get_bit(ciphertext, row);
            vector_set_bit(ciphertext, row, current_bit ^ 1);
        }
    }
}

// 矩阵向量乘法：result = mat * vec
void matrix_vector_multiply(const matrix_t *mat, const uint8_t *vec, uint8_t *result) {
    if (!mat || !vec || !result) return;
    
    int result_bytes = (mat->rows + 7) / 8;
    memset(result, 0, result_bytes);
    
    for (int row = 0; row < mat->rows; row++) {
        int sum = 0;

        for (int col = 0; col < mat->cols; col++) {
            int mat_bit = matrix_get_bit(mat, row, col);
            int vec_bit = vector_get_bit(vec, col);
            sum ^= (mat_bit & vec_bit);
        }

        // 去掉 if 语句，直接用 sum 作为 value 参数
        vector_set_bit(result, row, sum);
    }
}

// 构造完整的校验矩阵H = (I_mt | T)
matrix_t* construct_parity_check_matrix(const matrix_t *T) {
    if (!T) return NULL;
    
    int mt = T->rows;
    int total_cols = mt + T->cols;
    
    matrix_t *H = matrix_create(mt, total_cols);
    if (!H) return NULL;
    
    // 设置单位矩阵部分 I_mt
    for (int i = 0; i < mt; i++) {
        matrix_set_bit(H, i, i, 1);
    }
    
    // 复制T矩阵到右侧
    for (int row = 0; row < mt; row++) {
        for (int col = 0; col < T->cols; col++) {
            int bit = matrix_get_bit(T, row, col);
            matrix_set_bit(H, row, mt + col, bit);
        }
    }
    
    return H;
}

// Calculate syndrome according to Classic McEliece specification
void compute_syndrome(const uint8_t *received, const polynomial_t *g, 
                     const gf_elem_t *alpha, gf_elem_t *syndrome) {
    if (!received || !g || !alpha || !syndrome) return;
    
    // According to section 3.3.1: s_j = Σ_{i∈I} α_i^j / g(α_i)^2
    // where I is the set of error positions
    
    for (int j = 0; j < 2 * MCELIECE_T; j++) {
        syndrome[j] = 0;
        
        for (int i = 0; i < MCELIECE_N; i++) {
            if (vector_get_bit(received, i)) {
                gf_elem_t alpha_i = alpha[i];
                gf_elem_t g_alpha_i = polynomial_eval(g, alpha_i);
                
                if (g_alpha_i != 0) {
                    gf_elem_t alpha_power = gf_pow(alpha_i, j);
                    gf_elem_t g_squared = gf_mul(g_alpha_i, g_alpha_i);
                    gf_elem_t term = gf_div(alpha_power, g_squared);
                    syndrome[j] = gf_add(syndrome[j], term);
                }
            }
        }
    }
}

// 打印矩阵（调试用）
void print_matrix(const char *label, const matrix_t *mat) {
    printf("%s (%dx%d):\n", label, mat->rows, mat->cols);
    
    for (int row = 0; row < mat->rows && row < 16; row++) {  // 打印前 16 行
        printf("    row %2d: ", row);
        for (int col = 0; col < mat->cols && col < 32; col++) {  // 打印前 32 列
            printf("%d", matrix_get_bit(mat, row, col));
        }
        if (mat->cols > 32) printf("...");
        printf("\n");
    }

    if (mat->rows > 16) printf("...\n");
    printf("\n");
}

// 矩阵拷贝
void matrix_copy(matrix_t *dst, const matrix_t *src) {
    if (!dst || !src) return;
    
    int min_rows = (dst->rows < src->rows) ? dst->rows : src->rows;
    int min_cols_bytes = (dst->cols_bytes < src->cols_bytes) ? dst->cols_bytes : src->cols_bytes;
    
    for (int row = 0; row < min_rows; row++) {
        memcpy(dst->data + row * dst->cols_bytes,
               src->data + row * src->cols_bytes,
               min_cols_bytes);
    }
}