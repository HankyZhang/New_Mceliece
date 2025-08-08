#include "mceliece_types.h"
#include "matrix.h"
#include "gf.h"
#include "vector.h"
#include "matrix.h"
#include "mceliece.h"
#include <stdio.h>
// 矩阵行交换
void matrix_swap_rows(matrix_t *mat, int row1, int row2) {
    if (row1 == row2 || row1 >= mat->rows || row2 >= mat->rows) return;
    
    for (int col = 0; col < mat->cols_bytes; col++) {
        uint8_t temp = mat->data[row1 * mat->cols_bytes + col];
        mat->data[row1 * mat->cols_bytes + col] = mat->data[row2 * mat->cols_bytes + col];
        mat->data[row2 * mat->cols_bytes + col] = temp;
    }
}
static int reduce_to_systematic_form(matrix_t *H) {
    int mt = H->rows;

    int i, j;

    // --- 正向消元，形成上三角矩阵 ---
    for (i = 0; i < mt; i++) {
        // 1. 寻找主元 (在第i列，从第i行开始找)
        int pivot_row = i;
        while (pivot_row < mt && matrix_get_bit(H, pivot_row, i) == 0) {
            pivot_row++;
        }

        if (pivot_row == mt) {
            // 在这一列找不到主元，矩阵是奇异的，无法转换为系统形式
            return -1; // 失败
        }

        // 2. 将主元行交换到当前行i
        if (pivot_row != i) {
            matrix_swap_rows(H, i, pivot_row);
        }

        // 3. 将当前列(i)的其他所有行的元素消为0 (在主元下方)
        for (j = i + 1; j < mt; j++) {
            if (matrix_get_bit(H, j, i) == 1) {
                matrix_xor_rows(H, j, i);
            }
        }
    }

    // --- 反向消元，形成单位矩阵 ---
    for (i = mt - 1; i >= 0; i--) {
        for (j = 0; j < i; j++) {
            if (matrix_get_bit(H, j, i) == 1) {
                matrix_xor_rows(H, j, i);
            }
        }
    }

    return 0; // 成功
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

    matrix_t *H = matrix_create(mt, n);
    if (!H) return MCELIECE_ERROR_MEMORY;

    // 构建原始校验矩阵 H (这部分逻辑您已经写对了)
    for (int j = 0; j < n; j++) {
        gf_elem_t eval_g = polynomial_eval(g, alpha[j]);
        if (eval_g == 0) {
            matrix_free(H);
            return MCELIECE_ERROR_KEYGEN_FAIL; // 理论上不应发生，因为alpha是支持集
        }
        gf_elem_t inv_g = gf_inv(eval_g);
        gf_elem_t alpha_pow_j = 1;

        for (int row_t = 0; row_t < t; row_t++) {
            gf_elem_t element = gf_mul(alpha_pow_j, inv_g);
            for (int bit_m = 0; bit_m < m; bit_m++) {
                int row = row_t * m + bit_m;
                if ((element >> bit_m) & 1) {
                    matrix_set_bit(H, row, j, 1);
                }
            }
            alpha_pow_j = gf_mul(alpha_pow_j, alpha[j]);
        }
    }

    // 调用新的、符合规范的高斯消元函数
    if (reduce_to_systematic_form(H) != 0) {
        matrix_free(H);
        return MCELIECE_ERROR_KEYGEN_FAIL; // 矩阵奇异，生成失败
    }

    // 此时 H 的形式是 [I_mt | T]
    // 从 H 的右侧提取公钥 T
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

// 计算伴随式（syndrome）
void compute_syndrome(const uint8_t *received, const polynomial_t *g, 
                     const gf_elem_t *alpha, gf_elem_t *syndrome) {
    if (!received || !g || !alpha || !syndrome) return;
    
    // 计算syndrome[j] = sum_{i in error_positions} alpha[i]^j / g(alpha[i])^2
    // 但我们不知道错误位置，所以需要从接收到的向量计算
    
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