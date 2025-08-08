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
// 用列交换实现的系统形式化，目标把左侧mt列变成单位阵
static int reduce_to_systematic_form(matrix_t *H) {
    int mt = H->rows;
    int n = H->cols;

    // 当前处理的行、列
    int row = 0;
    for (int col = 0; col < n && row < mt; col++) {
        // 1) 在当前列col寻找一个主元行 >= row
        int pivot_row = -1;
        for (int r = row; r < mt; r++) {
            if (matrix_get_bit(H, r, col)) { pivot_row = r; break; }
        }

        // 如果当前列没有主元，尝试在右侧列中找到一个含1的列并交换到当前列
        if (pivot_row == -1) {
            int pivot_col = -1;
            for (int c = col + 1; c < n; c++) {
                for (int r = row; r < mt; r++) {
                    if (matrix_get_bit(H, r, c)) { pivot_col = c; break; }
                }
                if (pivot_col != -1) break;
            }
            if (pivot_col == -1) {
                // 无法在后续列中找到主元，继续让col前进
                continue;
            }
            // 交换列，把可作主元的列移到当前位置
            matrix_swap_cols(H, col, pivot_col);
            // 重新在该列寻找主元行
            for (int r = row; r < mt; r++) {
                if (matrix_get_bit(H, r, col)) { pivot_row = r; break; }
            }
            if (pivot_row == -1) {
                // 理论上不会发生
                continue;
            }
        }

        // 2) 把主元行交换到当前row
        if (pivot_row != row) {
            matrix_swap_rows(H, row, pivot_row);
        }

        // 3) 将该列的其他行清零
        for (int r = 0; r < mt; r++) {
            if (r == row) continue;
            if (matrix_get_bit(H, r, col)) {
                matrix_xor_rows(H, r, row);
            }
        }

        // 推进到下一行
        row++;
    }

    // 检查是否成功在左侧形成了mt阶单位阵（即前mt列为单位阵）
    // 如未形成，返回失败
    int identity_cols_found = 0;
    // 逐列找单位列向量
    for (int c = 0; c < mt; c++) {
        int ones = 0, one_row = -1;
        for (int r = 0; r < mt; r++) {
            int bit = matrix_get_bit(H, r, c);
            if (bit) { ones++; one_row = r; if (ones > 1) break; }
        }
        if (ones == 1) {
            // 确认为单位列，且该1位应位于对角位置，但我们只是要求单位性
            identity_cols_found++;
        } else {
            // 尝试在后续列中寻找单位列并交换到c位
            int id_col = -1;
            for (int cc = c + 1; cc < n; cc++) {
                int ones2 = 0, one_row2 = -1;
                for (int r = 0; r < mt; r++) {
                    int bit = matrix_get_bit(H, r, cc);
                    if (bit) { ones2++; one_row2 = r; if (ones2 > 1) break; }
                }
                if (ones2 == 1) { id_col = cc; break; }
            }
            if (id_col == -1) {
                return -1; // 无法形成系统形式
            }
            matrix_swap_cols(H, c, id_col);
            identity_cols_found++;
        }
    }

    return (identity_cols_found == mt) ? 0 : -1;
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

    // 调用新的、包含列交换的系统化函数
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