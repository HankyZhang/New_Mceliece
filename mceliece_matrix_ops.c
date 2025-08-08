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

    printf("    -> reduce_to_systematic_form: Starting reduction of %dx%d matrix...\n", mt, n);

    // 创建列置换数组
    int *col_perm = malloc(n * sizeof(int));
    if (!col_perm) return -1;
    
    // 初始化列置换
    for (i = 0; i < n; i++) {
        col_perm[i] = i;
    }

    // --- 正向消元，形成上三角矩阵 ---
    for (i = 0; i < mt; i++) {
        if (i % 100 == 0) {
            printf("    -> reduce_to_systematic_form: Processing row %d/%d...\n", i+1, mt);
        }
        
        // 1. 寻找主元 (从第i列开始，从第i行开始找)
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
            // 找不到主元，矩阵是奇异的
            printf("    -> reduce_to_systematic_form: ERROR - No pivot found at row %d, matrix is singular!\n", i);
            free(col_perm);
            return -1;
        }

        // 2. 将主元行交换到当前行i
        if (pivot_row != i) {
            matrix_swap_rows(H, i, pivot_row);
        }

        // 3. 将主元列交换到当前列i
        if (pivot_col != i) {
            matrix_swap_cols(H, i, pivot_col);
            // 更新列置换
            int temp = col_perm[i];
            col_perm[i] = col_perm[pivot_col];
            col_perm[pivot_col] = temp;
        }

        // 4. 将当前列(i)的其他所有行的元素消为0 (在主元下方)
        for (j = i + 1; j < mt; j++) {
            if (matrix_get_bit(H, j, i) == 1) {
                matrix_xor_rows(H, j, i);
            }
        }
    }

    printf("    -> reduce_to_systematic_form: Forward elimination completed.\n");

    // --- 反向消元，形成单位矩阵 ---
    for (i = mt - 1; i >= 0; i--) {
        for (j = 0; j < i; j++) {
            if (matrix_get_bit(H, j, i) == 1) {
                matrix_xor_rows(H, j, i);
            }
        }
    }

    printf("    -> reduce_to_systematic_form: Backward elimination completed.\n");

    free(col_perm);
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

    printf("    -> mat_gen: Creating Goppa parity-check matrix H (%dx%d)...\n", mt, n);
    matrix_t *H = matrix_create(mt, n);
    if (!H) return MCELIECE_ERROR_MEMORY;
    printf("    -> mat_gen: Matrix H created successfully.\n");

    // 根据规范 1.2.7，构造 Goppa 码的校验矩阵
    // M[i,j] = alpha[j]^i / g(alpha[j]) for i=0,...,t-1 and j=0,...,n-1
    printf("    -> mat_gen: Computing Goppa parity-check matrix entries...\n");
    
    // 检查一些样本值
    printf("    -> mat_gen: Sample values - g(alpha[0])=%04x, g(alpha[1])=%04x, g(alpha[2])=%04x\n", 
           polynomial_eval(g, alpha[0]), polynomial_eval(g, alpha[1]), polynomial_eval(g, alpha[2]));
    
    // 检查alpha值
    printf("    -> mat_gen: Sample alpha values - alpha[0]=%04x, alpha[1]=%04x, alpha[2]=%04x\n", 
           alpha[0], alpha[1], alpha[2]);
    
    // 检查多项式系数
    printf("    -> mat_gen: Goppa polynomial g(x) = ");
    for (int i = g->degree; i >= 0; i--) {
        if (g->coeffs[i] != 0) {
            if (i == g->degree) {
                printf("x^%d", i);
            } else if (i == 1) {
                printf(" + %04x*x", g->coeffs[i]);
            } else if (i == 0) {
                printf(" + %04x", g->coeffs[i]);
            } else {
                printf(" + %04x*x^%d", g->coeffs[i], i);
            }
        }
    }
    printf("\n");
    
    // 测试多项式求值
    printf("    -> mat_gen: Testing polynomial evaluation...\n");
    gf_elem_t test_x = 1;
    gf_elem_t test_result = polynomial_eval(g, test_x);
    printf("    -> mat_gen: g(%04x) = %04x\n", test_x, test_result);
    
    // 手动计算 g(1) = 1^128 + 1^7 + 1^2 + 1^1 + 1^0 = 1 + 1 + 1 + 1 + 1 = 5 (0x0005)
    printf("    -> mat_gen: Expected g(1) = 1 + 1 + 1 + 1 + 1 = 5 (0x0005)\n");
    
    // 测试基本GF运算
    printf("    -> mat_gen: Testing basic GF arithmetic...\n");
    printf("    -> mat_gen: gf_add(1, 1) = %04x\n", gf_add(1, 1));
    printf("    -> mat_gen: gf_mul(1, 1) = %04x\n", gf_mul(1, 1));
    printf("    -> mat_gen: gf_mul(2, 3) = %04x\n", gf_mul(2, 3));
    printf("    -> mat_gen: gf_pow(2, 3) = %04x\n", gf_pow(2, 3));
    
    for (int i = 0; i < t; i++) {
        if (i % 32 == 0) {
            printf("    -> mat_gen: Processing row %d/%d...\n", i+1, t);
        }
        for (int j = 0; j < n; j++) {
            // 计算 alpha[j]^i
            gf_elem_t alpha_power = 1;
            for (int pow = 0; pow < i; pow++) {
                alpha_power = gf_mul(alpha_power, alpha[j]);
            }
            
            // 计算 g(alpha[j])
            gf_elem_t g_alpha = polynomial_eval(g, alpha[j]);
            
            // 检查 g(alpha[j]) 是否为零（这会导致除零错误）
            if (g_alpha == 0) {
                printf("    -> mat_gen: ERROR - g(alpha[%d]) = 0, which would cause division by zero!\n", j);
                matrix_free(H);
                return MCELIECE_ERROR_KEYGEN_FAIL;
            }
            
            // 计算 M[i,j] = alpha[j]^i / g(alpha[j])
            gf_elem_t M_ij = gf_div(alpha_power, g_alpha);
            
            // 将 GF(2^m) 元素展开为 m 个二进制位
            for (int bit = 0; bit < m; bit++) {
                int bit_value = (M_ij >> bit) & 1;
                matrix_set_bit(H, i * m + bit, j, bit_value);
            }
        }
    }
    
    printf("    -> mat_gen: Goppa parity-check matrix computed successfully.\n");
    
    // 检查矩阵的一些统计信息
    int non_zero_count = 0;
    for (int i = 0; i < mt; i++) {
        for (int j = 0; j < n; j++) {
            if (matrix_get_bit(H, i, j) == 1) {
                non_zero_count++;
            }
        }
    }
    printf("    -> mat_gen: Matrix statistics - %d non-zero elements out of %d total (%.2f%%)\n", 
           non_zero_count, mt * n, (double)non_zero_count * 100 / (mt * n));

    // 将 H 转换为系统形式 [I_mt | T]
    printf("    -> mat_gen: Reducing H to systematic form...\n");
    if (reduce_to_systematic_form(H) != 0) {
        printf("    -> mat_gen: ERROR - Matrix is singular, key generation failed.\n");
        matrix_free(H);
        return MCELIECE_ERROR_KEYGEN_FAIL; // 矩阵奇异，生成失败
    }
    printf("    -> mat_gen: H reduced to systematic form successfully.\n");

    // 此时 H 的形式是 [I_mt | T]
    // 从 H 的右侧提取公钥 T
    printf("    -> mat_gen: Extracting public key T from H...\n");
    for (int i = 0; i < mt; i++) {
        for (int j = 0; j < k; j++) {
            int bit = matrix_get_bit(H, i, mt + j);
            matrix_set_bit(T_out, i, j, bit);
        }
    }
    printf("    -> mat_gen: Public key T extracted successfully.\n");

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