// test_gauss.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mceliece_types.h"
#include "gf.h"
#include "matrix.h"

// --- 从 mceliece_matrix_ops.c 中复制过来的函数 ---
// 我们需要在这里重新定义它们，或者创建一个共享的 .h 文件
// 为了简单，我们直接在这里定义

void matrix_swap_rows(matrix_t *mat, int row1, int row2) {
    if (row1 == row2) return;
    for (int col = 0; col < mat->cols_bytes; col++) {
        uint8_t temp = mat->data[row1 * mat->cols_bytes + col];
        mat->data[row1 * mat->cols_bytes + col] = mat->data[row2 * mat->cols_bytes + col];
        mat->data[row2 * mat->cols_bytes + col] = temp;
    }
}

void matrix_swap_cols(matrix_t *mat, int col1, int col2) {
    if (col1 == col2) return;
    for (int row = 0; row < mat->rows; row++) {
        int bit1 = matrix_get_bit(mat, row, col1);
        int bit2 = matrix_get_bit(mat, row, col2);
        matrix_set_bit(mat, row, col1, bit2);
        matrix_set_bit(mat, row, col2, bit1);
    }
}

void matrix_xor_rows(matrix_t *mat, int row_dst, int row_src) {
    for (int col = 0; col < mat->cols_bytes; col++) {
        mat->data[row_dst * mat->cols_bytes + col] ^= mat->data[row_src * mat->cols_bytes + col];
    }
}

void print_matrix(const char *label, const matrix_t *mat) {
    printf("%s (%dx%d):\n", label, mat->rows, mat->cols);
    for (int row = 0; row < mat->rows; row++) {
        printf("  ");
        for (int col = 0; col < mat->cols; col++) {
            printf("%d", matrix_get_bit(mat, row, col));
            if ((col + 1) % 8 == 0 && col < mat->cols - 1) printf(" ");
        }
        printf("\n");
    }
    printf("\n");
}


// 这是我们要测试的目标函数 - pk_gen
// (请确保这里的代码与 mceliece_matrix_ops.c 中的版本完全一致)
static int pk_gen(matrix_t *T, int *p, matrix_t *H) {
    int mt = H->rows;
    int n = H->cols;
    int k = n - mt;
    int i, j, row, col, temp;

    matrix_t *H_copy = matrix_create(mt, n);
    if (!H_copy) return -1;
    matrix_copy(H_copy, H);

    for (i = 0; i < n; i++) {
        p[i] = i;
    }

    for (i = 0; i < mt; i++) {
        int pivot_found = 0;
        for (j = i; j < n; j++) {
            for (row = i; row < mt; row++) {
                if (matrix_get_bit(H_copy, row, j) == 1) {
                    if (i != row) {
                        matrix_swap_rows(H_copy, i, row);
                    }
                    if (i != j) {
                        matrix_swap_cols(H_copy, i, j);
                        temp = p[i]; p[i] = p[j]; p[j] = temp;
                    }
                    pivot_found = 1;
                    break;
                }
            }
            if (pivot_found) {
                break;
            }
        }

        if (!pivot_found) {
            matrix_free(H_copy);
            return -1;
        }

        for (row = 0; row < mt; row++) {
            if (row != i && matrix_get_bit(H_copy, row, i) == 1) {
                matrix_xor_rows(H_copy, row, i);
            }
        }
    }

    for (i = 0; i < mt; i++) {
        for (j = 0; j < k; j++) {
            matrix_set_bit(T, i, j, matrix_get_bit(H_copy, i, mt + j));
        }
    }

    matrix_free(H_copy);
    return 0;
}


// --- 单元测试主函数 ---

// 创建一个随机的、保证满秩的 mt x n 矩阵
matrix_t* create_test_matrix(int mt, int n) {
    matrix_t *mat = matrix_create(mt, n);
    if (!mat) return NULL;

    // 1. 先创建一个系统形式的矩阵 [I_mt | R]，它必然是满秩的
    for (int i = 0; i < mt; i++) {
        matrix_set_bit(mat, i, i, 1);
    }
    for (int i = 0; i < mt; i++) {
        for (int j = mt; j < n; j++) {
            if (rand() % 2 == 1) {
                matrix_set_bit(mat, i, j, 1);
            }
        }
    }

    // 2. 对这个矩阵进行随机的行和列变换，打乱它
    for (int i = 0; i < 1000; i++) { // 1000 次随机变换
        int r1 = rand() % mt;
        int r2 = rand() % mt;
        matrix_swap_rows(mat, r1, r2);

        int c1 = rand() % n;
        int c2 = rand() % n;
        matrix_swap_cols(mat, c1, c2);

        if (rand() % 3 == 0) {
            matrix_xor_rows(mat, r1, r2);
        }
    }

    return mat;
}

int main(void) {
    srand(time(NULL));
    int num_tests = 100;
    int successes = 0;

    // 使用小一点的参数进行快速测试
    int test_mt = 16;
    int test_n = 40;
    int test_k = test_n - test_mt;

    printf("--- Running Gauss-Jordan Elimination Test ---\n");
    printf("Testing with %d random %dx%d full-rank matrices.\n\n", num_tests, test_mt, test_n);

    for (int i = 0; i < num_tests; i++) {
        printf("Test #%d:\n", i + 1);
        matrix_t *H = create_test_matrix(test_mt, test_n);
        if (!H) {
            printf("  Failed to create test matrix.\n");
            continue;
        }

        // print_matrix("  Original random matrix H:", H);

        matrix_t *T = matrix_create(test_mt, test_k);
        int *p = malloc(test_n * sizeof(int));
        if (!T || !p) {
            printf("  Memory allocation failed.\n");
            matrix_free(H);
            matrix_free(T);
            free(p);
            continue;
        }

        int result = pk_gen(T, p, H);

        if (result == 0) {
            printf("  pk_gen PASSED.\n");
            successes++;
        } else {
            printf("  pk_gen FAILED! Returned -1.\n");
            printf("  This matrix was incorrectly identified as singular:\n");
            print_matrix("  Failed Matrix H:", H);
        }

        matrix_free(H);
        matrix_free(T);
        free(p);
    }

    printf("\n--- Test Summary ---\n");
    printf("Successes: %d / %d\n", successes, num_tests);
    if (successes == num_tests) {
        printf("All tests passed! The pk_gen function seems correct.\n");
    } else {
        printf("There are failures! The pk_gen function has a bug.\n");
    }

    return 0;
}