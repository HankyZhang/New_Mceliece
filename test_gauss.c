// test_gauss.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mceliece_types.h"
#include "gf.h"
#include "matrix.h"

// Function declarations
extern int reduce_to_systematic_form(matrix_t *H);
extern int matrix_is_systematic(const matrix_t *mat);







// --- 单元测试主函数 ---



void test_random_matrix_reduction() {
    printf("=== Testing Random Matrix Reduction ===\n");
    
    // Test with McEliece parameters
    int mt = 1664;  // m * t = 13 * 128
    int n = 6688;   // 2^m = 2^13
    
    printf("Testing with McEliece parameters: %dx%d matrix\n", mt, n);
    
    matrix_t *H = matrix_create(mt, n);
    if (!H) {
        printf("Failed to create matrix\n");
        return;
    }
    
    // Fill matrix with random bits
    srand(time(NULL));
    for (int i = 0; i < mt; i++) {
        for (int j = 0; j < n; j++) {
            if (rand() % 2 == 1) {
                matrix_set_bit(H, i, j, 1);
            }
        }
    }
    
    printf("Random matrix created. Attempting reduction...\n");
    
    // Try to reduce to systematic form
    int result = reduce_to_systematic_form(H);
    printf("Reduction result: %d\n", result);
    
    if (result == 0) {
        printf("Reduction successful!\n");
        
        // Check if it's systematic
        int is_sys = matrix_is_systematic(H);
        printf("Is systematic: %d\n", is_sys);
        
        // Print first few rows and columns
        printf("First 4x8 elements of reduced matrix:\n");
        for (int i = 0; i < 4 && i < mt; i++) {
            for (int j = 0; j < 8 && j < n; j++) {
                printf("%d ", matrix_get_bit(H, i, j));
            }
            printf("\n");
        }
    } else {
        printf("Reduction failed: matrix is singular.\n");
    }
    
    matrix_free(H);
}

void test_matrix_construction() {
    printf("=== Testing Matrix Construction ===\n");
    
    // Test with a small matrix first
    int mt = 4;  // 4 rows
    int n = 8;   // 8 columns
    matrix_t *H = matrix_create(mt, n);
    
    if (!H) {
        printf("Failed to create matrix\n");
        return;
    }
    
    // Create a simple test matrix that should be full rank
    // [1 0 0 0 1 1 0 1]
    // [0 1 0 0 1 0 1 1]
    // [0 0 1 0 0 1 1 1]
    // [0 0 0 1 1 1 1 0]
    
    matrix_set_bit(H, 0, 0, 1);
    matrix_set_bit(H, 0, 4, 1);
    matrix_set_bit(H, 0, 5, 1);
    matrix_set_bit(H, 0, 7, 1);
    
    matrix_set_bit(H, 1, 1, 1);
    matrix_set_bit(H, 1, 4, 1);
    matrix_set_bit(H, 1, 6, 1);
    matrix_set_bit(H, 1, 7, 1);
    
    matrix_set_bit(H, 2, 2, 1);
    matrix_set_bit(H, 2, 5, 1);
    matrix_set_bit(H, 2, 6, 1);
    matrix_set_bit(H, 2, 7, 1);
    
    matrix_set_bit(H, 3, 3, 1);
    matrix_set_bit(H, 3, 4, 1);
    matrix_set_bit(H, 3, 5, 1);
    matrix_set_bit(H, 3, 6, 1);
    
    printf("Original matrix:\n");
    for (int i = 0; i < mt; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", matrix_get_bit(H, i, j));
        }
        printf("\n");
    }
    
    // Try to reduce to systematic form
    int result = reduce_to_systematic_form(H);
    printf("Reduction result: %d\n", result);
    
    if (result == 0) {
        printf("Reduced matrix:\n");
        for (int i = 0; i < mt; i++) {
            for (int j = 0; j < n; j++) {
                printf("%d ", matrix_get_bit(H, i, j));
            }
            printf("\n");
        }
        
        // Check if it's systematic
        int is_sys = matrix_is_systematic(H);
        printf("Is systematic: %d\n", is_sys);
    }
    
    matrix_free(H);
}

void test_basic_operations() {
    printf("=== Testing Basic Operations ===\n");
    
    // Test GF operations
    printf("Testing GF operations...\n");
    gf_elem_t a = 5, b = 3;
    gf_elem_t sum = gf_add(a, b);
    gf_elem_t prod = gf_mul(a, b);
    gf_elem_t inv_a = gf_inv(a);
    
    printf("a = %d, b = %d\n", a, b);
    printf("a + b = %d\n", sum);
    printf("a * b = %d\n", prod);
    printf("a^(-1) = %d\n", inv_a);
    printf("a * a^(-1) = %d\n", gf_mul(a, inv_a));
    
    // Test polynomial evaluation
    printf("\nTesting polynomial evaluation...\n");
    polynomial_t *poly = polynomial_create(3);
    polynomial_set_coeff(poly, 0, 1);  // constant term
    polynomial_set_coeff(poly, 1, 2);  // x term
    polynomial_set_coeff(poly, 2, 1);  // x^2 term
    polynomial_set_coeff(poly, 3, 1);  // x^3 term
    
    gf_elem_t eval_result = polynomial_eval(poly, 3);
    printf("f(x) = 1 + 2x + x^2 + x^3\n");
    printf("f(3) = %d\n", eval_result);
    
    polynomial_free(poly);
}

int main() {
    gf_init();  // Initialize the finite field
    
    test_basic_operations();
    printf("\n");
    test_matrix_construction();
    printf("\n");
    test_random_matrix_reduction();
    return 0;
}