#ifndef MCELIECE_TEST_H
#define MCELIECE_TEST_H

#include "mceliece_types.h"

#ifdef __cplusplus
extern "C" {
#endif

// Main test function
void run_all_tests(void);

// Individual test functions
void test_gf_operations(void);
void test_polynomial_operations(void);
void test_irreducible_polynomial(void);
void test_field_ordering(void);
void test_matrix_operations(void);
void test_fixed_weight_vector(void);
void test_berlekamp_massey(void);
void test_kem_operations(void);
void test_known_vectors_and_edge_cases(void);
void test_performance_and_stress(void);

// Test utility functions
void test_assert(int condition, const char *test_name);
void print_test_summary(void);

// Functions are declared in mceliece_random.h

#ifdef __cplusplus
}
#endif

#endif // MCELIECE_TEST_H
