// mceliece_random.h
#ifndef MCELIECE_RANDOM_H
#define MCELIECE_RANDOM_H

#include "mceliece_types.h" // 需要 mceliece_error_t 和 polynomial_t

#ifdef __cplusplus
extern "C" {
#endif

    // 函数声明
    void init_random(void); // 如果您还在用旧的随机源
    void random_bytes(uint8_t *output, int len);

    mceliece_error_t fixed_weight_vector(uint8_t *output, int vector_len, int target_weight);
    mceliece_error_t generate_field_ordering(gf_elem_t *alpha_output, const uint8_t *random_bits);
    mceliece_error_t generate_irreducible_poly_final(polynomial_t *g, const uint8_t *random_bits);
#ifdef __cplusplus
}
#endif

#endif // MCELIECE_RANDOM_H