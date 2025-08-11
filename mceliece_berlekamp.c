#include "mceliece_types.h"
#include "mceliece_decode.h"
#include "matrix.h"
#include "gf.h"
#include  "vector.h"

// Berlekamp-Massey Algorithm according to Classic McEliece specification
// Input: syndrome sequence s[0], s[1], ..., s[2t-1]
// Output: error locator polynomial sigma and error evaluator polynomial omega
mceliece_error_t berlekamp_massey(const gf_elem_t *syndrome, 
                                 polynomial_t *sigma, polynomial_t *omega) {
    if (!syndrome || !sigma || !omega) {
        return MCELIECE_ERROR_INVALID_PARAM;
    }
    
    // Initialize polynomials
    polynomial_t *C = polynomial_create(MCELIECE_T);  // Current connection polynomial
    polynomial_t *B = polynomial_create(MCELIECE_T);  // Backup polynomial
    polynomial_t *T = polynomial_create(MCELIECE_T);  // Temporary polynomial
    
    if (!C || !B || !T) {
        polynomial_free(C);
        polynomial_free(B);
        polynomial_free(T);
        return MCELIECE_ERROR_MEMORY;
    }
    
    // Initial state: C(x) = 1, B(x) = 1
    polynomial_set_coeff(C, 0, 1);
    polynomial_set_coeff(B, 0, 1);
    
    int L = 0;          // Current LFSR length
    int m = 1;          // Step counter since last L update
    gf_elem_t b = 1;    // Last best discrepancy
    
    // Iterate through each syndrome element
    for (int N = 0; N < 2 * MCELIECE_T; N++) {
        // Calculate discrepancy d_N = s_N + Σ C_i * s_{N-i}
        gf_elem_t d = syndrome[N];
        
        for (int i = 1; i <= L && (N - i) >= 0; i++) {
            if (i <= C->degree && C->coeffs[i] != 0) {
                d = gf_add(d, gf_mul(C->coeffs[i], syndrome[N - i]));
            }
        }
        
        if (d == 0) {
            // Discrepancy is 0, no correction needed
            m++;
        } else {
            // Discrepancy is non-zero, correction needed
            
            // Save current C to T: T(x) = C(x)
            polynomial_copy(T, C);
            
            // Correction: C(x) = C(x) - (d/b) * x^m * B(x)
            if (b != 0) {
                gf_elem_t correction_coeff = gf_div(d, b);
                
                for (int i = 0; i <= B->degree; i++) {
                    if (B->coeffs[i] != 0 && (i + m) <= C->max_degree) {
                        gf_elem_t term = gf_mul(correction_coeff, B->coeffs[i]);
                        gf_elem_t current_coeff = (i + m <= C->degree) ? C->coeffs[i + m] : 0;
                        gf_elem_t new_coeff = gf_add(current_coeff, term);
                        polynomial_set_coeff(C, i + m, new_coeff);
                    }
                }
            }
            
            // Check if L needs to be updated
            if (2 * L <= N) {
                L = N + 1 - L;
                polynomial_copy(B, T);  // B(x) = T(x) (the old C(x))
                b = d;
                m = 1;
            } else {
                m++;
            }
        }
    }
    
    // Output error locator polynomial
    polynomial_copy(sigma, C);
    
    // Calculate error evaluator polynomial omega
    // omega(x) = sigma(x) * S(x) mod x^(2t)
    // where S(x) = s[0] + s[1]*x + s[2]*x^2 + ...
    
    // Construct syndrome polynomial S(x)
    polynomial_t *S = polynomial_create(2 * MCELIECE_T - 1);
    if (!S) {
        polynomial_free(C);
        polynomial_free(B);
        polynomial_free(T);
        return MCELIECE_ERROR_MEMORY;
    }
    
    for (int i = 0; i < 2 * MCELIECE_T; i++) {
        polynomial_set_coeff(S, i, syndrome[i]);
    }
    
    // Calculate first 2t terms of sigma * S
    memset(omega->coeffs, 0, (omega->max_degree + 1) * sizeof(gf_elem_t));
    omega->degree = -1;
    
    for (int i = 0; i <= sigma->degree; i++) {
        if (sigma->coeffs[i] != 0) {
            for (int j = 0; j <= S->degree && (i + j) < 2 * MCELIECE_T; j++) {
                if (S->coeffs[j] != 0 && (i + j) <= omega->max_degree) {
                    gf_elem_t term = gf_mul(sigma->coeffs[i], S->coeffs[j]);
                    gf_elem_t current_coeff = ((i + j) <= omega->degree) ? omega->coeffs[i + j] : 0;
                    gf_elem_t new_coeff = gf_add(current_coeff, term);
                    polynomial_set_coeff(omega, i + j, new_coeff);
                }
            }
        }
    }
    
    // Cleanup
    polynomial_free(C);
    polynomial_free(B);
    polynomial_free(T);
    polynomial_free(S);
    
    return MCELIECE_SUCCESS;
}

// Chien Search: Find roots of error locator polynomial
// According to section 3.3.1: Check σ(α_j^{-1}) = 0 for error locations
mceliece_error_t chien_search(const polynomial_t *sigma, const gf_elem_t *alpha,
                             int *error_positions, int *num_errors) {
    if (!sigma || !alpha || !error_positions || !num_errors) {
        return MCELIECE_ERROR_INVALID_PARAM;
    }
    
    *num_errors = 0;
    
    // For each element alpha[j] in the support set, check if σ(1/α_j) = 0
    for (int j = 0; j < MCELIECE_N; j++) {
        if (alpha[j] == 0) continue;  // Skip zero elements
        
        gf_elem_t alpha_inv = gf_inv(alpha[j]);
        gf_elem_t result = polynomial_eval(sigma, alpha_inv);
        
        if (result == 0) {
            // Found a root, corresponding to error position
            error_positions[*num_errors] = j;
            (*num_errors)++;
            
            if (*num_errors >= MCELIECE_T) break;  // At most t errors
        }
    }
    
    return MCELIECE_SUCCESS;
}

// Forney算法：计算错误值（对于二元码，错误值总是1）
mceliece_error_t forney_algorithm(const polynomial_t *sigma, const polynomial_t *omega,
                                 const gf_elem_t *alpha, const int *error_positions, 
                                 int num_errors, gf_elem_t *error_values) {
    if (!sigma || !omega || !alpha || !error_positions || !error_values) {
        return MCELIECE_ERROR_INVALID_PARAM;
    }
    
    // 计算sigma的形式导数
    polynomial_t *sigma_prime = polynomial_create(sigma->max_degree);
    if (!sigma_prime) return MCELIECE_ERROR_MEMORY;
    
    // 在GF(2^m)中，多项式的导数只保留奇次项系数
    for (int i = 1; i <= sigma->degree; i += 2) {
        if (i - 1 <= sigma_prime->max_degree) {
            polynomial_set_coeff(sigma_prime, i - 1, sigma->coeffs[i]);
        }
    }
    
    // 对每个错误位置计算错误值
    for (int k = 0; k < num_errors; k++) {
        int pos = error_positions[k];
        gf_elem_t alpha_k = alpha[pos];
        
        if (alpha_k == 0) continue;
        
        gf_elem_t alpha_inv = gf_inv(alpha_k);
        
        // 计算omega(alpha_k^{-1})
        gf_elem_t omega_val = polynomial_eval(omega, alpha_inv);
        
        // 计算sigma'(alpha_k^{-1})
        gf_elem_t sigma_prime_val = polynomial_eval(sigma_prime, alpha_inv);
        
        if (sigma_prime_val == 0) {
            polynomial_free(sigma_prime);
            return MCELIECE_ERROR_DECODE_FAIL;
        }
        
        // 错误值 = alpha_k * omega(alpha_k^{-1}) / sigma'(alpha_k^{-1})
        gf_elem_t error_val = gf_div(gf_mul(alpha_k, omega_val), sigma_prime_val);
        error_values[k] = error_val;
    }
    
    polynomial_free(sigma_prime);
    return MCELIECE_SUCCESS;
}

// 完整的解码算法
mceliece_error_t decode_goppa(const uint8_t *received, const polynomial_t *g, 
                             const gf_elem_t *alpha, uint8_t *error_vector, 
                             int *decode_success) {
    if (!received || !g || !alpha || !error_vector || !decode_success) {
        return MCELIECE_ERROR_INVALID_PARAM;
    }
    
    *decode_success = 0;
    
    // 步骤1：计算伴随式
    gf_elem_t *syndrome = malloc(2 * MCELIECE_T * sizeof(gf_elem_t));
    if (!syndrome) return MCELIECE_ERROR_MEMORY;
    
    compute_syndrome(received, g, alpha, syndrome);
    
    // Check if syndrome is all zero (no errors)
    int has_error = 0;
    for (int i = 0; i < 2 * MCELIECE_T; i++) {
        if (syndrome[i] != 0) {
            has_error = 1;
            break;
        }
    }
    
    if (!has_error) {
        // No errors
        memset(error_vector, 0, MCELIECE_N_BYTES);
        *decode_success = 1;
        free(syndrome);
        return MCELIECE_SUCCESS;
    }
    
    // 步骤2：使用Berlekamp-Massey算法求解错误定位多项式
    polynomial_t *sigma = polynomial_create(MCELIECE_T);
    polynomial_t *omega = polynomial_create(MCELIECE_T - 1);
    
    if (!sigma || !omega) {
        free(syndrome);
        polynomial_free(sigma);
        polynomial_free(omega);
        return MCELIECE_ERROR_MEMORY;
    }
    
    mceliece_error_t ret = berlekamp_massey(syndrome, sigma, omega);
    if (ret != MCELIECE_SUCCESS) {
        free(syndrome);
        polynomial_free(sigma);
        polynomial_free(omega);
        return ret;
    }
    
    // 步骤3：使用Chien搜索找到错误位置
    int *error_positions = malloc(MCELIECE_T * sizeof(int));
    if (!error_positions) {
        free(syndrome);
        polynomial_free(sigma);
        polynomial_free(omega);
        return MCELIECE_ERROR_MEMORY;
    }
    
    int num_errors;
    ret = chien_search(sigma, alpha, error_positions, &num_errors);
    if (ret != MCELIECE_SUCCESS) {
        free(syndrome);
        polynomial_free(sigma);
        polynomial_free(omega);
        free(error_positions);
        return ret;
    }
    
    // Check if the number of errors matches the degree of sigma
    // For binary codes, this should be equal for successful decoding
    if (num_errors == 0 || num_errors != sigma->degree) {
        // Decoding failed - no errors found or mismatch
        *decode_success = 0;
        free(syndrome);
        polynomial_free(sigma);
        polynomial_free(omega);
        free(error_positions);
        return MCELIECE_SUCCESS;  // Not an error, just decoding failure
    }
    
    // Step 4: Construct error vector
    memset(error_vector, 0, MCELIECE_N_BYTES);

    for (int i = 0; i < num_errors; i++) {
        // Validate error position
        if (error_positions[i] >= 0 && error_positions[i] < MCELIECE_N) {
            vector_set_bit(error_vector, error_positions[i], 1);
        } else {
            // Invalid error position, decoding failed
            *decode_success = 0;
            free(syndrome);
            polynomial_free(sigma);
            polynomial_free(omega);
            free(error_positions);
            return MCELIECE_SUCCESS;
        }
    }
    
    // Final validation: check if we have exactly t errors
    int actual_weight = vector_weight(error_vector, MCELIECE_N_BYTES);
    if (actual_weight == num_errors && num_errors <= MCELIECE_T) {
        *decode_success = 1;
    } else {
        *decode_success = 0;
    }
    
    // 清理内存
    free(syndrome);
    polynomial_free(sigma);
    polynomial_free(omega);
    free(error_positions);
    
    return MCELIECE_SUCCESS;
}