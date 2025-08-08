#include "mceliece_types.h"
#include "mceliece_decode.h"
#include "matrix.h"
#include "gf.h"
#include  "vector.h"

// Berlekamp-Massey算法：求解线性反馈移位寄存器
// 输入：syndrome序列s[0], s[1], ..., s[2t-1]
// 输出：错误定位多项式sigma和错误求值多项式omega
mceliece_error_t berlekamp_massey(const gf_elem_t *syndrome, 
                                 polynomial_t *sigma, polynomial_t *omega) {
    if (!syndrome || !sigma || !omega) {
        return MCELIECE_ERROR_INVALID_PARAM;
    }
    
    // 初始化多项式
    polynomial_t *C = polynomial_create(MCELIECE_T);  // 当前连接多项式
    polynomial_t *B = polynomial_create(MCELIECE_T);  // 备用多项式
    polynomial_t *T = polynomial_create(MCELIECE_T);  // 临时多项式
    
    if (!C || !B || !T) {
        polynomial_free(C);
        polynomial_free(B);
        polynomial_free(T);
        return MCELIECE_ERROR_MEMORY;
    }
    
    // 初始状态：C(x) = 1, B(x) = 1
    polynomial_set_coeff(C, 0, 1);
    polynomial_set_coeff(B, 0, 1);
    
    int L = 0;          // 当前LFSR长度
    int m = 1;          // 修正量计数器
    gf_elem_t b = 1;    // 上一个最佳偏差
    
    // 对每个syndrome元素进行迭代
    for (int N = 0; N < 2 * MCELIECE_T; N++) {
        // 计算偏差d_N
        gf_elem_t d = syndrome[N];
        
        for (int i = 1; i <= L; i++) {
            if (C->coeffs[i] != 0) {
                d = gf_add(d, gf_mul(C->coeffs[i], syndrome[N - i]));
            }
        }
        
        if (d == 0) {
            // 偏差为0，不需要修正
            m++;
        } else {
            // 偏差非0，需要修正
            
            // 保存当前的C到T：T(x) = C(x)
            polynomial_copy(T, C);
            
            // 修正：C(x) = C(x) - (d/b) * x^m * B(x)
            if (b != 0) {
                gf_elem_t correction_coeff = gf_div(d, b);
                
                for (int i = 0; i <= B->degree && i + m <= C->max_degree; i++) {
                    if (B->coeffs[i] != 0) {
                        gf_elem_t term = gf_mul(correction_coeff, B->coeffs[i]);
                        gf_elem_t new_coeff = gf_add(C->coeffs[i + m], term);
                        polynomial_set_coeff(C, i + m, new_coeff);
                    }
                }
            }
            
            // 判断是否需要更新L
            if (2 * L <= N) {
                L = N + 1 - L;
                polynomial_copy(B, T);  // B(x) = T(x) (即原来的C(x))
                b = d;
                m = 1;
            } else {
                m++;
            }
        }
    }
    
    // 输出错误定位多项式
    polynomial_copy(sigma, C);
    
    // 计算错误求值多项式omega
    // omega(x) = sigma(x) * S(x) mod x^(2t)
    // 其中S(x) = s[0] + s[1]*x + s[2]*x^2 + ...
    
    // 构造syndrome多项式S(x)
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
    
    // 计算sigma * S的前2t项
    memset(omega->coeffs, 0, (omega->max_degree + 1) * sizeof(gf_elem_t));
    omega->degree = -1;
    
    for (int i = 0; i <= sigma->degree && i < 2 * MCELIECE_T; i++) {
        if (sigma->coeffs[i] != 0) {
            for (int j = 0; j <= S->degree && i + j < 2 * MCELIECE_T; j++) {
                if (S->coeffs[j] != 0 && i + j <= omega->max_degree) {
                    gf_elem_t term = gf_mul(sigma->coeffs[i], S->coeffs[j]);
                    gf_elem_t new_coeff = gf_add(omega->coeffs[i + j], term);
                    polynomial_set_coeff(omega, i + j, new_coeff);
                }
            }
        }
    }
    
    // 清理
    polynomial_free(C);
    polynomial_free(B);
    polynomial_free(T);
    polynomial_free(S);
    
    return MCELIECE_SUCCESS;
}

// Chien搜索：寻找错误定位多项式的根
mceliece_error_t chien_search(const polynomial_t *sigma, const gf_elem_t *alpha,
                             int *error_positions, int *num_errors) {
    if (!sigma || !alpha || !error_positions || !num_errors) {
        return MCELIECE_ERROR_INVALID_PARAM;
    }
    
    *num_errors = 0;
    
    // 对支撑集中的每个元素alpha[j]，检查sigma(1/alpha[j])是否为0
    for (int j = 0; j < MCELIECE_N; j++) {
        if (alpha[j] == 0) continue;  // 跳过0元素
        
        gf_elem_t alpha_inv = gf_inv(alpha[j]);
        gf_elem_t result = polynomial_eval(sigma, alpha_inv);
        
        if (result == 0) {
            // 找到一个根，对应错误位置
            error_positions[*num_errors] = j;
            (*num_errors)++;
            
            if (*num_errors >= MCELIECE_T) break;  // 最多t个错误
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
    
    // 检查syndrome是否全零（无错误）
    int has_error = 0;
    for (int i = 0; i < 2 * MCELIECE_T; i++) {
        if (syndrome[i] != 0) {
            has_error = 1;
            break;
        }
    }
    
    if (!has_error) {
        // 无错误
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
    
    // 检查错误个数是否正确
    if (num_errors != sigma->degree) {
        // 解码失败
        free(syndrome);
        polynomial_free(sigma);
        polynomial_free(omega);
        free(error_positions);
        return MCELIECE_SUCCESS;  // 不是错误，只是解码失败
    }
    
    // 步骤4：构造错误向量
    memset(error_vector, 0, MCELIECE_N_BYTES);

    for (int i = 0; i < num_errors; i++) {
        // 将旧的调用: vector_set_bit(error_vector, error_positions[i]);
        // 修改为新的调用，并添加第三个参数 '1'
        vector_set_bit(error_vector, error_positions[i], 1);
    }
    
    *decode_success = 1;
    
    // 清理内存
    free(syndrome);
    polynomial_free(sigma);
    polynomial_free(omega);
    free(error_positions);
    
    return MCELIECE_SUCCESS;
}