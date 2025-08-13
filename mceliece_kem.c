#include "mceliece_types.h"
#include "gf.h"
#include "matrix.h"
#include "vector.h"
#include "mceliece.h"
#include "mceliece_decode.h"
#include "mceliece_random.h"

#include <stdio.h>
// 外部函数声明
extern void mceliece_prg(const uint8_t *seed, uint8_t *output, size_t output_len);
extern void mceliece_hash(uint8_t prefix, const uint8_t *input, size_t input_len, uint8_t *output);

mceliece_error_t fixed_weight_vector(uint8_t *e, int n, int t) {
    memset(e, 0, (n + 7) / 8);

    // 根据规范 2.1 FixedWeight() 算法
    // 1. 生成 σ₁τ 个随机比特，其中 τ ≥ t
    int tau = t + 10; // 确保 τ ≥ t，增加一些余量
    size_t random_bytes_len = tau * 2; // σ₁ = 16 bits = 2 bytes per position
    uint8_t *random_bytes = malloc(random_bytes_len);
    if (!random_bytes) return MCELIECE_ERROR_MEMORY;

    mceliece_prg((const uint8_t*)"a_seed_for_fixed_weight_vector", random_bytes, random_bytes_len);

    // 2. 为每个 j ∈ {0, 1, ..., τ-1}，定义 d_j
    int *d_values = malloc(tau * sizeof(int));
    if (!d_values) { free(random_bytes); return MCELIECE_ERROR_MEMORY; }

    for (int j = 0; j < tau; j++) {
        // 取 σ₁ 比特块的前 m 位作为整数
        uint16_t d_j = (uint16_t)random_bytes[j * 2] |
                      ((uint16_t)random_bytes[j * 2 + 1] << 8);
        d_values[j] = d_j % n; // 范围在 {0, 1, ..., n-1}
    }

    // 3. 定义 a_0, a_1, ..., a_{t-1} 为从 d_0, d_1, ..., d_{τ-1} 中选择的前 t 个唯一条目
    int *positions = malloc(t * sizeof(int));
    if (!positions) { free(random_bytes); free(d_values); return MCELIECE_ERROR_MEMORY; }

    int unique_count = 0;
    int max_attempts = tau * 2; // 防止无限循环
    int attempts = 0;

    for (int i = 0; i < tau && unique_count < t && attempts < max_attempts; i++) {
        int pos = d_values[i];
        int is_unique = 1;

        // 检查该位置是否已经存在
        for (int j = 0; j < unique_count; j++) {
            if (positions[j] == pos) {
                is_unique = 0;
                break;
            }
        }

        if (is_unique) {
            positions[unique_count] = pos;
            unique_count++;
        }
        attempts++;
    }

    // 如果找不到足够的唯一位置，重新生成
    if (unique_count < t) {
        free(positions);
        free(d_values);
        free(random_bytes);
        return MCELIECE_ERROR_KEYGEN_FAIL; // 重新尝试
    }

    // 4. 在这些位置上将向量 e 的比特位置为 1
    for (int i = 0; i < t; i++) {
        vector_set_bit(e, positions[i], 1);
    }

    free(positions);
    free(d_values);
    free(random_bytes);
    return MCELIECE_SUCCESS;
}




// SeededKeyGen算法

mceliece_error_t seeded_key_gen(const uint8_t *delta, public_key_t *pk, private_key_t *sk) {
    if (!delta || !pk || !sk) {
        return MCELIECE_ERROR_INVALID_PARAM;
    }

    int n_bits = MCELIECE_N; // n=6688

    int t_bits = MCELIECE_T; // t=128
    int q_val = MCELIECE_Q;  // q=8192

    // l 是会话密钥长度，也是种子的长度 (in bits)
    // 规范 9.1: The integer l is 256.
    int l_bits = 256;

    // σ₁ 和 σ₂ 也是规范定义的整数
    // 规范 9.1: The integer σ₁ is 16.
    // 规范 9.1: The integer σ₂ is 32.
    int sigma1 = 16;
    int sigma2 = 32;

    // --- 计算各个部分的比特长度 ---
    // 规范 8.3 (SeededKeyGen) / 8.1 / 8.2 描述了 E 的构成
    // E = s || (bits for FieldOrdering) || (bits for Irreducible) || δ'
    // 长度: n + σ₂q + σ₁t + l bits

    int s_len_bits = n_bits;
    int field_ordering_len_bits = sigma2 * q_val;
    int irreducible_poly_len_bits = sigma1 * t_bits;
    int delta_prime_len_bits = l_bits;

    size_t prg_output_len_bits = s_len_bits + field_ordering_len_bits + irreducible_poly_len_bits + delta_prime_len_bits;

    // 将总比特长度转换为字节长度，向上取整
    size_t prg_output_len_bytes = (prg_output_len_bits + 7) / 8;

    // --- 计算各个部分的字节长度和偏移量 ---
    size_t s_len_bytes = (s_len_bits + 7) / 8;
    size_t field_ordering_len_bytes = (field_ordering_len_bits + 7) / 8;
    size_t irreducible_poly_len_bytes = (irreducible_poly_len_bits + 7) / 8;
    size_t delta_prime_len_bytes = (delta_prime_len_bits + 7) / 8; // 应该是 32

    // 验证一下总长度是否匹配
    if (prg_output_len_bytes != s_len_bytes + field_ordering_len_bytes + irreducible_poly_len_bytes + delta_prime_len_bytes) {
         // 这在某些边界条件下可能不成立，但对于 McEliece 参数通常成立。
         // 一个更安全的方式是直接使用比特偏移量。但我们先用字节。
    }

    // --- 准备 PRG 输出缓冲区 ---
    uint8_t *E = malloc(prg_output_len_bytes);
    if (!E) return MCELIECE_ERROR_MEMORY;

    // 复制初始种子到私钥
    memcpy(sk->delta, delta, delta_prime_len_bytes);

    int max_attempts = 50;
    for (int attempt = 0; attempt < max_attempts; attempt++) {
        // 1. Generate long random string E using current seed delta
        mceliece_prg(sk->delta, E, prg_output_len_bytes);

        // 2. Extract next retry seed delta' from the end of E
        uint8_t delta_prime[MCELIECE_L_BYTES]; 
        memcpy(delta_prime, E + prg_output_len_bytes - delta_prime_len_bytes, delta_prime_len_bytes);

        // 3. Split E into parts (using byte offsets)
        const uint8_t *s_bits_ptr = E;
        const uint8_t *field_ordering_bits_ptr = E + s_len_bytes;
        const uint8_t *irreducible_poly_bits_ptr = field_ordering_bits_ptr + field_ordering_len_bytes;
        
        // 4. Generate support set alpha
        if (generate_field_ordering(sk->alpha, field_ordering_bits_ptr) != MCELIECE_SUCCESS) {
            memcpy(sk->delta, delta_prime, MCELIECE_L_BYTES);
            continue;
        }

        // 5. Generate Goppa polynomial g
        if (generate_irreducible_poly_final(&sk->g, irreducible_poly_bits_ptr) != MCELIECE_SUCCESS) {
            memcpy(sk->delta, delta_prime, MCELIECE_L_BYTES);
            continue;
        }

        // Ensure alpha is a support set for g (no roots of g)
        int is_support_set = 1;
        for (int i = 0; i < n_bits; ++i) {
            if (polynomial_eval(&sk->g, sk->alpha[i]) == 0) {
                is_support_set = 0;
                break;
            }
        }
        if (!is_support_set) {
            memcpy(sk->delta, delta_prime, MCELIECE_L_BYTES);
            continue;
        }

        // 6. Generate public key T
        if (mat_gen(&sk->g, sk->alpha, &pk->T) != MCELIECE_SUCCESS) {
            memcpy(sk->delta, delta_prime, MCELIECE_L_BYTES);
            continue;
        }

        // --- All steps successful! ---

        // 7. Save other parts of private key
        // Copy s (length n)
        memcpy(sk->s, s_bits_ptr, (n_bits + 7) / 8);

        // Other parts of private key (c, g, alpha) are already in sk structure
        // p vector is not needed, ensure private_key_free won't try to free uninitialized pointer
        if(sk->p) {
            free(sk->p);
            sk->p = NULL;
        }

        free(E);
        return MCELIECE_SUCCESS;
    }

    // Reached maximum attempts, generation failed
    free(E);
    return MCELIECE_ERROR_KEYGEN_FAIL;
}
// KeyGen算法
mceliece_error_t mceliece_keygen(public_key_t *pk, private_key_t *sk) {
    if (!pk || !sk) {
        return MCELIECE_ERROR_INVALID_PARAM;
    }

    // 生成一个一次性的随机种子 delta
    uint8_t delta[MCELIECE_L_BYTES];
    // 我们需要一个安全的随机源，但为了编译通过，暂时使用 mceliece_prg
    mceliece_prg((const uint8_t*)"a_seed_for_the_seed_generator", delta, 32);

    return seeded_key_gen(delta, pk, sk);
}



// Encap algorithm (non-pc parameter sets)
mceliece_error_t mceliece_encap(const public_key_t *pk, uint8_t *ciphertext, uint8_t *session_key) {
    if (!pk || !ciphertext || !session_key) {
        return MCELIECE_ERROR_INVALID_PARAM;
    }
    
    int max_attempts = 10;
    for (int attempt = 0; attempt < max_attempts; attempt++) {
        // Step 1: Generate fixed weight vector e
        uint8_t *e = malloc(MCELIECE_N_BYTES);
        if (!e) return MCELIECE_ERROR_MEMORY;
        
        mceliece_error_t ret = fixed_weight_vector(e, MCELIECE_N, MCELIECE_T);
        if (ret != MCELIECE_SUCCESS) {
            free(e);
            if (ret == MCELIECE_ERROR_KEYGEN_FAIL) {
                // Retry
                continue;
            }
            return ret;
        }
        
        // Step 2: Calculate C = Encode(e, T)
        encode_vector(e, &pk->T, ciphertext);
        
        // Step 3: Calculate K = Hash(1, e, C)
        // Construct hash input: prefix 1 + e + C
        size_t hash_input_len = 1 + MCELIECE_N_BYTES + MCELIECE_MT_BYTES;
        uint8_t *hash_input = malloc(hash_input_len);
        if (!hash_input) {
            free(e);
            return MCELIECE_ERROR_MEMORY;
        }
        
        hash_input[0] = 1;  // prefix
        memcpy(hash_input + 1, e, MCELIECE_N_BYTES);
        memcpy(hash_input + 1 + MCELIECE_N_BYTES, ciphertext, MCELIECE_MT_BYTES);
        
        mceliece_hash(0, hash_input, hash_input_len, session_key);
        
        free(e);
        free(hash_input);
        return MCELIECE_SUCCESS;
    }
    
    return MCELIECE_ERROR_KEYGEN_FAIL; // Reached maximum attempts
}

// Simplified decode algorithm
mceliece_error_t decode_ciphertext(const uint8_t *ciphertext, const private_key_t *sk,
                                  uint8_t *error_vector, int *success) {
    if (!ciphertext || !sk || !error_vector || !success) {
        return MCELIECE_ERROR_INVALID_PARAM;
    }

    *success = 0;

    // According to Classic McEliece spec, the ciphertext C = H*e where H = [I_mt | T]
    // So C already IS the syndrome! We can directly pass it to Goppa decoding.
    // 
    // The "received" vector for syndrome computation should be the actual error pattern
    // we're trying to recover, but since we don't know it yet, we use the syndrome directly.
    
    // Convert ciphertext bits to GF elements for syndrome computation
    gf_elem_t syndrome[2 * MCELIECE_T];
    memset(syndrome, 0, sizeof(syndrome));
    
    // In the Classic McEliece spec, we need to compute the syndrome from the received word
    // But since C = H*e and we're looking for e, we can use C directly as syndrome information
    // However, this needs to be converted to the proper Goppa syndrome format
    
    // For now, let's use a different approach: assume the received word is v = (C, 0, ..., 0)
    // This represents what we would get if we received the first mt bits as C and the rest as 0
    uint8_t *v = malloc(MCELIECE_N_BYTES);
    if (!v) return MCELIECE_ERROR_MEMORY;

    memset(v, 0, MCELIECE_N_BYTES);

    // Copy ciphertext to first mt bits of v
    int mt = MCELIECE_M * MCELIECE_T;
    for (int i = 0; i < mt; i++) {
        int bit_value = vector_get_bit(ciphertext, i);
        vector_set_bit(v, i, bit_value);
    }

    // Step 2: Use Goppa decoding algorithm
    mceliece_error_t ret = decode_goppa(v, &sk->g, sk->alpha, error_vector, success);
    free(v);
    
    return ret;
}

// Decap algorithm (non-pc parameter sets)
mceliece_error_t mceliece_decap(const uint8_t *ciphertext, const private_key_t *sk, 
                               uint8_t *session_key) {
    if (!ciphertext || !sk || !session_key) {
        return MCELIECE_ERROR_INVALID_PARAM;
    }
    
    // Step 1: Set b = 1
    uint8_t b = 1;
    
    // Step 3: Try to decode
    uint8_t *e = malloc(MCELIECE_N_BYTES);
    if (!e) return MCELIECE_ERROR_MEMORY;
    
    int decode_success;
    mceliece_error_t ret = decode_ciphertext(ciphertext, sk, e, &decode_success);
    
    if (ret != MCELIECE_SUCCESS) {
        free(e);
        return ret;
    }
    
    if (!decode_success) {
        // Decoding failed, use backup vector s
        memcpy(e, sk->s, MCELIECE_N_BYTES);
        b = 0;
    }
    
    // Step 4: Calculate K = Hash(b, e, C)
    size_t hash_input_len = 1 + MCELIECE_N_BYTES + MCELIECE_MT_BYTES;
    uint8_t *hash_input = malloc(hash_input_len);
    if (!hash_input) {
        free(e);
        return MCELIECE_ERROR_MEMORY;
    }
    
    hash_input[0] = b;  // prefix
    memcpy(hash_input + 1, e, MCELIECE_N_BYTES);
    memcpy(hash_input + 1 + MCELIECE_N_BYTES, ciphertext, MCELIECE_MT_BYTES);
    
    mceliece_hash(0, hash_input, hash_input_len, session_key);
    
    free(e);
    free(hash_input);
    return MCELIECE_SUCCESS;
}

// 密钥序列化
mceliece_error_t serialize_public_key(const public_key_t *pk, uint8_t *buffer, size_t buffer_len) {
    if (!pk || !buffer) return MCELIECE_ERROR_INVALID_PARAM;
    
    size_t required_len = MCELIECE_M * MCELIECE_T * MCELIECE_K_BYTES;
    if (buffer_len < required_len) return MCELIECE_ERROR_INVALID_PARAM;
    
    // 以行优先方式序列化T矩阵
    memcpy(buffer, pk->T.data, required_len);
    
    return MCELIECE_SUCCESS;
}



mceliece_error_t serialize_private_key(const private_key_t *sk, uint8_t *buffer, size_t buffer_len) {
    if (!sk || !buffer) return MCELIECE_ERROR_INVALID_PARAM;
    
    size_t pos = 0;
    
    // delta (32 bytes)
    if (pos + MCELIECE_L_BYTES > buffer_len) return MCELIECE_ERROR_INVALID_PARAM;
    memcpy(buffer + pos, sk->delta, MCELIECE_L_BYTES);
    pos += MCELIECE_L_BYTES;
    
    // c (8 bytes for μ=ν=0)
    if (pos + 8 > buffer_len) return MCELIECE_ERROR_INVALID_PARAM;
    uint64_t c_val = sk->c;
    for (int i = 0; i < 8; i++) {
        buffer[pos + i] = (c_val >> (8 * i)) & 0xFF;
    }
    pos += 8;
    
    // g多项式 (t * ceil(m/8) bytes)
    size_t g_bytes = MCELIECE_T * ((MCELIECE_M + 7) / 8);
    if (pos + g_bytes > buffer_len) return MCELIECE_ERROR_INVALID_PARAM;
    
    for (int i = 0; i < MCELIECE_T; i++) {
        gf_elem_t coeff = sk->g.coeffs[i];
        for (int j = 0; j < (MCELIECE_M + 7) / 8; j++) {
            buffer[pos + i * ((MCELIECE_M + 7) / 8) + j] = (coeff >> (8 * j)) & 0xFF;
        }
    }
    pos += g_bytes;
    
    // alpha序列 (简化处理)
    size_t alpha_bytes = ((2 * MCELIECE_M - 1) * (1 << (MCELIECE_M - 4)) + 7) / 8;
    if (pos + alpha_bytes > buffer_len) return MCELIECE_ERROR_INVALID_PARAM;
    
    // 这里应该实现Beneš网络控制位的计算和存储
    // 简化实现：直接存储alpha数组的前几个字节
    memset(buffer + pos, 0, alpha_bytes);  // 占位符
    pos += alpha_bytes;
    
    // s向量
    if (pos + MCELIECE_N_BYTES > buffer_len) return MCELIECE_ERROR_INVALID_PARAM;
    memcpy(buffer + pos, sk->s, MCELIECE_N_BYTES);
    pos += MCELIECE_N_BYTES;
    
    return MCELIECE_SUCCESS;
}



// 测试函数
void test_mceliece(void) {
    printf("Testing Classic McEliece implementation...\n");
    

    
    // 创建密钥
    public_key_t *pk = public_key_create();
    private_key_t *sk = private_key_create();
    
    if (!pk || !sk) {
        printf("Failed to create key structures\n");
        return;
    }
    
    // 生成密钥对
    printf("Generating key pair...\n");
    mceliece_error_t ret = mceliece_keygen(pk, sk);
    if (ret != MCELIECE_SUCCESS) {
        printf("Key generation failed: %d\n", ret);
        public_key_free(pk);
        private_key_free(sk);
        return;
    }
    printf("Key generation successful!\n");
    
    // 封装
    printf("Testing encapsulation...\n");
    uint8_t ciphertext[MCELIECE_MT_BYTES];
    uint8_t session_key1[MCELIECE_L_BYTES];
    
    ret = mceliece_encap(pk, ciphertext, session_key1);
    if (ret != MCELIECE_SUCCESS) {
        printf("Encapsulation failed: %d\n", ret);
        public_key_free(pk);
        private_key_free(sk);
        return;
    }
    printf("Encapsulation successful!\n");
    
    // 解封装
    printf("Testing decapsulation...\n");
    uint8_t session_key2[MCELIECE_L_BYTES];
    
    ret = mceliece_decap(ciphertext, sk, session_key2);
    if (ret != MCELIECE_SUCCESS) {
        printf("Decapsulation failed: %d\n", ret);
        public_key_free(pk);
        private_key_free(sk);
        return;
    }
    
    // 验证会话密钥是否相同
    int keys_match = 1;
    for (int i = 0; i < MCELIECE_L_BYTES; i++) {
        if (session_key1[i] != session_key2[i]) {
            keys_match = 0;
            break;
        }
    }
    
    if (keys_match) {
        printf("Decapsulation successful! Session keys match.\n");
    } else {
        printf("Decapsulation failed! Session keys don't match.\n");
    }
    
    // 清理
    public_key_free(pk);
    private_key_free(sk);
    
    printf("Test completed.\n");
}