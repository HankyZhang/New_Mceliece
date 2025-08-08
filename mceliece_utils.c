#include "mceliece_types.h"
#include <stdio.h>
#include <assert.h>

// 多项式创建
polynomial_t* polynomial_create(int max_degree) {
    polynomial_t *poly = malloc(sizeof(polynomial_t));
    if (!poly) return NULL;
    
    poly->coeffs = calloc(max_degree + 1, sizeof(gf_elem_t));
    if (!poly->coeffs) {
        free(poly);
        return NULL;
    }
    
    poly->degree = -1;  // 表示零多项式
    poly->max_degree = max_degree;
    return poly;
}

// 多项式释放
void polynomial_free(polynomial_t *poly) {
    if (poly) {
        if (poly->coeffs) free(poly->coeffs);
        free(poly);
    }
}

// 矩阵创建
matrix_t* matrix_create(int rows, int cols) {
    matrix_t *mat = malloc(sizeof(matrix_t));
    if (!mat) return NULL;
    
    mat->rows = rows;
    mat->cols = cols;
    mat->cols_bytes = (cols + 7) / 8;  // 按字节对齐
    
    mat->data = calloc(rows * mat->cols_bytes, sizeof(uint8_t));
    if (!mat->data) {
        free(mat);
        return NULL;
    }
    
    return mat;
}

// 矩阵释放
void matrix_free(matrix_t *mat) {
    if (mat) {
        if (mat->data) free(mat->data);
        free(mat);
    }
}

// 私钥创建
private_key_t* private_key_create(void) {
    private_key_t *sk = malloc(sizeof(private_key_t));
    sk->p = NULL;
    if (!sk) return NULL;
    
    memset(sk, 0, sizeof(private_key_t));
    
    // 初始化Goppa多项式
    polynomial_t *g = polynomial_create(MCELIECE_T);
    if (!g) {
        free(sk);
        return NULL;
    }
    sk->g = *g;
    free(g);  // 只释放结构体，不释放coeffs
    
    // 分配alpha数组
    sk->alpha = calloc(MCELIECE_Q, sizeof(gf_elem_t));
    if (!sk->alpha) {
        free(sk->g.coeffs);
        free(sk);
        return NULL;
    }
    
    // 设置默认c值（对于μ=ν=0）
    sk->c = (1ULL << 32) - 1;
    
    return sk;
}

// 私钥释放
void private_key_free(private_key_t *sk) {
    if (sk->p) free(sk->p);
    if (sk) {
        if (sk->g.coeffs) free(sk->g.coeffs);
        if (sk->alpha) free(sk->alpha);
        free(sk);
    }
}

// 公钥创建
public_key_t* public_key_create(void) {
    public_key_t *pk = malloc(sizeof(public_key_t));
    if (!pk) return NULL;
    
    matrix_t *T = matrix_create(MCELIECE_M * MCELIECE_T, MCELIECE_K);
    if (!T) {
        free(pk);
        return NULL;
    }
    
    pk->T = *T;
    free(T);  // 只释放结构体，不释放data
    
    return pk;
}

// 公钥释放
void public_key_free(public_key_t *pk) {
    if (pk) {
        if (pk->T.data) free(pk->T.data);
        free(pk);
    }
}

// 打印字节数组（调试用）
void print_bytes(const char *label, const uint8_t *data, int len) {
    printf("%s: ", label);
    for (int i = 0; i < len; i++) {
        printf("%02x", data[i]);
        if (i < len - 1) printf(" ");
    }
    printf("\n");
}

// 字节异或
void xor_bytes(uint8_t *dst, const uint8_t *src1, const uint8_t *src2, int len) {
    for (int i = 0; i < len; i++) {
        dst[i] = src1[i] ^ src2[i];
    }
}

// 矩阵位设置
void vector_set_bit(uint8_t *vec, int bit_idx, int value) {
    int byte_idx = bit_idx / 8;
    int bit_pos = bit_idx % 8;
    if (value) {
        vec[byte_idx] |= (1 << bit_pos);
    } else {
        vec[byte_idx] &= ~(1 << bit_pos);
    }
}

// 矩阵位获取
int matrix_get_bit(const matrix_t *mat, int row, int col) {
    assert(row < mat->rows && col < mat->cols);
    int byte_idx = row * mat->cols_bytes + (col / 8);
    int bit_idx = col % 8;
    
    return (mat->data[byte_idx] >> bit_idx) & 1;
}

// 向量权重计算
int vector_weight(const uint8_t *vec, int len_bytes) {
    int weight = 0;
    for (int i = 0; i < len_bytes; i++) {
        uint8_t byte = vec[i];
        // 计算字节中1的个数（Brian Kernighan算法）
        while (byte) {
            byte &= byte - 1;
            weight++;
        }
    }
    return weight;
}

void matrix_set_bit(matrix_t *mat, int row, int col, int bit) {
    assert(row < mat->rows && col < mat->cols);
    int byte_idx = row * mat->cols_bytes + (col / 8);
    int bit_idx = col % 8;

    if (bit) {
        mat->data[byte_idx] |= (1 << bit_idx);
    } else {
        mat->data[byte_idx] &= ~(1 << bit_idx);
    }
}

// 向量位清除
void vector_clear_bit(uint8_t *vec, int bit_idx) {
    int byte_idx = bit_idx / 8;
    int bit_pos = bit_idx % 8;
    vec[byte_idx] &= ~(1 << bit_pos);
}

// 向量位获取
int vector_get_bit(const uint8_t *vec, int bit_idx) {
    int byte_idx = bit_idx / 8;
    int bit_pos = bit_idx % 8;
    return (vec[byte_idx] >> bit_pos) & 1;
}