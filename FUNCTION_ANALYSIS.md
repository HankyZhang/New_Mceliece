# Function Analysis and Cleanup Report

## Overview
This document provides a comprehensive analysis of all functions in the Classic McEliece implementation, categorizing them by usage patterns and identifying opportunities for cleanup.

## Function Categories

### 1. Core Used Functions (Keep)
These functions are essential and actively used in the implementation:

#### Main Entry Points
- `main()` - Entry point
- `print_usage()` - Command-line help
- `print_parameters()` - Parameter display
- `demo_keygen()` - Key generation demo
- `demo_complete()` - Complete demo
- `benchmark()` - Performance testing

#### Core KEM Operations
- `mceliece_keygen()` - Key generation
- `mceliece_encap()` - Encapsulation
- `mceliece_decap()` - Decapsulation
- `seeded_key_gen()` - Seeded key generation
- `fixed_weight_vector()` - Error vector generation
- `decode_ciphertext()` - Ciphertext decoding

#### Finite Field Operations
- `gf_init()` - GF table initialization
- `gf_add()` - GF addition
- `gf_mul()` - GF multiplication (optimized)
- `gf_mul_for_init()` - GF multiplication for initialization (static)
- `gf_inv()` - GF inverse
- `gf_div()` - GF division
- `gf_pow()` - GF power
- `bits_to_gf()` - Convert bits to GF element
- `gf_to_bits()` - Convert GF element to bits

#### Polynomial Operations (Used)
- `polynomial_eval()` - Polynomial evaluation (heavily used)
- `polynomial_set_coeff()` - Set polynomial coefficient
- `polynomial_copy()` - Copy polynomial
- `polynomial_is_zero()` - Check if polynomial is zero
- `polynomial_add()` - Polynomial addition
- `polynomial_mul()` - Polynomial multiplication
- `polynomial_div()` - Polynomial division

#### Matrix Operations (Used)
- `matrix_create()` - Create matrix
- `matrix_free()` - Free matrix
- `matrix_set_bit()` - Set matrix bit
- `matrix_get_bit()` - Get matrix bit
- `matrix_swap_rows()` - Swap matrix rows
- `matrix_swap_cols()` - Swap matrix columns
- `matrix_xor_rows()` - XOR matrix rows
- `matrix_is_systematic()` - Check systematic form
- `reduce_to_systematic_form()` - Convert to systematic form
- `mat_gen()` - Generate matrix from Goppa polynomial
- `encode_vector()` - Encode error vector
- `matrix_vector_multiply()` - Matrix-vector multiplication
- `compute_syndrome()` - Compute syndrome

#### Decoding Algorithms
- `berlekamp_massey()` - Berlekamp-Massey algorithm
- `chien_search()` - Chien search
- `decode_goppa()` - Goppa code decoding

#### Random Generation
- `generate_field_ordering()` - Generate field ordering
- `generate_irreducible_poly_final()` - Generate irreducible polynomial

#### Cryptographic Primitives (Used internally)
- `shake256()` - SHAKE256 one-shot function
- `mceliece_hash()` - McEliece hash function
- `mceliece_prg()` - Pseudorandom generator

#### Memory Management
- `polynomial_create()` - Create polynomial
- `polynomial_free()` - Free polynomial
- `public_key_create()` - Create public key
- `public_key_free()` - Free public key
- `private_key_create()` - Create private key
- `private_key_free()` - Free private key

#### Utility Functions (Used)
- `vector_set_bit()` - Set vector bit
- `vector_get_bit()` - Get vector bit
- `vector_weight()` - Calculate vector weight

#### Serialization (Used)
- `serialize_public_key()` - Serialize public key
- `serialize_private_key()` - Serialize private key

#### Test Functions (Used)
- `test_mceliece()` - Basic test
- `run_all_tests()` - Comprehensive test suite
- All individual test functions (test_gf_operations, etc.)

### 2. UNUSED FUNCTIONS (Can be removed)

#### Unused Serialization
- **`deserialize_public_key()`** - Public key deserialization (never called)
- **`deserialize_private_key()`** - Private key deserialization (never called)

#### Unused Polynomial Operations
- **`polynomial_gcd()`** - Polynomial GCD calculation (never called)
- **`polynomial_pow_mod()`** - Polynomial power modulo (never called)
- **`polynomial_eval_at_poly()`** - Polynomial evaluation at polynomial (never called)

#### Unused Matrix Operations
- **`matrix_find_pivot()`** - Find matrix pivot (never called)
- **`print_matrix()`** - Debug matrix printing (never called)
- **`matrix_copy()`** - Matrix copying (never called)
- **`construct_parity_check_matrix()`** - Construct H matrix (never called)

#### Unused Finite Field Operations
- **`matrix_fq_create()`** - Create GF matrix (never called)
- **`matrix_fq_free()`** - Free GF matrix (never called)
- **`solve_linear_system()`** - Solve linear system over GF (never called)

#### Unused Decoding Functions
- **`forney_algorithm()`** - Forney algorithm for error values (never called, binary codes don't need error values)

#### Unused Utility Functions
- **`print_bytes()`** - Debug byte printing (never called)
- **`xor_bytes()`** - Byte array XOR (never called)
- **`vector_clear_bit()`** - Clear vector bit (never called, use vector_set_bit with 0)

#### Unused Random Generation
- **`is_irreducible_simple()`** - Simple irreducibility check (only used by generate_irreducible_poly_final internally)
- **`set_random_seed()`** - Set random seed (never called)
- **`init_random()`** - Initialize random (never called)
- **`random_byte()`** - Generate random byte (never called)
- **`random_bytes()`** - Generate random bytes (never called)
- **`compare_pairs()`** - Compare function for sorting (used by qsort but could be inlined)

#### Unused Cryptographic Internal Functions
- **`shake256_init()`** - SHAKE256 initialization (used internally by shake256 and mceliece_*)
- **`shake256_absorb()`** - SHAKE256 absorb (used internally)
- **`shake256_finalize()`** - SHAKE256 finalize (used internally)
- **`shake256_squeeze()`** - SHAKE256 squeeze (used internally)
- **`rotl64()`** - 64-bit rotate left (used internally by keccak_f1600)
- **`keccak_f1600()`** - Keccak permutation (used internally by SHAKE256)

### 3. Internal/Static Functions (Keep but could be static)
These functions are only used internally within their modules:
- `gf_mul_for_init()` - Already static
- `is_irreducible_simple()` - Could be static
- `compare_pairs()` - Could be static
- `rotl64()` - Already static
- `keccak_f1600()` - Already static
- SHAKE256 internal functions - Used internally by public SHAKE256 functions

## Recommendations

### Phase 1: Remove Clearly Unused Functions
Remove these functions as they are never called:
1. `deserialize_public_key()`
2. `deserialize_private_key()`
3. `polynomial_gcd()`
4. `polynomial_pow_mod()`
5. `polynomial_eval_at_poly()`
6. `matrix_find_pivot()`
7. `print_matrix()`
8. `matrix_copy()`
9. `construct_parity_check_matrix()`
10. `matrix_fq_create()`
11. `matrix_fq_free()`
12. `solve_linear_system()`
13. `forney_algorithm()`
14. `print_bytes()`
15. `xor_bytes()`
16. `vector_clear_bit()`
17. `set_random_seed()`
18. `init_random()`
19. `random_byte()`
20. `random_bytes()`

### Phase 2: Make Internal Functions Static
Make these functions static to their compilation units:
1. `is_irreducible_simple()` in mceliece_random.c
2. `compare_pairs()` in mceliece_random.c

### Phase 3: Keep Essential Internal Functions
Keep SHAKE256 and Keccak internal functions as they are essential for the cryptographic implementation.

## Impact Analysis
Removing these unused functions will:
- Reduce binary size
- Improve code maintainability
- Eliminate dead code
- Reduce compilation time
- Make the codebase cleaner and easier to understand

The relationship between `gf_mul` and `gf_mul_for_init` is crucial:
- `gf_mul_for_init`: Bootstrap function for table generation during GF initialization
- `gf_mul`: Optimized runtime function using pre-computed lookup tables
- Both are essential for the GF arithmetic system
