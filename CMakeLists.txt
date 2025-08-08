cmake_minimum_required(VERSION 4.0)
project(untitled1)

set(CMAKE_CXX_STANDARD 20)

add_executable(untitled1 main.c
        mceliece_types.h
        mceliece_utils.c
        mceliece_gf.c
        mceliece_shake.c
        mceliece_random.c
        gf.h
        vector.h
        matrix.h
        mceliece_matrix_ops.c
        mceliece_decode.h
        mceliece.h
        mceliece_random.h
        test_gauss.c)
