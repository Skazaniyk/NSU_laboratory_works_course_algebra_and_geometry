cmake_minimum_required(VERSION 3.20)
project(hamming_code)

set(CMAKE_CXX_STANDARD 17)

add_executable(Matrix
        Mat.h)

add_executable(Hamming
        Hamming.h)

add_executable(Test
        test.cpp)