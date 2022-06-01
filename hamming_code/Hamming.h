#pragma once

#ifndef HAMMING_CODE
#define HAMMING_CODE

#include "../Mat.h"
#include <cmath>
#include <iostream>

class Hamming : public Mat<int> {
private:
    int P, R, N;
    Mat<int> h_parity_check_matrix, h_generator_matrix, left_side_parity_check_matrix, right_side_parity_check_matrix;

    bool is_collinear(const Mat<int> &vector_1, const Mat<int> &vector_2, int len_of_vector, int mod);

    bool is_vector_collinear_to_column(const Mat<int> &vector, Mat<int> matrix, int index, int mod);

    bool
    is_vector_collinear_to_columns_from_range(const Mat<int> &vector, const Mat<int> &matrix, int begin_id, int end_id,
                                              int mod);

    bool is_coefficient_valid(int coordinate_in_vector_1, int coordinate_in_vector_2, double coefficient,
                              bool coefficient_was_made_of_div_1_on_2, int mod);

    Mat<int> add_one_to_vector(Mat<int> vector, int vector_size, int mod);

public:
    Hamming(int p, int r);

    Mat<int> Party_check_matrix();

    Mat<int> Generator_matrix();

    bool syndrome_is_zero(Mat<int> syndrome);

    int find_index_of_mistake(Mat<int> parityCheckMatrix, const Mat<int> &syndrome);

    Mat<int> decode(Mat<int> &code_with_effords);

    Mat<int> find_syndrome(Mat<int> &parity_check_matrix, Mat<int> vector);
};

#endif //HAMMING_CODE

Hamming::Hamming(int p, int r) : Mat(p, r) {
    P = p;
    R = r;
    N = (int) ((pow(P, R) - 1) / (P - 1));

    h_generator_matrix = zeros(N - R, N);
    h_parity_check_matrix = zeros(R, N);
    left_side_parity_check_matrix = zeros(R, N - R);
    right_side_parity_check_matrix = zeros(N - R, R);

    int rows = R, cols = N;
    Mat<int> identity_matrix = eye(rows);

    for (int j = rows - cols; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            h_parity_check_matrix(i, j) = identity_matrix(i, j - (cols - rows));
        }
    }

    Mat<int> vector = Mat<int>(1, rows);
    vector(0, rows - 1) = 1;

    for (int j = cols - rows - 1; j >= 0; j--) {
        while (is_vector_collinear_to_columns_from_range(vector, h_parity_check_matrix, j + 1, cols, P)) {
            vector = add_one_to_vector(vector, R, P);
        }
        for (int i = 0; i < rows; i++) {
            h_parity_check_matrix(i, j) = vector(0, i);
        }
    }

    int K = N - R;
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < N; j++) {
            left_side_parity_check_matrix(i, j) = h_parity_check_matrix(i, j);
        }
    }

    right_side_parity_check_matrix = left_side_parity_check_matrix.transposition() * -1;

    Mat<int> identity_matrix_2 = eye(K);
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            h_generator_matrix(i, j) = identity_matrix_2(i, j);
        }
        for (int j = K; j < N; j++) {
            h_generator_matrix(i, j) = right_side_parity_check_matrix(i, j - K);
        }
    }
}

bool Hamming::is_coefficient_valid(int coordinate_in_vector_1, int coordinate_in_vector_2, double coefficient,
                                   bool coefficient_was_made_of_div_1_on_2, int mod) {
    if (coefficient_was_made_of_div_1_on_2) {
        if (((int) (coordinate_in_vector_1 * coefficient) % mod == coordinate_in_vector_2) &&
            (coordinate_in_vector_1 * coefficient - (int) coordinate_in_vector_1 * coefficient == 0.0)) {
            return true;
        }
        return false;
    } else {
        if (((int) (coordinate_in_vector_2 * coefficient) % mod == coordinate_in_vector_1) &&
            (coordinate_in_vector_2 * coefficient - (int) coordinate_in_vector_2 * coefficient == 0.0)) {
            return true;
        }
        return false;
    }
}

bool Hamming::is_collinear(const Mat<int> &vector_1, const Mat<int> &vector_2, int len_of_vector, int mod) {
    double coefficient = NULL;
    bool coefficient_was_made_of_div_1_on_2 = NULL;

    for (int i = 0; i < len_of_vector; i++) {
        if ((vector_1(0, i) + vector_2(0, i) != 0) && (vector_1(0, i) * vector_2(0, i) == 0)) {
            return false;
        } else {
            if (vector_1(0, i) == vector_2(0, i) && vector_1(0, i) == 0) {
                continue;
            } else {
                if (coefficient == NULL) {
                    if (vector_2(0, i) > vector_1(0, i)) {
                        coefficient = (double) vector_2(0, i) / vector_1(0, i);
                        coefficient_was_made_of_div_1_on_2 = true;
                    } else {
                        coefficient = (double) vector_1(0, i) / vector_2(0, i);
                        coefficient_was_made_of_div_1_on_2 = false;
                    }
                } else {
                    if (!is_coefficient_valid(vector_1(0, i), vector_2(0, i), coefficient,
                                              coefficient_was_made_of_div_1_on_2, mod)) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

bool Hamming::is_vector_collinear_to_column(const Mat<int> &vector, Mat<int> matrix, int index, int mod) {
    int len_of_vector = matrix.get_num_rows();
    Mat<int> remake_column = Mat<int>(1, len_of_vector);

    for (int i = 0; i < len_of_vector; i++) {
        remake_column(0, i) = matrix(i, index);
    }
    if (is_collinear(vector, remake_column, len_of_vector, mod)) {
        return true;
    }
    return false;
}

bool Hamming::is_vector_collinear_to_columns_from_range(const Mat<int> &vector, const Mat<int> &matrix, int begin_id,
                                                        int end_id, int mod) {
    for (int index = begin_id; index < end_id; index++) {
        if (is_vector_collinear_to_column(vector, matrix, index, mod)) {
            return true;
        }
    }
    return false;
}

Mat<int> Hamming::add_one_to_vector(Mat<int> vector, int vector_size, int mod) {
    for (int i = vector_size - 1; i >= 0; i--) {
        if (vector(0, i) + 1 != mod) {
            vector(0, i) = vector(0, i) + 1;
            return vector;
        } else {
            if (i == 0) {
                std::cout << "Could not increase vector!" << std::endl;
                exit(-1);
            }
            vector(0, i) = 0;
        }
    }
    return vector;
}

Mat<int> Hamming::Party_check_matrix() {
    return h_parity_check_matrix;
}

Mat<int> Hamming::Generator_matrix() {
    return h_generator_matrix;
}

Mat<int> Hamming::find_syndrome(Mat<int> &parity_check_matrix, Mat<int> vector){
    Mat<int> transposed_vector, syndrome;
    transposed_vector = vector.transposition();
    syndrome = (parity_check_matrix * transposed_vector);
    syndrome.matrix_by_mod(P);
    return syndrome;
}

bool Hamming::syndrome_is_zero(Mat<int> syndrome) {
    int rows = syndrome.get_num_rows();
    for (int i = 0; i < rows; ++i) {
        if (syndrome(i, 0) != 0) {
            return false;
        }
    }
    return true;
}


int Hamming::find_index_of_mistake(Mat<int> parityCheckMatrix, const Mat<int> &syndrome) {
    int cols = h_parity_check_matrix.get_num_cols();
    int rows = h_parity_check_matrix.get_num_rows();

    for (int columnIdx = 0; columnIdx < cols; ++columnIdx) {
        bool syndrome_is_zero = true;

        for (int rowIdx = 0; rowIdx < rows; ++rowIdx) {
            if (syndrome(rowIdx, 0) != parityCheckMatrix(rowIdx, columnIdx)) {
                syndrome_is_zero = false;
                break;
            }
        }
        if (syndrome_is_zero) {
            return columnIdx;
        }
    }

    return -1;
}

Mat<int> Hamming::decode(Mat<int> &code_with_effords) {
    Mat<int> syndrome = find_syndrome(h_parity_check_matrix, code_with_effords);
    if (syndrome_is_zero(syndrome)) {
        return code_with_effords;
    }

    int mistake_id = find_index_of_mistake(h_parity_check_matrix, syndrome), mistake_id_value = -1;
    Mat<int> code_without_mistakes;
    code_without_mistakes = code_with_effords;
    code_without_mistakes(0, mistake_id) = mistake_id_value;

    do {
        ++mistake_id_value;
        code_without_mistakes(0, mistake_id) = mistake_id_value;
        syndrome = find_syndrome(h_parity_check_matrix, code_without_mistakes);
    } while (!syndrome_is_zero(syndrome));

    return code_without_mistakes;
}

