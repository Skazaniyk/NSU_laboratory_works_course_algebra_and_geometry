#include <stdio.h>

#define SIZE_LIMIT 3

int matrix_rank(double matrix[SIZE_LIMIT][SIZE_LIMIT], int martix_size_column, int matrix_size_row) {
    int non_zero_row = 0, current_row = 0;

    for (int current_column = 0; current_column < martix_size_column; ++current_column) {

        int column_consists_of_zeros = 1;
        for (int row = current_row; row < matrix_size_row; ++row) {
            if (matrix[row][current_column] != 0) {
                column_consists_of_zeros = 0;
                non_zero_row = row;
                break;
            }
        }
        if (column_consists_of_zeros == 1) {
            continue;
        }

        double non_zero_element = matrix[non_zero_row][current_column];
        if (non_zero_row != current_row) {
            for (int column = 0; column < martix_size_column; ++column) {
                matrix[current_column][column] += matrix[non_zero_row][column];
            }
        }

        for (int row = current_row + 1; row < matrix_size_row; ++row) {
            double first_element_of_current_row = matrix[row][current_column];
            for (int column = 0; column < martix_size_column; ++column) {
                matrix[row][column] +=
                        (-(first_element_of_current_row / non_zero_element)) * matrix[current_row][column];
            }
        }
        ++current_row;
    }

    int rank = 0;
    for (int row = 0; row < matrix_size_row; ++row) {
        for (int column = 0; column < martix_size_column; ++column) {
            if (matrix[row][column] != 0) {
                ++rank;
                break;
            }
        }
    }

    return rank;
}

void vector(double matrix[SIZE_LIMIT][SIZE_LIMIT]) {
    int non_zero_row = 0, current_row = 0, martix_size_column = SIZE_LIMIT, matrix_size_row = SIZE_LIMIT;
    double duplicat_matrix[SIZE_LIMIT][SIZE_LIMIT] = {0}, duplicat_1_2_matrix[SIZE_LIMIT - 1][SIZE_LIMIT] = {0},
            duplicat_1_3_matrix[SIZE_LIMIT - 1][SIZE_LIMIT] = {0}, duplicat_2_3_matrix[
            SIZE_LIMIT - 1][SIZE_LIMIT] = {0};

    for (int row = 0; row < matrix_size_row; ++row) {
        for (int column = 0; column < martix_size_column; ++column) {
            duplicat_matrix[row][column] = matrix[row][column];
            switch (row) {
                case 0:
                    duplicat_1_2_matrix[row][column] = matrix[row][column];
                    duplicat_1_3_matrix[row][column] = matrix[row][column];
                    break;

                case 1:
                    duplicat_1_2_matrix[row][column] = matrix[row][column];
                    duplicat_2_3_matrix[row - 1][column] = matrix[row][column];
                    break;

                case 2:
                    duplicat_2_3_matrix[row - 1][column] = matrix[row][column];
                    duplicat_1_3_matrix[row - 1][column] = matrix[row][column];
                    break;
            }
        }
    }

    int rank_all = matrix_rank(duplicat_matrix, 3, 3);

    if (rank_all == 3) {
        printf("false");
        return;
    }
    if (rank_all == 1 || rank_all == 0) {
        printf("true\n");
        for (int row = 0; row < matrix_size_row; ++row) {
            for (int column = 0; column < martix_size_column; ++column) {
                printf("%0.2lf ", matrix[row][column]);
            }
            printf("\n");
        }
        return;
    }

    int rank_1_2 = matrix_rank(duplicat_1_2_matrix, 3, 2),
            rank_1_3 = matrix_rank(duplicat_1_3_matrix, 3, 2),
            rank_2_3 = matrix_rank(duplicat_2_3_matrix, 3, 2);
    if (rank_1_2 == 1) {
        printf("true\n");
        for (int row = 0; row < matrix_size_row; ++row) {
            if (row == 2) {
                continue;
            }
            for (int column = 0; column < martix_size_column; ++column) {
                printf("%0.2lf ", matrix[row][column]);
            }
            printf("\n");
        }
        return;
    }
    if (rank_1_3 == 1) {
        printf("true\n");
        for (int row = 0; row < matrix_size_row; ++row) {
            if (row == 1) {
                continue;
            }
            for (int column = 0; column < martix_size_column; ++column) {
                printf("%0.2lf ", matrix[row][column]);
            }
            printf("\n");
        }
        return;
    }
    if (rank_2_3 == 1) {
        printf("true\n");
        for (int row = 0; row < matrix_size_row; ++row) {
            if (row == 0) {
                continue;
            }
            for (int column = 0; column < martix_size_column; ++column) {
                printf("%0.2lf ", matrix[row][column]);
            }
            printf("\n");
        }
        return;
    }
}

int main() {
    int matrix_size_row = SIZE_LIMIT, matrix_size_column = SIZE_LIMIT;
    double matrix[SIZE_LIMIT][SIZE_LIMIT] = {0};

    for (int row = 0; row < matrix_size_row; ++row) {
        for (int column = 0; column < matrix_size_column; ++column) {
            scanf("%lf", &matrix[row][column]);
        }
    }

    vector(matrix);

    return 0;
}