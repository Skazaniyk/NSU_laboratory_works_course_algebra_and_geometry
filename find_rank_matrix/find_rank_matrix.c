#include <stdio.h>

#define SIZE_LIMIT 200

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

int main() {
    int matrix_size_row, matrix_size_column;
    double matrix[SIZE_LIMIT][SIZE_LIMIT] = {0};

    scanf("%d%d", &matrix_size_row, &matrix_size_column);
    for (int row = 0; row < matrix_size_row; ++row) {
        for (int column = 0; column < matrix_size_column; ++column) {
            scanf("%lf", &matrix[row][column]);
        }
    }

    int rank = matrix_rank(matrix, matrix_size_column, matrix_size_row);
    printf("%d ", rank);

    if (rank == matrix_size_row) {
        printf("true");
    } else {
        printf("false");
    }

    return 0;
}