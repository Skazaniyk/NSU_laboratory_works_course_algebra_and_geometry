#include <stdio.h>

#define SIZE_LIMIT 200

double determinant_of_matrix(double matrix[SIZE_LIMIT][SIZE_LIMIT], int matrix_size) {
    int non_zero_row = 0, current_row = 0;

    for (int current_column = 0; current_column < matrix_size; ++current_column) {

        int column_consists_of_zeros = 1;
        for (int row = current_row; row < matrix_size; ++row) {
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
            for (int column = 0; column < matrix_size; ++column) {
                matrix[current_column][column] += matrix[non_zero_row][column];
            }
        }

        for (int row = current_row + 1; row < matrix_size; ++row) {
            double first_element_of_current_row = matrix[row][current_column];
            for (int column = 0; column < matrix_size; ++column) {
                matrix[row][column] +=
                        (-(first_element_of_current_row / non_zero_element)) * matrix[current_row][column];
            }
        }
        ++current_row;
    }

    double determinant = 1;
    for (int diagonal_index = 0; diagonal_index < matrix_size; ++diagonal_index) {
        determinant *= matrix[diagonal_index][diagonal_index];
    }
    return determinant;
}

int main() {
    int matrix_size;
    double matrix[SIZE_LIMIT][SIZE_LIMIT] = {0};

    scanf("%d", &matrix_size);
    for (int row = 0; row < matrix_size; ++row) {
        for (int column = 0; column < matrix_size; ++column) {
            scanf("%lf", &matrix[row][column]);
        }
    }

    printf("%0.2lf", determinant_of_matrix(matrix, matrix_size));

    return 0;
}