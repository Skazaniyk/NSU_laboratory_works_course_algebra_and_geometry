#include <stdio.h>

#define SIZE_LIMIT 200
double identity_matrix[SIZE_LIMIT][SIZE_LIMIT] = {0};

double inverse_matrix(double matrix[SIZE_LIMIT][SIZE_LIMIT], int matrix_size) {
    int non_zero_row = 0, current_row = 0;
    for (int diagonal_index = 0; diagonal_index < matrix_size; ++diagonal_index) {
        identity_matrix[diagonal_index][diagonal_index] = 1;
    }

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
                identity_matrix[current_column][column] += identity_matrix[non_zero_row][column];
            }
        }

        for (int row = current_row + 1; row < matrix_size; ++row) {
            double first_element_of_current_row = matrix[row][current_column];
            double alpha_for_transvection = (-(first_element_of_current_row / non_zero_element));
            for (int column = 0; column < matrix_size; ++column) {
                matrix[row][column] += alpha_for_transvection * matrix[current_row][column];
                identity_matrix[row][column] += alpha_for_transvection * identity_matrix[current_row][column];
            }
        }
        ++current_row;
    }

    double determinant = 1;
    for (int diagonal_index = 0; diagonal_index < matrix_size; ++diagonal_index) {
        determinant *= matrix[diagonal_index][diagonal_index];
    }

    if (determinant == 0) {
        return 0;
    } else {
        for (int column = matrix_size - 1; column >= 0; --column) {
            int fixed_row = column;
            for (int row = fixed_row - 1; row >= 0; --row) {
                double alpha_for_transvection = (-(matrix[row][column] / matrix[fixed_row][column]));
                for (int column_2 = 0; column_2 < matrix_size; ++column_2) {
                    identity_matrix[row][column_2] += alpha_for_transvection * identity_matrix[column][column_2];
                    matrix[row][column_2] += alpha_for_transvection * matrix[column][column_2];
                }
            }
            for (int column_3 = 0; column_3 < matrix_size; ++column_3) {
                identity_matrix[column][column_3] *= 1 / matrix[column][column];
            }
        }
        return 1;
    }
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

    if (inverse_matrix(matrix, matrix_size) == 0) {
        printf("DOES NOT EXIST");
    } else {
        for (int row = 0; row < matrix_size; ++row) {
            for (int column = 0; column < matrix_size; ++column) {
                printf("%0.2lf ", identity_matrix[row][column]);
            }
            printf("\n");
        }
    }

    return 0;
}