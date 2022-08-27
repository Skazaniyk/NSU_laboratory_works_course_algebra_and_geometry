#include <stdio.h>

#define SIZE_LIMIT 200

void system_of_linear_equations(double matrix[SIZE_LIMIT][SIZE_LIMIT], int matrix_size, double free_terms[SIZE_LIMIT]) {
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
            free_terms[current_column] += free_terms[non_zero_row];
        }

        for (int row = current_row + 1; row < matrix_size; ++row) {
            double first_element_of_current_row = matrix[row][current_column];
            double alpha_for_transvection = (-(first_element_of_current_row / non_zero_element));
            for (int column = 0; column < matrix_size; ++column) {
                matrix[row][column] += alpha_for_transvection * matrix[current_row][column];
            }
            free_terms[row] += alpha_for_transvection * free_terms[current_row];

        }
        ++current_row;
    }

    for (int column = matrix_size - 1; column >= 0; --column) {
        int fixed_row = column;
        for (int row = fixed_row - 1; row >= 0; --row) {
            double alpha_for_transvection = (-(matrix[row][column] / matrix[fixed_row][column]));
            for (int column_2 = 0; column_2 < matrix_size; ++column_2) {
                matrix[row][column_2] += alpha_for_transvection * matrix[column][column_2];
            }
            free_terms[row] += alpha_for_transvection * free_terms[column];
        }
        free_terms[column] *= 1 / matrix[column][column];
    }

}

int main() {
    int matrix_size;
    double matrix[SIZE_LIMIT][SIZE_LIMIT] = {0}, free_terms[SIZE_LIMIT] = {0};

    scanf("%d", &matrix_size);
    for (int row = 0; row < matrix_size; ++row) {
        for (int column = 0; column < matrix_size; ++column) {
            scanf("%lf", &matrix[row][column]);
        }
    }
    for (int row = 0; row < matrix_size; ++row) {
        scanf("%lf", &free_terms[row]);
    }

    system_of_linear_equations(matrix, matrix_size, free_terms);

    for (int row = 0; row < matrix_size; ++row) {
        printf("%0.2lf ", free_terms[row]);
    }
    return 0;
}