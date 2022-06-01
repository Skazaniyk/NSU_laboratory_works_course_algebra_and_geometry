#pragma once

#ifndef HAMMING_CODE
#define HAMMING_CODE

#include <cstring>
#include <iomanip>

template<typename T>
class Mat {
private:
    int cols;
    int rows;
    T *M;

public:
    Mat(int rows, int cols);

    Mat();

    Mat(int rows, int cols, const T &value);

    Mat(const Mat<T> &m);

    ~Mat();

    T &operator()(int row, int col);

    const T &operator()(int row, int col) const;

    Mat<T> &operator=(const Mat<T> &m);

    Mat<T> &operator=(const T &s);

    Mat operator+(const Mat<T> &b);

    Mat operator+(const T &s);

    Mat operator-(const Mat<T> &b);

    Mat operator-(const T &s);

    Mat operator*(const T &s);

    Mat operator*(const Mat<T> &a);

    bool operator==(const Mat<T> &b);

    Mat<T> &operator=(Mat<T> &m);

    Mat<T> &operator=(T &s);

    Mat operator+(Mat<T> &b);

    Mat operator+(T &s);

    Mat operator-(Mat<T> &b);

    Mat operator-(T &s);

    Mat operator*(T &s);

    Mat operator*(Mat<T> &a);

    bool operator==(Mat<T> &b);

    Mat<T> transposition();

    double determination();

    static Mat<T> eye(int rows_and_cols);

    static Mat<T> zeros(int rows, int cols);

    int get_num_rows();

    int get_num_cols();

    bool is_square();

    void print_matrix();

    void matrix_by_mod(int mod);
};

#endif //HAMMING_CODE

template<class T>
Mat<T>::Mat() {
    rows = 1;
    cols = 1;
    M = new T[1];
    M[0] = 0.0;
}


template<typename T>
Mat<T>::Mat(int rows, int cols): Mat(rows, cols, 0) {}

template<typename T>
Mat<T>::Mat(int rows, int cols, const T &value) {
    M = (T *) new T[rows * cols];
    for (int i = 0; i < rows * cols; i++) {
        M[i] = value;
    }
    this->rows = rows;
    this->cols = cols;
}

template<typename T>
Mat<T>::Mat(const Mat &m):cols(m.cols), rows(m.rows) {
    M = (T *) new T[rows * cols];
    for (int i = 0; i < cols * rows; ++i) {
        M[i] = m.M[i];
    }
}


template<typename T>
Mat<T>::~Mat() {
    delete[] M;
}

template<typename T>
T &Mat<T>::operator()(int row, int col) {
    return M[row * cols + col];
}

template<typename T>
const T &Mat<T>::operator()(int row, int col) const {
    return M[row * cols + col];
}

template<typename T>
Mat<T> Mat<T>::zeros(int rows, int cols) {
    return Mat<int>(rows, cols, 0);
}

template<typename T>
Mat<T> Mat<T>::eye(int rows_and_cols) {
    Mat<int> Eye(rows_and_cols, rows_and_cols, 0);
    for (int i = 0; i < rows_and_cols; ++i) {
        Eye(i, i) = 1;
    }
    return Eye;
}

template<typename T>
bool Mat<T>::operator==(const Mat<T> &b) {
    if (cols != b.cols || rows != b.rows) {
        return false;
    }
    for (int i = 0; i < b.cols * b.rows; i++) {
        if (M[i] != b.M[i]) {
            return false;
        }
    }
    return true;
}

template<typename T>
bool Mat<T>::operator==(Mat<T> &b) {
    if (cols != b.cols || rows != b.rows) {
        return false;
    }
    for (int i = 0; i < b.cols * b.rows; i++) {
        if (M[i] != b.M[i]) {
            return false;
        }
    }
    return true;
}

template<typename T>
Mat<T> Mat<T>::operator*(const Mat<T> &b) {
    Mat<T> this_m(*this);
    Mat<T> result_m(rows, b.cols, 0);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < b.cols; j++) {
            for (int p = 0; p < cols; p++) {
                result_m(i, j) += this_m(i, p) * b(p, j);
            }
        }
    }
    return result_m;
}

template<typename T>
Mat<T> Mat<T>::operator*(Mat<T> &b) {
    Mat<T> this_m(*this);
    Mat<T> result_m(rows, b.cols, 0);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < b.cols; j++) {
            for (int p = 0; p < cols; p++) {
                result_m(i, j) += this_m(i, p) * b(p, j);
            }
        }
    }
    return result_m;
}

template<typename T>
Mat<T> Mat<T>::operator*(const T &s) {
    Mat<T> this_m(*this);
    Mat<T> result_m(this_m);
    for (int i = 0; i < rows * cols; i++) {
        result_m[i] *= s;
    }
    return result_m;
}

template<typename T>
Mat<T> Mat<T>::operator*(T &s) {
    Mat<T> this_m(*this);
    Mat<T> result_m(this_m);
    for (int i = 0; i < rows * cols; i++) {
        result_m[i] *= s;
    }
    return result_m;
}

template<typename T>
Mat<T> &Mat<T>::operator=(const Mat<T> &m) {
    rows = m.rows;
    cols = m.cols;
    M = new T[rows * cols];
    for (int i = 0; i < rows * cols; ++i) {
        M[i] = m.M[i];
    }
    return *this;
}

template<typename T>
Mat<T> &Mat<T>::operator=(Mat<T> &m) {
    rows = m.rows;
    cols = m.cols;
    M = new T[rows * cols];
    for (int i = 0; i < rows * cols; ++i) {
        M[i] = m.M[i];
    }
    return *this;
}

template<typename T>
Mat<T> &Mat<T>::operator=(const T &s) {
    for (int i = 0; i < rows * cols; ++i) {
        M[i] = s;
    }
    return *this;
}

template<typename T>
Mat<T> &Mat<T>::operator=(T &s) {
    for (int i = 0; i < rows * cols; ++i) {
        M[i] = s;
    }
    return *this;
}

template<typename T>
Mat<T> Mat<T>::operator+(const Mat<T> &b) {
    Mat<T> result_m(b);
    if (rows == b.rows && cols == b.cols) {
        for (int i = 0; i < cols * rows; ++i) {
            result_m.M[i] += M[i];
        }
    }
    return result_m;
}

template<typename T>
Mat<T> Mat<T>::operator+(Mat<T> &b) {
    Mat<T> result_m(b);
    if (rows == b.rows && cols == b.cols) {
        for (int i = 0; i < cols * rows; ++i) {
            result_m.M[i] += M[i];
        }
    }
    return result_m;
}

template<typename T>
Mat<T> Mat<T>::operator-(const Mat<T> &b) {
    Mat<T> result_m(b);
    if (rows == b.rows && cols == b.cols) {
        for (int i = 0; i < cols * rows; ++i) {
            result_m.M[i] += M[i];
        }
    }
    return result_m;
}

template<typename T>
Mat<T> Mat<T>::operator-(Mat<T> &b) {
    Mat<T> result_m(b);
    if (rows == b.rows && cols == b.cols) {
        for (int i = 0; i < cols * rows; ++i) {
            result_m.M[i] += M[i];
        }
    }
    return result_m;
}

template<typename T>
Mat<T> Mat<T>::operator+(const T &s) {
    Mat<T> result_m(rows, cols, s);
    for (int i = 0; i < cols * rows; ++i) {
        result_m.M[i] += M[i];
    }
    return result_m;
}

template<typename T>
Mat<T> Mat<T>::operator+(T &s) {
    Mat<T> result_m(rows, cols, s);
    for (int i = 0; i < cols * rows; ++i) {
        result_m.M[i] += M[i];
    }
    return result_m;
}

template<typename T>
Mat<T> Mat<T>::operator-(const T &s) {
    Mat<T> result_m(rows, cols, s);
    for (int i = 0; i < cols * rows; ++i) {
        result_m.M[i] = M[i] - result_m.M[i];
    }

    return result_m;
}

template<typename T>
Mat<T> Mat<T>::operator-(T &s) {
    Mat<T> result_m(rows, cols, s);
    for (int i = 0; i < cols * rows; ++i) {
        result_m.M[i] = M[i] - result_m.M[i];
    }

    return result_m;
}

template<typename T>
Mat<T> Mat<T>::transposition() {
    int _rows = this->rows;
    int _cols = this->cols;
    Mat<T> result_m = Mat(_cols, _rows);
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _cols; ++j) {
            result_m.M[j * _rows + i] = M[j * _rows + i];
        }
    }
    return result_m;
}


template<typename T>
int Mat<T>::get_num_cols() {
    return this->cols;
}

template<typename T>
int Mat<T>::get_num_rows() {
    return this->rows;
}

template<typename T>
bool Mat<T>::is_square() {
    if (cols == rows) {
        return true;
    } else {
        return false;
    }
}

template<typename T>
double Mat<T>::determination() {
    if (!is_square()) {
        throw std::invalid_argument("Cannot compute the determinant of matrix that is not square.");
    } else {
        Mat<T> matrix(M);
        int number_of_cols = get_num_cols(M);
        int number_of_rows = get_num_rows(M);

        int non_zero_row = 0, current_row = 0;
        for (int current_column = 0; current_column < number_of_cols; ++current_column) {

            int column_consists_of_zeros = 1;
            for (int row = current_row; row < number_of_rows; ++row) {
                if (matrix[row * number_of_rows + current_column] != 0) {
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
                for (int column = 0; column < number_of_cols; ++column) {
                    matrix[current_column][column] += matrix[non_zero_row][column];
                }
            }

            for (int row = current_row + 1; row < number_of_rows; ++row) {
                double first_element_of_current_row = matrix[row][current_column];
                for (int column = 0; column < number_of_cols; ++column) {
                    matrix[row][column] +=
                            (-(first_element_of_current_row / non_zero_element)) * matrix[current_row][column];
                }
            }
            ++current_row;
        }

        double determinant = 1;
        for (int diagonal_index = 0; diagonal_index < number_of_rows; ++diagonal_index) {
            determinant *= matrix[diagonal_index][diagonal_index];
        }
        return determinant;
    }
}


template<class T>
void Mat<T>::print_matrix() {
    int nRows = this->GetNumRows();
    int nCols = this->GetNumCols();

    for (int row = 0; row < nRows; ++row) {
        for (int col = 0; col < nCols; ++col) {
            std::cout << std::fixed << std::setprecision(2) << this->GetElement(row, col) << " ";
        }
        std::cout << std::endl;
    }
}

template<class T>
void Mat<T>::matrix_by_mod(int mod) {
    for (int i = 0; i < cols * rows; i++) {
        M[i] = (mod + (M[i] % mod)) % mod;
    }
}