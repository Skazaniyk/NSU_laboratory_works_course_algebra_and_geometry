#include <iostream>
#include "../Mat.h"
#include "../Hamming.h"

using namespace std;

int main() {
    int x = 2;
    int y = 2;
    Mat<int> matrix_2(x, y, 5);
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            cout << matrix_2(i, j) << " ";
        }
        cout << endl;
    }

    double a = matrix_2.determination();

    return 0;
}
