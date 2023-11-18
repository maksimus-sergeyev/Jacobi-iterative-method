#pragma once
#include "matrix.h"
#include <iostream>

// method Jacobi
// 
// A * X = B; A, B, X - square matrices
// sufficient condition - diagonal dominance: abs(a[i][i]) > abs(a[i][1]) + ... + abs(a[i][i-1]) + abs(a[i][i+1]) + ... + abs(a[i][n]), forall i = 1,...,n
// 
// A * x = b -> x = C * x + d ; x(i+1) = A * x(i) + B;
//
//                    0         | -A[1][2] / A[1][1] | -A[1][3] / A[1][1] | ... | -A[1][n] / A[1][1] |
//          -A[2][1] / A[2][2]  |          0         | -A[2][3] / A[2][2] | ... | -A[2][n] / A[2][2] |
// A  :=    -A[3][1] / A[3][3]  | -A[3][2] / A[3][3] |          0         | ... | -A[3][n] / A[3][3] |
//          ..........................................................................................
//          -A[n][1] / A[n][n]  | -A[n][2] / A[n][n] | -A[n][3] / A[n][n] | ... |          0         |
//
//
//           B[1][1] / A[1][1]  |  B[1][2] / A[1][1] |  B[1][3] / A[1][1] | ... |  B[1][n] / A[1][1] |
//           B[2][1] / A[2][2]  |  B[2][2] / A[2][2] |  B[2][3] / A[2][2] | ... |  B[2][n] / A[2][2] |
// B  :=     B[3][1] / A[3][3]  |  B[3][2] / A[3][3] |  B[3][3] / A[3][3] | ... |  B[3][n] / A[3][3] |
//          ..........................................................................................
//           B[n][1] / A[n][n]  |  B[n][2] / A[n][n] |  B[n][3] / A[n][n] | ... |  B[n][n] / A[n][n] |
//
// x(0) := B

const double min = -1000;
const double max = 1000;
const int max_iter = 100;

template <class T>
class solver
{
private:
    int size;
    matrix<T> B;
    matrix<T> A;
    matrix<T> X_PREV;
    matrix<T> X_NEXT;
public:

    solver(matrix<T>& _A, matrix<T>& _B) : size(_A.get_row()), B(_B), A(_A), X_PREV(_A.get_row()), X_NEXT(_A.get_row())
    {
        if ((_A.get_row() != _A.get_col())|| (_B.get_row() != _B.get_col()) || (_A.get_row() != _B.get_col())) throw std::invalid_argument("not equal sizes");
    }

    solver(int _size) : size(_size), B(size, size), A(size, size), X_PREV(size, size), X_NEXT(size, size)
    {
        if (size <= 0) throw std::invalid_argument("size should be positive");
    }

    void random_diag_dom_fill()
    {
        B.randomfill();
        A.randomfill();

        double summ = 0;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
                if (i != j) summ += std::abs(static_cast<double>(A[i * size + j]));

            if (summ > std::abs(static_cast<double>(A[i * size + i])))
                A[i * size + i] = summ + static_cast<T>(std::abs((static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * (max - min) + min));
        }
    }

    bool check()// sufficient condition - diagonal dominance : abs(a[i][i]) > abs(a[i][1]) + ... + abs(a[i][i - 1]) + abs(a[i][i + 1]) + ... + abs(a[i][n]), forall i = 1, ..., n
    {
        for (int i = 0; i < size; i++)
        {
            double summ = 0;

            for (int j = 0; j < size; j++)
                if (i != j) summ += std::abs(static_cast<double>(A[i * size + j]));

            if (summ > std::abs(static_cast<double>(A[i * size + i]))) return 0;

        }

        return 1;
    }
    int solution(double eps, bool mod) //mod = 1 => need to check sufficient condition
    {
        if (mod) { if (check()) return 0; }
        
#pragma omp parallel for
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                B[i * size + j] /= A[i * size + i];

#pragma omp parallel for
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                if (i != j) A[i * size + j] = -A[i * size + j] / A[i * size + i];

        for (int i = 0; i < size; i++)
            A[i * size + i] = static_cast<T>(0);

        X_PREV = B;

        parallel_block_mult(A, X_PREV, X_NEXT);

        X_NEXT += B;

        int iter = 0;

        while ((X_NEXT - X_PREV).norm() > static_cast<T>(eps))
        {
            X_PREV = X_NEXT;
            parallel_block_mult(A, X_PREV, X_NEXT);
            X_NEXT += B;
            iter++;
            if (iter > max_iter) return 2;
        }

        //std::cout << iter << "\n";

        return 1;
    }
    matrix<T>& getX() noexcept
    {
        return X_NEXT;
    }
    matrix<T>& getA() noexcept
    {
        return A;
    }
    matrix<T>& getB() noexcept
    {
        return B;
    }
};

