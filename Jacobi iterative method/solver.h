#pragma once
#include "vector.h"
#include "square_matrix.h"
#include <iostream>

// method Jacobi
// A * x = b
// sufficient condition: abs(a[i][i]) > abs(a[i][1]) + ... + abs(a[i][i-1]) + abs(a[i][i+1]) + ... + abs(a[i][n]), forall i = 1,...,n
// A * x = b -> x = C * x + d ; x(i+1) = C * x(i) + d;
//
//                    0         | -A[1][2] / A[1][1] | -A[1][3] / A[1][1] | ... | -A[1][n] / A[1][1] |
//          -A[2][1] / A[2][2]  |          0         | -A[2][3] / A[2][2] | ... | -A[2][n] / A[2][2] |
// C   =    -A[3][1] / A[3][3]  | -A[3][2] / A[3][3] |          0         | ... | -A[3][n] / A[3][3] |
//          ..........................................................................................
//          -A[n][1] / A[n][n]  | -A[n][2] / A[n][n] | -A[n][3] / A[n][n] | ... |          0         |
//
//
// d = (b[1] / A[1][1] | b[2] / A[2][2] | .... | b[n]/A[n][n]) // <- transp
//
// x(0) = (b[1] / A[1][1] | b[2] / A[2][2] | .... | b[n]/A[n][n]) // <- transp

const double min = -1E10;
const double max = 1E10;
const int max_iter = 100;

template <class T>
class solver
{
private:
    int size;
    vector<T> v; // in solution v := d
    square_matrix<T> m; //in solution m:= C
    vector<T> x_prev;
    vector<T> x_next;
public:

    solver(square_matrix<T>& _m, vector<T>& _v) : size(_m.get_size()), v(_v), m(_m), x_prev(_m.get_size()), x_next(_m.get_size()),
    {
        if (_m.get_size() != _v.get_size()) throw std::exception("not equal sizes");
    }

    solver(int _size) : v(_size), m(_size), x_prev(_size), x_next(_size), size(_size)
    {
        if (_size <= 0) throw std::exception("size should be positive");
    }

    void random_diag_dom_fill()
    {
        v.randomfill();
        m.randomfill();

        double summ = 0;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
                if (i != j) summ += std::abs(static_cast<double>(m[i * size + j]));

            if (summ > std::abs(static_cast<double>(m[i * size + i])))
                m[i * size + i] = summ + static_cast<T>(std::abs((static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * (max - min) + min));
        }
    }

    bool check() // 
    {
        for (int i = 0; i < size; i++)
        {
            double summ = 0;

            for (int j = 0; j < size; j++)
                if (i != j) summ += std::abs(static_cast<double>(m[i * size + j]));

            if (summ > std::abs(static_cast<double>(m[i * size + i]))) return 0;

        }

        return 1;
    }
    int solution(double eps, bool mod) //mod = 1 => need to check sufficient condition
    {
        if (mod) { if (check()) return 0; }

        for (int i = 0; i < size; i++)
            v[i] /= m[i * size + i];

#pragma omp parallel for
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                if (i != j) m[i * size + j] = -m[i * size + j] / m[i * size + i];

        for (int i = 0; i < size; i++)
            m[i * size + i] = static_cast<T>(0);

        x_prev = v;

        x_next = m * x_prev + v;

        //int iter = 0;

        while ((x_next - x_prev).abs() > static_cast<T>(eps))
        {
            x_prev = x_next;
            x_next = m * x_prev + v;
            //iter++;
            //if (iter > max_iter) return 2;
        }

        //std::cout << iter << "\n";

        return 1;
    }
    vector<T>& get_x()
    {
        return x_next;
    }

    square_matrix<T>& get_m()
    {
        return m;
    }

    vector<T>& get_v()
    {
        return v;
    }
};

