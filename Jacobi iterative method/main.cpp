#include <iostream>
#include <chrono>
#include <iomanip>
#include <time.h>
#include "solver.h"
#include "vector.h"
#include "matrix.h"

const double EPS = 1E-300;

// EPS = 1E-300; size = 4000; max = MAXM = 1000; min = MINM = -1000;
// time of solving ~ 20000 ms, precision ~1e-17, iter ~7
// ICC COMPILER {CPU 11th Gen Intel(R) Core(TM) i7 - 11700F @ 2.50GHz, 16gb RAM, windows 10 64}

int main()
{
    srand(time(NULL));     
    
    int size = 4000;
                  
    using T = double;

    matrix<T> X(size, size);

    matrix<T> A(size, size);
    matrix<T> B(size, size);

    solver<T> s(size);

    s.random_diag_dom_fill();

    A = s.getA();

    B = s.getB();

    auto begin = std::chrono::steady_clock::now();

    int flag = s.solution(EPS, 0);

    auto end = std::chrono::steady_clock::now();

    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    if (flag == 1) 
    {
        matrix<T> tmp(size, size);

        X = s.getX();

        parallel_block_mult(A, X, tmp);

        std::cout << "The time: " << elapsed_ms.count() << " ms\n";

        std::cout << (tmp - B).norm() << "\n";
    }

    return 0;
}