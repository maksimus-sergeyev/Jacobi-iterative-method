#include <iostream>
#include <chrono>
#include <iomanip>
#include <time.h>
#include "solver.h"
#include "vector.h"
#include "square_matrix.h"

const double EPS = 1E-300;

// EPS = 1E-300; size = 30000; headers: max = MAXM = MAXV = 1E10; min = MINM = MINV = -1E10;
// time of solving ~ 1500 ms, precision ~1e-8, iter ~5
// DPCPP COMPILER {CPU 11th Gen Intel(R) Core(TM) i7 - 11700F @ 2.50GHz, 16gb RAM, windows 10 64}

int main()
{
    srand(time(NULL));     
    
    int size = 30000;
                  
    using T = double;

    vector<T> x(size);

    square_matrix<T> A(size);
    vector<T>b(size);

    solver<T>s(size);

    s.random_diag_dom_fill();

    A = s.get_m();

    b = s.get_v();

    //std::cout << std::setprecision(10) << A;

    //std::cout << std::setprecision(10) << b;

    auto begin = std::chrono::steady_clock::now();

    int flag = s.solution(EPS, 0);

    auto end = std::chrono::steady_clock::now();

    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    std::cout << "The time: " << elapsed_ms.count() << " ms\n";

    x = s.get_x();

    std::cout << (A * x - b).abs() << "\n";

    return 0;
}