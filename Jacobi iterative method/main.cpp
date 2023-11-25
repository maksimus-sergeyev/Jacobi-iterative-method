#include <iostream>
#include <chrono>
#include <iomanip>
#include <time.h>
#include "solver.h"
#include "matrix.h"

const double EPS = 1E-300;

/*
CPU:
	11th Gen Intel(R) Core(TM) i7-11700F @ 2.50GHz
	Base speed:	2,50 GHz
	Sockets:	1
	Cores:	8
	Logical processors:	16
	Virtualization:	Enabled
	L1 cache:	640 KB
	L2 cache:	4,0 MB
	L3 cache:	16,0 MB

	//micro architecture - rocket lake
ICC
    options: / permissive - / GS / Qopenmp / Qftz / W3 / QxHost / Gy / Zc:wchar_t / Zi / O3 / Fd"x64\Release\vc143.pdb" / Zc : inline / Quse - intel - optimized - headers / D "NDEBUG" / D "_CONSOLE" / D "_UNICODE" / D "UNICODE" / Qipo / Zc : forScope / Oi / MD / FC / Fa"x64\Release\" /EHsc /nologo /Qparallel /Fo"x64\Release\" /Ot /Fp"x64\Release\Jacobi iterative method.pch" 


// EPS = 1E-300; size = 8000; max = MAXM = 1000; min = MINM = -1000;
// time of solving ~ 66000 ms, precision ~1e-17, iter ~7

*/

int main()
{
    srand(time(NULL));

    int size = 8000;

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