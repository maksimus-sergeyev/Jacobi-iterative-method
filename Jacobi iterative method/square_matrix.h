#pragma once

#include <iostream>

const double MAXM = 1E10;
const double MINM = -1E10;

template <class T>
class square_matrix
{
private:
    int size;
    T* arr;
public:
    square_matrix(int _size = 1000) : size(_size)
    {
        if (size <= 0) throw std::exception("size should be positive");

        arr = new T[size * size]();
    }
    square_matrix(square_matrix& m)
    {
        size = m.size;
        arr = new T[size * size]();
#pragma omp parallel for
        for (int i = 0; i < size * size; i++)
            arr[i] = m[i];
    }
    square_matrix(square_matrix&& m)
    {
        arr = nullptr;
        std::swap(size, m.size);
        std::swap(arr, m.arr);
    }
    ~square_matrix()
    {
        delete[] arr;
        arr = nullptr;
        size = 0;
    }
    T& operator[](int i)
    {
        return arr[i];
    }
    const T& operator[](int i) const
    {
        return arr[i];
    }
    void randomfill()noexcept
    {
        for (int i = 0; i < size * size; i++)
            arr[i] = static_cast<T>((static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * (MAXM - MINM) + MINM);
    }
    int getsize() const noexcept
    {
        return size;
    }
    T* getarr() const noexcept
    {
        return arr;
    }
    square_matrix& operator=(const square_matrix& m)
    {
        if (this == &m) return *this;

        if (size != m.size)
        {
            T* tmp = new T[m.size]();
            delete[] arr;
            arr = tmp;
            size = m.size;
        }

#pragma omp parallel for
        for (int i = 0; i < size * size; i++)
            arr[i] = m[i];

        return *this;
    }
    square_matrix& operator=(square_matrix&& m) noexcept
    {
        std::swap(arr, m.arr);
        std::swap(size, m.size);

        return *this;
    }
    bool operator == (const square_matrix& m) const noexcept
    {
        if (size != m.size) return false;

        for (int i = 0; i < size * size; i++)
            if (arr[i] != m[i]) return false;

        return true;
    }
    bool operator!= (const square_matrix& m) const
    {
        return !(*this == m);
    }
    square_matrix operator*(const square_matrix& m)
    {
        if (size != m.size) throw std::out_of_range("matrices sizes should be equal!");

        square_matrix res(size);

#pragma omp parallel for
        for (int i = 0; i < size; i++)
            for (int k = 0; k < size; k++)
                for (int j = 0; j < size; j++)
                    res[i * size + j] += arr[i * size + k] * m[k * size + j];

        return res;
    }
    square_matrix operator-(const square_matrix& m)
    {
        if (size != m.size) throw std::out_of_range("matrices sizes should be equal!");

        square_matrix res(size);
#pragma omp parallel for
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                res[i * size + j] = arr[i * size + j] - m[i * size + j];

        return res;
    }
    T abs() {
        T res = 0;
        for (int i = 0; i < size * size; i++)
            res += arr[i] * arr[i];
        return res;
    }
};

template <class T>
std::istream& operator>> (std::istream& in, square_matrix<T>& m)
{
    for (int i = 0; i < m.getsize(); i++)
        for (int j = 0; j < m.getsize(); j++)
            in >> m[i * m.getsize() + j];
    return in;
}

template <class T>
std::ostream& operator<< (std::ostream& out, square_matrix<T>& m)
{
    for (int i = 0; i < m.getsize(); i++)
    {
        for (int j = 0; j < m.getsize(); j++)
            out << m[i * m.getsize() + j] << " ";

        out << "\n";
    }
    return out;
}



template <class T>
vector<T> operator*(square_matrix<T>& m, vector<T>& v)
{
    int size = v.getsize();
    if (size != m.getsize()) throw std::exception("not equal sizes");

    vector<T> vec(size);
#pragma omp parallel for
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            vec[i] += m[i * size + j] * v[j];

    return vec;
}
