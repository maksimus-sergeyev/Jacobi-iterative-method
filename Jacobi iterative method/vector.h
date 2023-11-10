#pragma once
#include <iostream>

const double MAXV = 1E10;
const double MINV = -1E10;

template <class T>
class vector
{
private:
    int size;
    T* arr;
public:
    vector(int _size = 3) : size(_size)
    {
        if (size <= 0) throw std::exception("size should be positive");

        arr = new T[size]();
    }
    vector(vector& v)
    {
        size = v.size;
        arr = new T[size]();
        for (int i = 0; i < size; i++)
            arr[i] = v[i];
    }
    vector(vector&& v)
    {
        arr = nullptr;
        std::swap(size, v.size);
        std::swap(arr, v.arr);
    }
    ~vector()
    {
        delete[] arr;
        arr = nullptr;
        size = 0;
    }
    T& operator [](int i)
    {
        return arr[i];
    }
    const T& operator [](int i) const
    {
        return arr[i];
    }
    void randomfill()
    {
        for (int i = 0; i < size; i++)
            arr[i] = static_cast<T>((static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * (MAXV - MINV) + MINV);
    }
    int getsize()
    {
        return size;
    }
    bool operator==(const vector& v) 
    {
        if (size != v.size) return false;

        for (int i = 0; i < size; i++)
            if (arr[i] != v[i]) return false;

        return true;
    }
    bool operator!=(const vector& v) 
    {
        return !(*this == v);
    }
    vector& operator=(const vector& v)
    {
        if (this == &v) return *this;

        if (size != v.size)
        {
            T* tmp = new T[v.size]();
            delete[] arr;
            arr = tmp;
            size = v.size;
        }

#pragma omp parallel for
        for (int i = 0; i < size; i++)
            arr[i] = v[i];

        return *this;
    }
    vector& operator=(vector&& v) 
    {
        std::swap(arr, v.arr);
        std::swap(size, v.size);
    }
    vector operator-(vector& v)
    {
        if (size != v.size) throw std::exception("not equal sizes");

        vector<T> res(size);

#pragma omp parallel for
        for (int i = 0; i < size; i++)
            res[i] = arr[i] - v[i];

        return res;
    }

    vector operator+(vector& v)
    {
        if (size != v.size) throw std::exception("not equal sizes");

        vector<T> res(size);

#pragma omp parallel for
        for (int i = 0; i < size; i++)
            res[i] = arr[i] + v[i];

        return res;
    }
    vector& operator-=(vector& v)
    {
        if (size != v.size) throw std::exception("not equal sizes");

#pragma omp parallel for
        for (int i = 0; i < size; i++)
            arr[i] -= v[i];

        return *this;
    }
    T abs()
    {
        T abs = 0;
        for (int i = 0; i < size; i++)
            abs += arr[i] * arr[i];

        return abs;
    }
};

template <class T>
std::istream& operator>> (std::istream& in, vector<T>& v)
{
    for (int i = 0; i < v.getsize(); i++)
        in >> v[i];

    return in;
}

template <class T>
std::ostream& operator<< (std::ostream& out, vector<T>& v)
{
    for (int i = 0; i < v.getsize(); i++)
        out << v[i] << "\n";

    return out;
}
