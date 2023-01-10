#include <vector>
#include <map>
#include <iostream>
#include <fvhyper/parallel.h>


// Sparse vector and matrix implementation
// Supports indexing, adding elements
// Operations allowed: matrix - vector product



namespace fvhyper {




// Class for a distributed sparse vector
template<class T>
class vector {
private:
    std::map<uint, T> x;
    mpi_wrapper* pool;

public:
    vector();
    vector(uint n);

    void set_pool(mpi_wrapper* pool);

    T get(const uint i);
    T& set(const uint i);

    T& fast_get(const uint i);

    typename std::map<uint, T>::iterator begin();
    typename std::map<uint, T>::iterator end();

    std::map<uint, T>& iterate();

    void print();

};


template<class T>
vector<T>::vector() {}

template<class T>
vector<T>::vector(uint n) {
    x.reserve(n);
}

template<class T>
void vector<T>::set_pool(mpi_wrapper* pool_i) {
    pool = pool_i;
}

template<class T>
T vector<T>::get(const uint i) {
    if (x.find(i) != x.end()) {
        return x[i];
    } else {
        return 0;
    }
}

template<class T>
T& vector<T>::set(const uint i) {
    return x[i];
}

template<class T>
T& vector<T>::fast_get(const uint i) {
    return x.at(i);
}

template<class T>
typename std::map<uint, T>::iterator vector<T>::begin() {
    return x.begin();
}

template<class T>
typename std::map<uint, T>::iterator vector<T>::end() {
    return x.end();
}

template<class T>
std::map<uint, T>& vector<T>::iterate() {
    return x;
}

template<class T>
void vector<T>::print() {
    for (auto& keyval : x) {
        std::cout << "(" << keyval.first << ") -> " << keyval.second << std::endl;
    }
}




// Class for a distributed sparse matrix
template<class T>
class matrix {
private:
    std::vector<T> V;
    std::vector<uint> I;
    std::vector<uint> J;
    std::map<uint, uint> row_indices;
    mpi_wrapper* pool;

public:
    
    matrix();
    matrix(uint n);

    void set_pool(mpi_wrapper* pool);

    void push_back(T v, uint i, uint j);
    void place_new(T v, uint i, uint j);

    void mult(vector<T>& ans, vector<T>& x);

    void print();

};

template<class T>
matrix<T>::matrix() {};

template<class T>
matrix<T>::matrix(uint n) {
    V.reserve(n);
    I.reserve(n);
    J.reserve(n);
};


template<class T>
void matrix<T>::push_back(T v, uint i, uint j) {
    V.push_back(v);
    I.push_back(i);
    J.push_back(j);
}

template<class T>
void matrix<T>::set_pool(mpi_wrapper* pool_i) {
    pool = pool_i;
}



template<class T>
void matrix<T>::place_new(T v, uint i, uint j) {
    uint k = 0;
    if (row_indices.find(i) != row_indices.end()) {
        k = row_indices[i];
        while (k < J.size()) {
            if (J[k] > j) {break;} else {k += 1;}
        }
        // +1 to all higher row indices
        for (auto& keyval : row_indices) {
            if (keyval.first > i) {
                keyval.second += 1;
            }
        }
    } else {
        k = V.size();
        row_indices[i] = k;
    }
    // We are now at position
    V.insert(V.begin() + k, v);
    I.insert(I.begin() + k, i);
    J.insert(J.begin() + k, j);
}



template<class T>
void matrix<T>::print() {
    for (uint k=0; k<V.size(); ++k) {
        std::cout << "(" << I[k] << ", " << J[k] << ") -> " << V[k] << "\n";
    }
}


template<class T>
void matrix<T>::mult(vector<T>& ans, vector<T>& x) {
    // Reset ans
    for (auto& keyval : ans.iterate()) {
        keyval.second = 0;
    }
    // Loop over each non zero lines of matrix
    for (const auto& rowval : row_indices) {
        uint i = rowval.first;
        uint k = rowval.second;
        while (I[k] == i) {
            uint j = J[k];
            ans.set(j) += V[k] * x.get(j);
            k += 1;
        }
    }
}



}


