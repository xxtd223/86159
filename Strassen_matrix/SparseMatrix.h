#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <numeric> 
#include <chrono>
#include <functional>
#include <cmath>
#include <queue>
#include <unordered_map>
#include "Matrix.h" 
using namespace std;

template <class T>
class SparseMatrix {
public:
    size_t rows, cols;
    unordered_map<size_t, unordered_map<size_t, T>> data;

    SparseMatrix(size_t r, size_t c) : rows(r), cols(c) {}

    Matrix<T> toMatrix() const;
    SparseMatrix<T> operator*(const SparseMatrix<T>& other) const;

    void set(size_t r, size_t c, T value) {
        if (value != 0) {
            data[r][c] = value;
        }
    }

    T get(size_t r, size_t c) const {
        auto row_it = data.find(r);
        if (row_it != data.end()) {
            auto col_it = row_it->second.find(c);
            if (col_it != row_it->second.end()) return col_it->second;
        }
        return 0;
    }

    void print() const {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                cout << get(i, j) << " ";
            }
            cout << endl;
        }
    }
};

//将稀疏矩阵转换为普通矩阵
template <class T>
Matrix<T> SparseMatrix<T>::toMatrix() const {
    Matrix<T> dense(rows, cols);
    for (const auto& row_pair : data) {
        for (const auto& col_pair : row_pair.second) {
            dense.data[row_pair.first][col_pair.first] = col_pair.second;
        }
    }
    return dense;
}

//优化乘法
template <class T>
SparseMatrix<T> SparseMatrix<T>::operator*(const SparseMatrix<T>& other) const {
    if (cols != other.rows) {
        throw invalid_argument("乘法行列不匹配");
    }

    SparseMatrix<T> result(rows, other.cols);

    for (const auto& row : data) {
        for (const auto& elem : row.second) {
            size_t k = elem.first;  
            if (other.data.count(k)) {  
                for (const auto& other_elem : other.data.at(k)) {
                    size_t j = other_elem.first;  
                    T product = elem.second * other_elem.second;
                    if (product != 0) {
                        result.data[row.first][j] += product;
                    }
                }
            }
        }
    }

    return result;
}
