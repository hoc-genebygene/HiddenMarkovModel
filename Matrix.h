#pragma once

#include <cstdint>
#include <vector>

template <class T>
class Matrix {
    struct MatrixRow; // forward declare this to be used by operator[]

public:
    Matrix(size_t num_rows, size_t num_cols, const T & value) : num_rows_(num_rows), num_cols_(num_cols), vec_(std::vector<T>(num_rows * num_cols, value)) {}

    MatrixRow operator[] (size_t row_index) { return MatrixRow(row_index, *this); }

    size_t GetNumRows() { return num_rows_; }
    size_t GetNumCols() { return num_cols_; }

private:
    struct MatrixRow {
        Matrix & mat;
        size_t row_index;

        MatrixRow(size_t row_index, Matrix & mat) : mat(mat), row_index(row_index) {}

        T& operator[] (size_t col_index) { return mat.vec_[col_index * mat.num_rows_ + row_index]; }
    };

    size_t num_rows_;
    size_t num_cols_;

    std::vector<T> vec_;
};
