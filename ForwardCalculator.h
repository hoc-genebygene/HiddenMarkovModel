#pragma once

#include "Matrix.h"
#include "PairHMM.h"

namespace PairHMM {
class ForwardCalculator {
public:
    ForwardCalculator(PairHMM hmm);

    void Calculate();

    void Print();
private:
    void CalculateForwardMatch(size_t r, size_t c);
    void CalculateForwardX(size_t r, size_t c);
    void CalculateForwardY(size_t r, size_t c);

    void SetupInitialConditions();
private:
    PairHMM hmm_;

    size_t num_rows_;
    size_t num_cols_;

    Matrix<double> f_match_;
    Matrix<double> f_x_;
    Matrix<double> f_y_;
};
} // namespace PairHMM
