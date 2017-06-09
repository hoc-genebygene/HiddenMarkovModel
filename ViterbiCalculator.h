#pragma once

#include "Matrix.h"
#include "PairHMM.h"

namespace PairHMM {
class ViterbiCalculator {
public:
    ViterbiCalculator(PairHMM hmm);

    void Calculate();

    void Print();

private:
    void SetupInitialConditions();

    double delta();
    double epsilon();
    double tau();

    void CalculateMatchProbability(size_t r, size_t c);
    void CalculateInsertXProbability(size_t r, size_t c);
    void CalculateInsertYProbability(size_t r, size_t c);

private:
    PairHMM hmm_;

    size_t num_rows_;
    size_t num_cols_;

    Matrix<double> v_match_;
    Matrix<double> v_x_;
    Matrix<double> v_y_;
};
}
