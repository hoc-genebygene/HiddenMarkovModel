#pragma once

#include "Matrix.h"
#include "PairHMM.h"
#include "STATE.h"

namespace PairHMM {
class BackwardCalculator {
public:
    BackwardCalculator(PairHMM hmm);

    void Calculate();
    void Print();

private:
    void CalculateBackwardMatch(size_t r, size_t c);
    void CalculateBackwardX(size_t r, size_t c);
    void CalculateBackwardY(size_t r, size_t c);

    void SetupInitialConditions();

    double delta() {
        return hmm_.GetTransitionMatrix().GetTransitionProbability(MATCH, INSERT_X);
    }

    double epsilon() {
        return hmm_.GetTransitionMatrix().GetTransitionProbability(INSERT_X, INSERT_X);
    }

    double tau() {
        return hmm_.GetTransitionMatrix().GetTransitionProbability(MATCH, END);
    }
    
private:
    PairHMM hmm_;
    
    size_t num_rows_;
    size_t num_cols_;
    
    Matrix<double> b_match_;
    Matrix<double> b_x_;
    Matrix<double> b_y_;
};
} // namespace PairHMM
