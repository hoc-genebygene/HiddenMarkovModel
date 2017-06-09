#include "BackwardCalculator.h"

#include <iostream>

namespace PairHMM {

BackwardCalculator::BackwardCalculator(PairHMM hmm) :
    hmm_(std::move(hmm)),
    num_rows_(hmm_.GetSeq2().size() + 1),
    num_cols_(hmm_.GetSeq1().size() + 1),
    b_match_(Matrix<double>(num_rows_, num_cols_, 0.0)),
    b_x_(Matrix<double>(num_rows_, num_cols_, 0.0)),
    b_y_(Matrix<double>(num_rows_, num_cols_, 0.0))
    {
        SetupInitialConditions();
    }

void BackwardCalculator::Calculate() {
    for (size_t j = 1; j < num_rows_; ++j) {
        for (size_t i = 1; i < num_cols_; ++i) {
            size_t r = num_rows_ - 1 - j;
            size_t c = num_cols_ - 1 - i;

            if (r == num_rows_ - 2 && c == num_cols_ - 2) {
                continue;
            }

            CalculateBackwardMatch(r, c);
            CalculateBackwardX(r, c);
            CalculateBackwardY(r, c);
        }
    }
}

void BackwardCalculator::Print() {
    std::cout << "Match: " << std::endl;
    for (size_t j = 1; j < num_rows_; ++j) {
        for (size_t i = 1; i < num_cols_; ++i) {
            size_t r = num_rows_ - 1 - j;
            size_t c = num_cols_ - 1 - i;

            std::cout << b_match_[r][c] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "InsertX: " << std::endl;
    for (size_t j = 1; j < num_rows_; ++j) {
        for (size_t i = 1; i < num_cols_; ++i) {
            size_t r = num_rows_ - 1 - j;
            size_t c = num_cols_ - 1 - i;

            std::cout << b_x_[r][c] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "InsertY: " << std::endl;
    for (size_t j = 1; j < num_rows_; ++j) {
        for (size_t i = 1; i < num_cols_; ++i) {
            size_t r = num_rows_ - 1 - j;
            size_t c = num_cols_ - 1 - i;

            std::cout << b_y_[r][c] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void BackwardCalculator::CalculateBackwardMatch(size_t r, size_t c) {
    double term1 = (1.0 - 2.0 * delta() - tau()) * hmm_.GetEmissionMatrix().GetEmissionProbability(ConvertCharToNucleotide(hmm_.GetSeq1()[c]), ConvertCharToNucleotide(hmm_.GetSeq2()[r])) * b_match_[r+1][c+1];
    double term2 = delta() * hmm_.GetGapEmissionProbability() * (b_x_[r][c+1] + b_y_[r+1][c]);

    double score = term1 + term2;

    b_match_[r][c] = score;
}

void BackwardCalculator::CalculateBackwardX(size_t r, size_t c) {
    double score = (1.0 - epsilon() - tau()) * hmm_.GetEmissionMatrix().GetEmissionProbability(ConvertCharToNucleotide(hmm_.GetSeq1()[c]), ConvertCharToNucleotide(hmm_.GetSeq2()[r])) * b_match_[r+1][c+1] + epsilon() * hmm_.GetGapEmissionProbability() * b_x_[r][c+1];
    b_x_[r][c] = score;
}

void BackwardCalculator::CalculateBackwardY(size_t r, size_t c) {
    double score = (1.0 - epsilon() - tau()) * hmm_.GetEmissionMatrix().GetEmissionProbability(ConvertCharToNucleotide(hmm_.GetSeq1()[c]), ConvertCharToNucleotide(hmm_.GetSeq2()[r])) * b_match_[r+1][c+1] + epsilon() * hmm_.GetGapEmissionProbability() * b_y_[r+1][c];
    b_y_[r][c] = score;
}

void BackwardCalculator::SetupInitialConditions() {
    for (size_t j = 0; j < num_cols_; --j) {
        for (size_t i = 0; i < num_rows_; --i) {
            size_t c = num_cols_ - 1 - j;
            size_t r = num_rows_ - 1 - i;

            b_match_[r][c] = 0.0;
            b_y_[r][c] = 0.0;
            b_x_[r][c] = 0.0;
        }
    }

    b_match_[num_rows_-2][num_cols_-2] = tau();
    b_x_[num_rows_-2][num_cols_-2] = tau();
    b_y_[num_rows_-2][num_cols_-2] = tau();
}

} // namespace PairHMM
