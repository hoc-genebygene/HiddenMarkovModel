#include "ViterbiCalculator.h"

#include <cmath>
#include <iostream>

namespace PairHMM {
ViterbiCalculator::ViterbiCalculator(PairHMM hmm) :
    hmm_(std::move(hmm)),
    num_rows_(hmm_.GetSeq2().size() + 1),
    num_cols_(hmm_.GetSeq1().size() + 1),
    v_match_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
    v_x_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
    v_y_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0)))
{
    SetupInitialConditions();
}

void ViterbiCalculator::Calculate() {
    // Traverse column by column
    for (size_t c = 1; c < num_cols_; ++c) {
        for (size_t r = 1; r < num_rows_; ++r) {
            CalculateMatchProbability(r, c);
            CalculateInsertXProbability(r, c);
            CalculateInsertYProbability(r, c);
        }
    }
}

void ViterbiCalculator::Print() {
    std::cout << "Match: " << std::endl;
    for (size_t r = 0; r < num_rows_; ++r) {
        for (size_t c = 0; c < num_cols_; ++c) {
            std::cout << std::pow(10, v_match_[r][c]) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "InsertX: " << std::endl;
    for (size_t r = 0; r < num_rows_; ++r) {
        for (size_t c = 0; c < num_cols_; ++c) {
            std::cout << std::pow(10, v_x_[r][c]) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "InsertY: " << std::endl;
    for (size_t r = 0; r < num_rows_; ++r) {
        for (size_t c = 0; c < num_cols_; ++c) {
            std::cout << std::pow(10, v_y_[r][c]) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void ViterbiCalculator::SetupInitialConditions() {
    v_match_[0][0] = std::log10(1.0);
}

double ViterbiCalculator::delta() {
    return hmm_.GetTransitionMatrix().GetTransitionProbability(MATCH, INSERT_X);
}

double ViterbiCalculator::epsilon() {
    return hmm_.GetTransitionMatrix().GetTransitionProbability(INSERT_X, INSERT_X);
}

double ViterbiCalculator::tau() {
    return hmm_.GetTransitionMatrix().GetTransitionProbability(MATCH, END);
}

void ViterbiCalculator::CalculateMatchProbability(size_t r, size_t c) {
    double best_score = std::log10(hmm_.GetTransitionMatrix().GetTransitionProbability(MATCH, MATCH)) + v_match_[r-1][c-1];

    double insert_x_score = std::log10(hmm_.GetTransitionMatrix().GetTransitionProbability(INSERT_X, MATCH)) + v_x_[r-1][c-1];
    if (insert_x_score > best_score) {
        best_score = insert_x_score;
    }

    double insert_y_score = std::log10(hmm_.GetTransitionMatrix().GetTransitionProbability(INSERT_Y, MATCH)) + v_y_[r-1][c-1];
    if (insert_y_score > best_score) {
        best_score = insert_y_score;
    }

    NUCLEOTIDE x_i = ConvertCharToNucleotide(hmm_.GetSeq1()[c-1]);
    NUCLEOTIDE y_j = ConvertCharToNucleotide(hmm_.GetSeq2()[r-1]);

    auto transition_score = std::log10(hmm_.GetEmissionMatrix().GetEmissionProbability(x_i, y_j));

    v_match_[r][c] = transition_score + best_score;
}

void ViterbiCalculator::CalculateInsertXProbability(size_t r, size_t c) {
    double best_score = std::log10(hmm_.GetTransitionMatrix().GetTransitionProbability(MATCH, INSERT_X)) + v_match_[r][c-1];

    double insert_x_score = std::log10(hmm_.GetTransitionMatrix().GetTransitionProbability(INSERT_X, INSERT_X)) + v_x_[r][c-1];
    if (insert_x_score > best_score) {
        best_score = insert_x_score;
    }

    v_x_[r][c] = std::log10(hmm_.GetGapEmissionProbability()) + best_score;
}

void ViterbiCalculator::CalculateInsertYProbability(size_t r, size_t c) {
    double best_score = std::log10(hmm_.GetTransitionMatrix().GetTransitionProbability(MATCH, INSERT_Y)) + v_match_[r-1][c];
    double insert_y_score = std::log10(hmm_.GetTransitionMatrix().GetTransitionProbability(INSERT_Y, INSERT_Y)) + v_y_[r-1][c];
    if (insert_y_score > best_score) {
        best_score = insert_y_score;
    }

    v_y_[r][c] = std::log10(hmm_.GetGapEmissionProbability()) + best_score;
}
} // namespace PairHMM
