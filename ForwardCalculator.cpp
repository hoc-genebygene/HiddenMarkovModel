#include "ForwardCalculator.h"

#include <iostream>

namespace PairHMM {
ForwardCalculator::ForwardCalculator(PairHMM hmm) :
    hmm_(std::move(hmm)),
    num_rows_(hmm_.GetSeq2().size() + 1),
    num_cols_(hmm_.GetSeq1().size() + 1),
    f_match_(Matrix<double>(num_rows_, num_cols_, 0.0)),
    f_x_(Matrix<double>(num_rows_, num_cols_, 0.0)),
    f_y_(Matrix<double>(num_rows_, num_cols_, 0.0))
{
    SetupInitialConditions();
}

void ForwardCalculator::Calculate() {
    for (size_t c = 0; c < num_cols_; ++c) {
        for (size_t r = 0; r < num_rows_; ++r) {
            if (r == 0 && c == 0) {
                continue;
            }

            CalculateForwardMatch(r, c);
            CalculateForwardX(r, c);
            CalculateForwardY(r, c);
        }
    }
}

void ForwardCalculator::Print() {
    std::cout << "Match: " << std::endl;
    for (size_t r = 0; r < num_rows_; ++r) {
        for (size_t c = 0; c < num_cols_; ++c) {
            std::cout << f_match_[r][c] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "InsertX: " << std::endl;
    for (size_t r = 0; r < num_rows_; ++r) {
        for (size_t c = 0; c < num_cols_; ++c) {
            std::cout << f_x_[r][c] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "InsertY: " << std::endl;
    for (size_t r = 0; r < num_rows_; ++r) {
        for (size_t c = 0; c < num_cols_; ++c) {
            std::cout << f_y_[r][c] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void ForwardCalculator::CalculateForwardMatch(size_t r, size_t c) {
        if (r == 0 || c == 0) {
            f_match_[r][c] = 0.0;
            return;
        }

        NUCLEOTIDE x_i = ConvertCharToNucleotide(hmm_.GetSeq1()[c-1]);
        NUCLEOTIDE y_j = ConvertCharToNucleotide(hmm_.GetSeq2()[r-1]);

        double score =
            hmm_.GetEmissionMatrix().GetEmissionProbability(x_i, y_j)
            * (
                hmm_.GetTransitionMatrix().GetTransitionProbability(MATCH, MATCH) * f_match_[r-1][c-1]
            +   hmm_.GetTransitionMatrix().GetTransitionProbability(INSERT_X, MATCH) * f_x_[r-1][c-1]
            +   hmm_.GetTransitionMatrix().GetTransitionProbability(INSERT_Y, MATCH) * f_y_[r-1][c-1]
        );

        f_match_[r][c] = score;
    }

void ForwardCalculator::CalculateForwardX(size_t r, size_t c) {
    if (c == 0) {
        f_x_[r][c] = 0.0;
        return;
    }
    f_x_[r][c] =
        hmm_.GetGapEmissionProbability()
        * (
            hmm_.GetTransitionMatrix().GetTransitionProbability(MATCH, INSERT_X) * f_match_[r][c-1]
        +   hmm_.GetTransitionMatrix().GetTransitionProbability(INSERT_X, INSERT_X) * f_x_[r][c-1]
    );
}

void ForwardCalculator::CalculateForwardY(size_t r, size_t c) {
    if (r == 0) {
        f_y_[r][c] = 0.0;
        return;
    }
    f_y_[r][c] =
    hmm_.GetGapEmissionProbability()
        * (
            hmm_.GetTransitionMatrix().GetTransitionProbability(MATCH, INSERT_Y) * f_match_[r-1][c]
        +   hmm_.GetTransitionMatrix().GetTransitionProbability(INSERT_Y, INSERT_Y) * f_y_[r-1][c]
    );
}

void ForwardCalculator::SetupInitialConditions() {
    f_match_[0][0] = 1.0;
}
} // namespace PairHMM
