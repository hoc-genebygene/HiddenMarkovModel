#include <random>
#include <string>

#include "BackwardCalculator.h"
#include "EmissionMatrix.h"
#include "ForwardCalculator.h"
#include "NUCLEOTIDE.h"
#include "PairHMM.h"
#include "STATE.h"
#include "TransitionMatrix.h"
#include "ViterbiCalculator.h"

int main() {
    PairHMM::EmissionMatrix emission_mat(
        0.5, 0.15, 0.05, 0.3,
        0.15, 0.5, 0.3, 0.05,
        0.05, 0.3, 0.5, 0.15,
        0.3, 0.05, 0.15, 0.5
    );

    PairHMM::TransitionMatrix trans_mat(
        0.2,
        0.1,
        0.1
    );

    std::string seq1 = "TTACG";
    std::string seq2 = "TAG";

    PairHMM::PairHMM hmm(std::move(emission_mat), std::move(trans_mat), seq1, seq2);

    PairHMM::ViterbiCalculator viterbi_calc(hmm);
    viterbi_calc.Calculate();
    viterbi_calc.Print();

    PairHMM::ForwardCalculator forward_calc(hmm);
    forward_calc.Calculate();
    forward_calc.Print();

    PairHMM::BackwardCalculator backward_calc(hmm);
    backward_calc.Calculate();
    backward_calc.Print();
}
