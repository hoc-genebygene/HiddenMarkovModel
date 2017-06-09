#pragma once

#include "EmissionMatrix.h"
#include "NUCLEOTIDE.h"
#include "TransitionMatrix.h"

namespace PairHMM {
class PairHMM {
public:
    PairHMM(EmissionMatrix pair_emission_probabilities, TransitionMatrix transition_probabilities, const std::string & seq1, const std::string & seq2) : pair_emission_probabilities_(std::move(pair_emission_probabilities)), transition_probabilities_(std::move(transition_probabilities)), seq1_(seq1), seq2_(seq2)
    {
    }

    std::string GetSeq1() const { return seq1_; }
    std::string GetSeq2() const { return seq2_; }

    const TransitionMatrix & GetTransitionMatrix() const { return transition_probabilities_; }
    const EmissionMatrix & GetEmissionMatrix() const { return pair_emission_probabilities_; }

    double GetGapEmissionProbability() { return gap_emission_probability_; }

private:
    EmissionMatrix pair_emission_probabilities_;
    double gap_emission_probability_ = 0.25; // Gap emission for A,_, C,_, G,_, T,_ should be equal and total to 1 I think?? Anyways can turn into LUT if this turns out to be untrue...

    TransitionMatrix transition_probabilities_;

    std::string seq1_;
    std::string seq2_;
};

} // namespace PairHMM
