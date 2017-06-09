#pragma once

#include <array>

#include "NUCLEOTIDE.h"

namespace PairHMM {
class EmissionMatrix {
public:
    EmissionMatrix(double p_aa, double p_ac, double p_ag, double p_at,
                   double p_ca, double p_cc, double p_cg, double p_ct,
                   double p_ga, double p_gc, double p_gg, double p_gt,
                   double p_ta, double p_tc, double p_tg, double p_tt);

    double GetEmissionProbability(NUCLEOTIDE first_nuc, NUCLEOTIDE second_nuc) const { return emission_probs_[first_nuc * 4 + second_nuc]; }

private:
    std::array<double, 16> emission_probs_; // Probability of a pair of nucleotide being emitted (16 of them, AA, AC, AG, AT, CA, CC... etc)
};
} // namespace PairHMM
