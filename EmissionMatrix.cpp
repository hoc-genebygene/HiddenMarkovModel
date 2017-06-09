#include "EmissionMatrix.h"

namespace PairHMM {
EmissionMatrix::EmissionMatrix(double p_aa, double p_ac, double p_ag, double p_at,
                               double p_ca, double p_cc, double p_cg, double p_ct,
                               double p_ga, double p_gc, double p_gg, double p_gt,
                               double p_ta, double p_tc, double p_tg, double p_tt) {

    auto SetEmissionProbability = [this](NUCLEOTIDE first_nuc, NUCLEOTIDE second_nuc, double probability) {
        emission_probs_[first_nuc * 4 + second_nuc] = probability;
    };

    SetEmissionProbability(A, A, p_aa);
    SetEmissionProbability(A, C, p_ac);
    SetEmissionProbability(A, G, p_ag);
    SetEmissionProbability(A, T, p_at);
    SetEmissionProbability(C, A, p_ca);
    SetEmissionProbability(C, C, p_cc);
    SetEmissionProbability(C, G, p_cg);
    SetEmissionProbability(C, T, p_ct);
    SetEmissionProbability(G, A, p_ga);
    SetEmissionProbability(G, C, p_gc);
    SetEmissionProbability(G, G, p_gg);
    SetEmissionProbability(G, T, p_gt);
    SetEmissionProbability(T, A, p_ta);
    SetEmissionProbability(T, C, p_tc);
    SetEmissionProbability(T, G, p_tg);
    SetEmissionProbability(T, T, p_tt);
}
} // namespace PairHMM