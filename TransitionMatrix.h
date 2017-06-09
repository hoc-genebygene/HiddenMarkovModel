#pragma once

#include <array>

#include "STATE.h"

namespace PairHMM {
class TransitionMatrix {
public:
    TransitionMatrix(double delta, double epsilon, double tau);

    double GetTransitionProbability(STATE before_state, STATE after_state) const {
        return transition_probs_[before_state * 5 + after_state];
    }

private:
    std::array<double, 25> transition_probs_;

};
} // namespace PairHMM
