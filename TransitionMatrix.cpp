#include "TransitionMatrix.h"

namespace PairHMM {
TransitionMatrix::TransitionMatrix(double delta, double epsilon, double tau) {
    for (int k = 0; k < transition_probs_.size(); ++k) {
        transition_probs_[k] = 0.0;
    }

    auto SetTransitionProbability = [this](STATE before_state, STATE after_state, double probability) {
        transition_probs_[before_state * 5 + after_state] = probability;
    };

    SetTransitionProbability(BEGIN, MATCH, 1.0 - 2.0 * delta - tau);
    SetTransitionProbability(BEGIN, INSERT_X, 1.0 - 2.0 * delta - tau);
    SetTransitionProbability(BEGIN, INSERT_Y, 1.0 - 2.0 * delta - tau);

    SetTransitionProbability(MATCH, MATCH, 1.0 - 2.0 * delta - tau);
    SetTransitionProbability(MATCH, INSERT_X, delta);
    SetTransitionProbability(MATCH, INSERT_Y, delta);
    SetTransitionProbability(MATCH, END, tau);

    SetTransitionProbability(INSERT_X, MATCH, 1.0 - epsilon - tau);
    SetTransitionProbability(INSERT_X, INSERT_X, epsilon);
    SetTransitionProbability(INSERT_X, END, tau);

    SetTransitionProbability(INSERT_Y, MATCH, 1.0 - epsilon - tau);
    SetTransitionProbability(INSERT_Y, INSERT_Y, epsilon);
    SetTransitionProbability(INSERT_Y, END, tau);
}
} // namespace PairHMM
