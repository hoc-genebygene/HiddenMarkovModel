#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stack>
#include <vector>

namespace PairHMM {


namespace {
enum NUCLEOTIDE {
    A,
    C,
    G,
    T
};

enum STATE {
    BEGINHIDDEN1,
    RX1,
    BEGINHIDDEN2,
    RY1,
    BEGINHIDDEN3,
    MATCH,
    INSERT_X,
    INSERT_Y,
    ENDHIDDEN1,
    RX2,
    ENDHIDDEN2,
    RY2,
    ENDHIDDEN3
};

std::string StateToString(STATE state) {
    switch (state) {
        case INSERT_X:
            return "DELETION";
        case INSERT_Y:
            return "INSERTION";
        case MATCH:
            return "MATCH/MISMATCH";
        default:
            throw std::logic_error("Unhandled state in StateToString");
    }
}

NUCLEOTIDE ConvertCharToNucleotide(char c) {
    if (c == 'a' || c == 'A')
        return A;

    if (c == 'c' || c == 'C')
        return C;

    if (c == 'g' || c == 'G')
        return G;

    if (c == 't' || c == 'T')
        return T;

    throw std::runtime_error("Unrecognized nucleotide");
}

} // anonymous namespace

// Matrix class but implemented as linear vector due to performance considerations.
// Note that most of the HMM algorithms access things in row major order. We will keep the mat[r][c] convention but internally organize the vector so that accessing the next element in the column is contiguous.
template <class T>
class Matrix {
    struct MatrixRow; // forward declare this to be used by operator[]

public:
    Matrix(size_t num_rows, size_t num_cols, const T & value) : num_rows_(num_rows), num_cols_(num_cols), vec_(std::vector<T>(num_rows * num_cols, value)) {}

    MatrixRow operator[] (size_t row_index) { return MatrixRow(row_index, *this); }

    size_t GetNumRows() { return num_rows_; }
    size_t GetNumCols() { return num_cols_; }

private:
    struct MatrixRow {
        Matrix & mat;
        size_t row_index;

        MatrixRow(size_t row_index, Matrix & mat) : mat(mat), row_index(row_index) {}

        T& operator[] (size_t col_index) { return mat.vec_[col_index * mat.num_rows_ + row_index]; }
    };

    size_t num_rows_;
    size_t num_cols_;

    std::vector<T> vec_;
};

// Store the EmissionMatrix as a log matrix.
// This allows multiplication to become addition, solving issues with precision from multiplying probabilities repeatedly.
class EmissionMatrix {
public:
    EmissionMatrix(double p_aa, double p_ac, double p_ag, double p_at,
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

    double GetEmissionProbability(NUCLEOTIDE first_nuc, NUCLEOTIDE second_nuc) const { return emission_probs_[first_nuc * 4 + second_nuc]; }

private:
    std::array<double, 16> emission_probs_; // Probability of a pair of nucleotide being emitted (16 of them, AA, AC, AG, AT, CA, CC... etc)
};

class TransitionMatrix {
public:
    TransitionMatrix(double delta, double epsilon, double eta, double tau) {
        for (int k = 0; k < transition_probs_.size(); ++k) {
            transition_probs_[k] = 0.0;
        }


        auto SetTransitionProbability = [this](STATE before_state, STATE after_state, double probability) {
            transition_probs_[before_state * 5 + after_state] = probability;
        };

        SetTransitionProbability(BEGINHIDDEN1, RX1, 1.0 - eta);
        SetTransitionProbability(BEGINHIDDEN1, BEGINHIDDEN2, eta);

        SetTransitionProbability(RX1, RX1, 1.0 - eta);
        SetTransitionProbability(RX1, BEGINHIDDEN2, eta);

        SetTransitionProbability(BEGINHIDDEN2, RY1, 1.0-eta);
        SetTransitionProbability(BEGINHIDDEN2, BEGINHIDDEN3, eta);

        SetTransitionProbability(RY1, RY1, 1.0 - eta);
        SetTransitionProbability(RY1, BEGINHIDDEN3, eta);

        SetTransitionProbability(BEGINHIDDEN3, MATCH, 1.0-2.0*delta-tau);
        SetTransitionProbability(BEGINHIDDEN3, INSERT_X, delta);
        SetTransitionProbability(BEGINHIDDEN3, INSERT_Y, delta);
        SetTransitionProbability(BEGINHIDDEN3, ENDHIDDEN1, tau);

//        SetTransitionProbability(BEGINHIDDEN3, MATCH, 1.0 - tau);
//        SetTransitionProbability(BEGINHIDDEN3, ENDHIDDEN1, tau);

        SetTransitionProbability(MATCH, MATCH, 1.0 - 2.0 * delta - tau);
        SetTransitionProbability(MATCH, INSERT_X, delta);
        SetTransitionProbability(MATCH, INSERT_Y, delta);
        SetTransitionProbability(MATCH, ENDHIDDEN1, tau);

        SetTransitionProbability(INSERT_X, MATCH, 1.0 - epsilon - tau);
        SetTransitionProbability(INSERT_X, INSERT_X, epsilon);
        SetTransitionProbability(INSERT_X, ENDHIDDEN1, tau);

//        SetTransitionProbability(INSERT_X, MATCH, 1.0 - epsilon);
//        SetTransitionProbability(INSERT_X, INSERT_X, epsilon);


        SetTransitionProbability(INSERT_Y, MATCH, 1.0 - epsilon - tau);
        SetTransitionProbability(INSERT_Y, INSERT_Y, epsilon);
        SetTransitionProbability(INSERT_Y, ENDHIDDEN1, tau);

//        SetTransitionProbability(INSERT_Y, MATCH, 1.0 - epsilon);
//        SetTransitionProbability(INSERT_Y, INSERT_Y, epsilon);


        SetTransitionProbability(ENDHIDDEN1, RX2, 1.0 - eta);
        SetTransitionProbability(ENDHIDDEN1, ENDHIDDEN2, eta);

        SetTransitionProbability(RX2, RX2, 1.0 - eta);
        SetTransitionProbability(RX2, ENDHIDDEN2, eta);

        SetTransitionProbability(ENDHIDDEN2, RY2, 1.0 - eta);
        SetTransitionProbability(ENDHIDDEN2, ENDHIDDEN3, eta);

        SetTransitionProbability(RY2, RY2, 1.0 - eta);
        SetTransitionProbability(RY2, ENDHIDDEN3, eta);
    }

    double GetTransitionProbability(STATE before_state, STATE after_state) const {
        return transition_probs_[before_state * 5 + after_state];
    }

private:
    std::array<double, 169>  transition_probs_;

};

class PairHMM {
public:
    PairHMM(EmissionMatrix pair_emission_probabilities, TransitionMatrix transition_probabilities, const std::string & seq1, const std::string & seq2) : pair_emission_probabilities_(std::move(pair_emission_probabilities)), transition_probabilities_(std::move(transition_probabilities)), seq1_(seq1), seq2_(seq2)
    {
        // TODO: Make sure rows and columns of emission probabilities sum to 1
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

class ViterbiPathCalculator {
public:
    ViterbiPathCalculator(PairHMM hmm) :
        hmm_(std::move(hmm)),
        num_rows_(hmm_.GetSeq2().size() + 1),
        num_cols_(hmm_.GetSeq1().size() + 1),
        v_b1_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_x1_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_b2_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_y1_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_b3_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_match_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_x_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_y_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_e1_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_x2_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_e2_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_y2_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_e3_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),

        backtrack_b1_(Matrix<STATE>(num_rows_, num_cols_, BEGINHIDDEN1)),
        backtrack_x1_(Matrix<STATE>(num_rows_, num_cols_, BEGINHIDDEN1)),
        backtrack_b2_(Matrix<STATE>(num_rows_, num_cols_, BEGINHIDDEN1)),
        backtrack_y1_(Matrix<STATE>(num_rows_, num_cols_, BEGINHIDDEN2)),
        backtrack_b3_(Matrix<STATE>(num_rows_, num_cols_, BEGINHIDDEN2)),
        backtrack_match_(Matrix<STATE>(num_rows_, num_cols_, MATCH)),
        backtrack_x_(Matrix<STATE>(num_rows_, num_cols_, BEGINHIDDEN3)),
        backtrack_y_(Matrix<STATE>(num_rows_, num_cols_, BEGINHIDDEN3)),
        backtrack_e1_(Matrix<STATE>(num_rows_, num_cols_, BEGINHIDDEN3)),
        backtrack_x2_(Matrix<STATE>(num_rows_, num_cols_, ENDHIDDEN1)),
        backtrack_e2_(Matrix<STATE>(num_rows_, num_cols_, ENDHIDDEN1)),
        backtrack_y2_(Matrix<STATE>(num_rows_, num_cols_, ENDHIDDEN2)),
        backtrack_e3_(Matrix<STATE>(num_rows_, num_cols_, ENDHIDDEN2))
    {
        SetupInitialConditions();
    }

    void Calculate() {
        // Traverse column by column
        for (size_t c = 1; c < num_cols_; ++c) {
            for (size_t r = 1; r < num_rows_; ++r) {
//                CalculateBeginHidden1Probability(r, c);
//                CalculateX1Probability(r, c);
//                CalculateBeginHidden2Probability(r, c);
//                CalculateY1Probability(r, c);
//                CalculateBeginHidden3Probability(r, c);
                CalculateMatchProbability(r, c);
                CalculateInsertXProbability(r, c);
                CalculateInsertYProbability(r, c);
//                CalculateEndHidden1Probability(r, c);
//                CalculateX2Probability(r, c);
//                CalculateEndHidden2Probability(r, c);
//                CalculateY2Probability(r, c);
//                CalculateEndHidden3Probability(r, c);
            }
        }

        std::cout
//            << std::setiosflags(std::ios::scientific)
            << std::setprecision(3)
//            << std::setw(2)
        ;

        std::cout << "Match matrix: " << std::endl;
        for (size_t r = 0; r < num_rows_; ++r) {
            for (size_t c = 0; c < num_cols_; ++c) {
                std::cout << std::pow(10, v_match_[r][c]) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "InsertX matrix: " << std::endl;
        for (size_t r = 0; r < num_rows_; ++r) {
            for (size_t c = 0; c < num_cols_; ++c) {
                std::cout << std::pow(10, v_x_[r][c]) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "InsertY matrix: " << std::endl;
        for (size_t r = 0; r < num_rows_; ++r) {
            for (size_t c = 0; c < num_cols_; ++c) {
                std::cout << std::pow(10, v_y_[r][c]) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        Backtrack();
    }

private:
    void SetupInitialConditions() {
//        for (int r = 0; r < num_rows_; ++r)
//            v_match_[r][0] = -2.0*std::log10(eta());
//
//        for (int c = 1; c < num_cols_; ++c)
//            v_match_[0][c] = -2.0*std::log10(eta());

        v_match_[0][0] = -2.0*std::log10(eta());
    }

    double delta() {
        return hmm_.GetTransitionMatrix().GetTransitionProbability(MATCH, INSERT_X);
    }

    double epsilon() {
        return hmm_.GetTransitionMatrix().GetTransitionProbability(INSERT_X, INSERT_X);
    }

    double tau() {
        return hmm_.GetTransitionMatrix().GetTransitionProbability(MATCH, ENDHIDDEN1);
    }

    double eta() {
        return hmm_.GetTransitionMatrix().GetTransitionProbability(BEGINHIDDEN1, BEGINHIDDEN2);
    }

    double d() {
        return
            -(
                    std::log10(delta())
                +   std::log10(1.0 - epsilon() - tau())
                -   std::log10(1.0 - eta())
                -   std::log10(1.0 - 2 * delta() - tau())
            );
    }

    double e() {
        return -(std::log10(epsilon()) - std::log10(1.0 - eta()));
    }

    double c() {
        return std::log10(1.0 - 2.0 * delta() - tau()) - std::log10(1.0 - epsilon() - tau());
    }

    double s(NUCLEOTIDE x_i, NUCLEOTIDE y_j) {
//        return
//                std::log10(hmm_.GetEmissionMatrix().GetEmissionProbability(x_i, y_j))
//            -   std::log10(2.0 * hmm_.GetGapEmissionProbability())
//            +   std::log10(1.0 - 2.0 * delta() - tau())
//            -   std::log10(std::pow(1.0 - eta(), 2));

        double a = std::log10(hmm_.GetEmissionMatrix().GetEmissionProbability(x_i, y_j));
        double b = -std::log10(2.0 * hmm_.GetGapEmissionProbability());
        double c = std::log10(1.0 - 2.0 * delta() - tau());
        double d = -std::log10(std::pow(1.0 - eta(), 2));

        return a - b + c - d;
    }

    void CalculateX1Probability(size_t r, size_t c) {
    }

    void CalculateBeginHidden2Probability(size_t r, size_t c) {
    }

    void CalculateY1Probability(size_t r, size_t c) {
    }

    void CalculateBeginHidden3Probability(size_t r, size_t c) {
    }

    void CalculateMatchProbability(size_t r, size_t c) {
        NUCLEOTIDE x_i = ConvertCharToNucleotide(hmm_.GetSeq1()[c-1]);
        NUCLEOTIDE y_j = ConvertCharToNucleotide(hmm_.GetSeq2()[r-1]);
        double substitution_score = s(x_i, y_j);

        double best_score = v_match_[r-1][c-1];
        STATE best_previous_state = MATCH;

        double insert_x_score = v_x_[r-1][c-1];
        if (insert_x_score > best_score) {
            best_score = insert_x_score;
            best_previous_state = INSERT_X;
        }

        double insert_y_score = v_y_[r-1][c-1];
        if (insert_y_score > best_score) {
            best_score = insert_y_score;
            best_previous_state = INSERT_Y;
        }

        v_match_[r][c] = best_score + substitution_score;
        backtrack_match_[r][c] = best_previous_state;
    }

    void CalculateInsertXProbability(size_t r, size_t c) {
        double best_score = v_match_[r][c-1] - d();

        STATE best_previous_state = MATCH;

        double insert_x_score = v_x_[r][c-1] - e();
        if (insert_x_score > best_score) {
            best_score = insert_x_score;
            best_previous_state = INSERT_X;
        }

        v_x_[r][c] = best_score;
        backtrack_x_[r][c] = best_previous_state;
    }

    void CalculateInsertYProbability(size_t r, size_t c) {
        double best_score = v_match_[r-1][c] - d();
        STATE best_previous_state = MATCH;

        double insert_y_score = v_y_[r-1][c] - e();
        if (insert_y_score > best_score) {
            best_score = insert_y_score;
            best_previous_state = INSERT_Y;
        }

        v_y_[r][c] = best_score;
        backtrack_y_[r][c] = best_previous_state;
    }

    void CalculateEndHidden1Probability(size_t r, size_t c) {
    }

    void CalculateX2Probability(size_t r, size_t c) {
    }

    void CalculateEndHidden2Probability(size_t r, size_t c) {
    }

    void CalculateY2Probability(size_t r, size_t c) {
    }

    void CalculateEndHidden3Probability(size_t r, size_t c) {
    }

    void Backtrack() {
        // Scan last column for best prob...
        size_t best_col_index = num_cols_ - 1;
        size_t best_row_index = num_rows_-1;
        double best_prob = v_match_[best_row_index][best_col_index];
        STATE best_previous_state = MATCH;
        STATE best_current_state = MATCH;

        for (int r = num_rows_-1; r < num_rows_; ++r) {
            double match_prob = v_match_[r][best_col_index];
            if (match_prob > best_prob) {
                best_row_index = r;
                best_previous_state = backtrack_match_[r][best_col_index];
                best_current_state = MATCH;
                best_prob = match_prob;
            }

            double x_prob = v_x_[r][best_col_index] + c();
            if (x_prob > best_prob) {
                best_row_index = r;
                best_previous_state = backtrack_x_[r][best_col_index];
                best_current_state = INSERT_X;
                best_prob = x_prob;
            }

            double y_prob = v_y_[r][best_col_index] + c();
            if (y_prob > best_prob) {
                best_row_index = r;
                best_previous_state = backtrack_y_[r][best_col_index];
                best_current_state = INSERT_Y;
                best_prob = y_prob;
            }
        }

        std::stack<STATE> state_stack;

        // Now do actual backtracking
        while (best_row_index != 0 && best_col_index != 0) {

            // Best current state onto stack
            state_stack.push(best_current_state);
            std::cout << "(" << best_row_index << ", " << best_col_index << ")" << std::endl;

            // Update indicies to previous state
            switch (best_previous_state) {
                case MATCH:
                    --best_row_index;
                    --best_col_index;
                    break;
                case INSERT_X:
                    --best_col_index;
                    break;
                case INSERT_Y:
                    --best_row_index;
                    break;
                default:
                    throw std::logic_error("Best previous state impossible value");
            }

            // Update current state to previous state
            best_current_state = best_previous_state;

            // Update previous state from the appropriate backtrack matrix
            switch (best_previous_state) {
                case MATCH:
                    best_previous_state = backtrack_match_[best_row_index][best_col_index];
                    break;
                case INSERT_X:
                    best_previous_state = backtrack_x_[best_row_index][best_col_index];
                    break;
                case INSERT_Y:
                    best_previous_state = backtrack_y_[best_row_index][best_col_index];
                    break;
                default:
                    throw std::logic_error("Best previous state impossible value");
            }
        }

        std::cout << "\n(" << best_row_index << ", " << best_col_index << ")" << std::endl;

        // Print the path
        while (!state_stack.empty()) {
            std::cout << StateToString(state_stack.top()) << std::endl;
            state_stack.pop();
        }
    }

    

private:
    PairHMM hmm_;

    size_t num_rows_;
    size_t num_cols_;

    Matrix<double> v_b1_;
    Matrix<double> v_x1_;
    Matrix<double> v_b2_;
    Matrix<double> v_y1_;
    Matrix<double> v_b3_;
    Matrix<double> v_match_;
    Matrix<double> v_x_;
    Matrix<double> v_y_;
    Matrix<double> v_e1_;
    Matrix<double> v_x2_;
    Matrix<double> v_e2_;
    Matrix<double> v_y2_;
    Matrix<double> v_e3_;

    Matrix<STATE> backtrack_b1_;
    Matrix<STATE> backtrack_x1_;
    Matrix<STATE> backtrack_b2_;
    Matrix<STATE> backtrack_y1_;
    Matrix<STATE> backtrack_b3_;
    Matrix<STATE> backtrack_match_;
    Matrix<STATE> backtrack_x_;
    Matrix<STATE> backtrack_y_;
    Matrix<STATE> backtrack_e1_;
    Matrix<STATE> backtrack_x2_;
    Matrix<STATE> backtrack_e2_;
    Matrix<STATE> backtrack_y2_;
    Matrix<STATE> backtrack_e3_;
};

} // namespace PairHMM

int main() {
    PairHMM::TransitionMatrix trans_mat(
        0.2,    // delta
        0.1,    // epsilon
        0.001,   // eta
        0.001     // tau
    );

//    // Single insertion and deletion only!
//    PairHMM::TransitionMatrix trans_mat(
//                                        0.0, 0.2, 0.2, 0.5, 0.1, // B
//                                        0.0, 0.0, 0.1, 0.8, 0.1, // X
//                                        0.0, 0.1, 0.0, 0.8, 0.1, // Y
//                                        0.0, 0.2, 0.2, 0.5, 0.1, // M
//                                        0.0, 0.0, 0.0, 0.0, 1.0  // H
//                                        );0

//    PairHMM::EmissionMatrix emission_mat(
//        0.5, 0.15, 0.05, 0.3,
//        0.15, 0.5, 0.3, 0.05,
//        0.05, 0.3, 0.5, 0.15,
//        0.3, 0.05, 0.15, 0.5
//    );

    PairHMM::EmissionMatrix emission_mat(
                                         0.997, 0.001, 0.001, 0.001,
                                         0.001, 0.997, 0.001, 0.001,
                                         0.001, 0.001, 0.997, 0.001,
                                         0.001, 0.001, 0.001, 0.997
                                         );

//    PairHMM::EmissionMatrix emission_mat(
//                                         1.0, 0.0, 0.0, 0.0,
//                                         0.0, 1.0, 0.0, 0.0,
//                                         0.0, 0.0, 1.0 ,0.0,
//                                         0.0, 0.0, 0.0, 1.0
//                                         );

    std::string seq1 = "ACAT";
    std::string seq2 = "CAT";

    PairHMM::PairHMM hmm(std::move(emission_mat), std::move(trans_mat), seq1, seq2);

    PairHMM::ViterbiPathCalculator ver_path_calc(hmm);
    ver_path_calc.Calculate();
}
