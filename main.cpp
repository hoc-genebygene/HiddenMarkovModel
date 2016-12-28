#include <array>
#include <cmath>
#include <iostream>
#include <stack>
#include <vector>


/* PairHMM implementation, symbols based on loosely on Problem 4.7 in Problems and Solutions in Biological Sequence Analysis by Mark Borodovsky and Svetlana Ekisheva */

namespace PairHMM {


namespace {
enum NUCLEOTIDE {
    A,
    C,
    G,
    T
};

enum STATE {
    BEGIN,
    INSERT_X,
    INSERT_Y,
    MATCH,
    END // Denoted as H in Borovsky, also note that there is a typo for πYO in the text, should be πYH
};

std::string StateToString(STATE state) {
    switch (state) {
        case BEGIN:
            return "BEGIN";
        case INSERT_X:
            return "DELETION";
        case INSERT_Y:
            return "INSERTION";
        case MATCH:
            return "MATCH/MISMATCH";
        case END:
            return "END";
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

using Probability = double;
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

        auto SetEmissionProbability = [this](NUCLEOTIDE first_nuc, NUCLEOTIDE second_nuc, Probability probability) {
            emission_probs_[first_nuc * 4 + second_nuc] = std::log10(probability);
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

    Probability GetEmissionLogProbability(NUCLEOTIDE first_nuc, NUCLEOTIDE second_nuc) const { return emission_probs_[first_nuc * 4 + second_nuc]; }

private:
    std::array<Probability, 16> emission_probs_; // Probability of a pair of nucleotide being emitted (16 of them, AA, AC, AG, AT, CA, CC... etc)
};

class TransitionMatrix {
public:
    TransitionMatrix(double p_bb, double p_bx, double p_by, double p_bm, double p_bh,
                     double p_xb, double p_xx, double p_xy, double p_xm, double p_xh,
                     double p_yb, double p_yx, double p_yy, double p_ym, double p_yh,
                     double p_mb, double p_mx, double p_my, double p_mm, double p_mh,
                     double p_hb, double p_hx, double p_hy, double p_hm, double p_hh) {

        auto SetTransitionProbability = [this](STATE before_state, STATE after_state, Probability probability) {
            transition_probs_[before_state * 5 + after_state] = std::log10(probability);
        };

        SetTransitionProbability(BEGIN, BEGIN, p_mb);
        SetTransitionProbability(BEGIN, INSERT_X, p_bx);
        SetTransitionProbability(BEGIN, INSERT_Y, p_by);
        SetTransitionProbability(BEGIN, MATCH, p_bm);
        SetTransitionProbability(BEGIN, END, p_bh);
        SetTransitionProbability(INSERT_X, BEGIN, p_xb);
        SetTransitionProbability(INSERT_X, INSERT_X, p_xx);
        SetTransitionProbability(INSERT_X, INSERT_Y, p_xy);
        SetTransitionProbability(INSERT_X, MATCH, p_xm);
        SetTransitionProbability(INSERT_X, END, p_xh);
        SetTransitionProbability(INSERT_Y, BEGIN, p_yb);
        SetTransitionProbability(INSERT_Y, INSERT_X, p_yx);
        SetTransitionProbability(INSERT_Y, INSERT_Y, p_yy);
        SetTransitionProbability(INSERT_Y, MATCH, p_ym);
        SetTransitionProbability(INSERT_Y, END, p_yh);
        SetTransitionProbability(MATCH, BEGIN, p_mb);
        SetTransitionProbability(MATCH, INSERT_X, p_mx);
        SetTransitionProbability(MATCH, INSERT_Y, p_my);
        SetTransitionProbability(MATCH, MATCH, p_mm);
        SetTransitionProbability(MATCH, END, p_mh);
        SetTransitionProbability(END, BEGIN, p_hb);
        SetTransitionProbability(END, INSERT_X, p_hx);
        SetTransitionProbability(END, INSERT_Y, p_hy);
        SetTransitionProbability(END, MATCH, p_hm);
        SetTransitionProbability(END, END, p_hh);
    }

    Probability GetTransitionLogProbability(STATE before_state, STATE after_state) const {
        return transition_probs_[before_state * 5 + after_state];
    }

private:
    std::array<Probability, 25>  transition_probs_;

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

    Probability GetGapEmissionLogProbability() { return gap_emission_log_probability_; }

private:
    EmissionMatrix pair_emission_probabilities_;
    Probability gap_emission_log_probability_ = std::log10(0.25); // Gap emission for A,_, C,_, G,_, T,_ should be equal and total to 1 I think?? Anyways can turn into LUT if this turns out to be untrue...

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
        v_match_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_x_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        v_y_(Matrix<double>(num_rows_, num_cols_, std::log10(0.0))),
        backtrack_match_(Matrix<STATE>(num_rows_, num_cols_, MATCH)),
        backtrack_x_(Matrix<STATE>(num_rows_, num_cols_, MATCH)),
        backtrack_y_(Matrix<STATE>(num_rows_, num_cols_, MATCH))
    {
        SetupInitialConditions();
    }

    void Calculate() {
        // Traverse column by column
        for (size_t c = 1; c < num_cols_; ++c) {
            for (size_t r = 1; r < num_rows_; ++r) {
                CalculateViterbiMatchProbability(r, c);
                CalculateViterbiInsertXProbability(r, c);
                CalculateViterbiInsertYProbability(r, c);
            }
        }

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
        v_match_[0][0] = std::log10(1.0);
    }

    void CalculateViterbiMatchProbability(size_t r, size_t c) {
        double best_score = hmm_.GetTransitionMatrix().GetTransitionLogProbability(MATCH, MATCH) + v_match_[r-1][c-1];
        STATE best_previous_state = MATCH;

        double insert_x_score = hmm_.GetTransitionMatrix().GetTransitionLogProbability(INSERT_X, MATCH) + v_x_[r-1][c-1];
        if (insert_x_score > best_score) {
            best_score = insert_x_score;
            best_previous_state = INSERT_X;
        }

        double insert_y_score = hmm_.GetTransitionMatrix().GetTransitionLogProbability(INSERT_Y, MATCH) + v_y_[r-1][c-1];
        if (insert_y_score > best_score) {
            best_score = insert_y_score;
            best_previous_state = INSERT_Y;
        }

        NUCLEOTIDE x_i = ConvertCharToNucleotide(hmm_.GetSeq1()[c-1]);
        NUCLEOTIDE y_j = ConvertCharToNucleotide(hmm_.GetSeq2()[r-1]);
        double v_match_prob = hmm_.GetEmissionMatrix().GetEmissionLogProbability(x_i, y_j);
        v_match_[r][c] = v_match_prob + best_score;
        backtrack_match_[r][c] = best_previous_state;
    }

    void CalculateViterbiInsertXProbability(size_t r, size_t c) {
        double best_score = hmm_.GetTransitionMatrix().GetTransitionLogProbability(MATCH, INSERT_X) + v_match_[r][c-1];
        STATE best_previous_state = MATCH;

        double insert_x_score = hmm_.GetTransitionMatrix().GetTransitionLogProbability(INSERT_X, INSERT_X) + v_x_[r][c-1];
        if (insert_x_score > best_score) {
            best_score = insert_x_score;
            best_previous_state = INSERT_X;
        }

        v_x_[r][c] = best_score + hmm_.GetGapEmissionLogProbability();
        backtrack_x_[r][c] = best_previous_state;
    }

    void CalculateViterbiInsertYProbability(size_t r, size_t c) {
        double best_score = hmm_.GetTransitionMatrix().GetTransitionLogProbability(MATCH, INSERT_Y) + v_match_[r-1][c];
        STATE best_previous_state = MATCH;

        double insert_y_score = hmm_.GetTransitionMatrix().GetTransitionLogProbability(INSERT_Y, INSERT_Y) + v_y_[r-1][c];
        if (insert_y_score > best_score) {
            best_score = insert_y_score;
            best_previous_state = INSERT_Y;
        }

        v_y_[r][c] = best_score + hmm_.GetGapEmissionLogProbability();
        backtrack_y_[r][c] = best_previous_state;
    }

    void Backtrack() {
        // Scan last column for best prob...
        size_t best_col_index = num_cols_ - 1;
        size_t best_row_index = 0;
        double best_prob = v_match_[best_row_index][best_col_index];
        STATE best_previous_state = MATCH;
        STATE best_current_state = MATCH;

        for (int r = 0; r < num_rows_; ++r) {
            double match_prob = v_match_[r][best_col_index];
            if (match_prob > best_prob) {
                best_row_index = r;
                best_previous_state = backtrack_match_[r][best_col_index];
                best_current_state = MATCH;
                best_prob = match_prob;
            }

            double x_prob = v_x_[r][best_col_index];
            if (match_prob > best_prob) {
                best_row_index = r;
                best_previous_state = backtrack_x_[r][best_col_index];
                best_current_state = INSERT_X;
                best_prob = x_prob;
            }

            double y_prob = v_y_[r][best_col_index];
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
            state_stack.push(best_current_state);
            std::cout << "(" << best_row_index << ", " << best_col_index << ")" << std::endl;

            switch (best_current_state) {
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

            best_current_state = best_previous_state;

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

    Matrix<double> v_match_;
    Matrix<double> v_x_;
    Matrix<double> v_y_;

    Matrix<STATE> backtrack_match_;
    Matrix<STATE> backtrack_x_;
    Matrix<STATE> backtrack_y_;
};

} // namespace PairHMM

int main() {
    PairHMM::TransitionMatrix trans_mat(
        0.0, 0.2, 0.2, 0.5, 0.1,
        0.0, 0.1, 0.0, 0.8, 0.1,
        0.0, 0.0, 0.1, 0.8, 0.1,
        0.0, 0.2, 0.2, 0.5, 0.1,
        0.0, 0.0, 0.0, 0.0, 1.0
    );

//    PairHMM::EmissionMatrix emission_mat(
//        0.5, 0.15, 0.05, 0.3,
//        0.15, 0.5, 0.3, 0.05,
//        0.05, 0.3, 0.5, 0.15,
//        0.3, 0.05, 0.15, 0.5
//    );

//    PairHMM::EmissionMatrix emission_mat(
//                                         0.997, 0.001, 0.001, 0.001,
//                                         0.001, 0.997, 0.001, 0.001,
//                                         0.001, 0.001, 0.997, 0.001,
//                                         0.001, 0.001, 0.001, 0.997
//                                         );

    PairHMM::EmissionMatrix emission_mat(
                                         1.0, 0.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0, 0.0,
                                         0.0, 0.0, 1.0 ,0.0,
                                         0.0, 0.0, 0.0, 1.0
                                         );

    std::string seq1 = "CGTCAT";
    std::string seq2 = "CAT";

    PairHMM::PairHMM hmm(std::move(emission_mat), std::move(trans_mat), seq1, seq2);

    PairHMM::ViterbiPathCalculator ver_path_calc(hmm);
    ver_path_calc.Calculate();
}
