#pragma once

namespace PairHMM {
enum NUCLEOTIDE {
    A,
    C,
    G,
    T
};

NUCLEOTIDE ConvertCharToNucleotide(char c);
} // namespace PairHMM
