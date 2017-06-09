#include "NUCLEOTIDE.h"

#include <stdexcept>

namespace PairHMM {
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
} // namespace PairHMM
