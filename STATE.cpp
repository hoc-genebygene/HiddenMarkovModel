#include "STATE.h"

namespace PairHMM {
std::string StateToString(STATE state) {
    switch (state) {
        case INSERT_X:
            return "DELETION";
        case INSERT_Y:
            return "INSERTION";
        case MATCH:
            return "MATCH/MISMATCH";
        case BEGIN:
        case END:

        default:
            throw std::logic_error("Unhandled state in StateToString");
    }
}
}