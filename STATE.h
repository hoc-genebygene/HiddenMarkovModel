#pragma once

#include <string>

namespace PairHMM {
enum STATE {
    BEGIN,
    MATCH,
    INSERT_X,
    INSERT_Y,
    END
};

std::string StateToString(STATE state);
}
