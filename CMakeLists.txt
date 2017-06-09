project(HiddenMarkovModel CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HMM_HDRS
    BackwardCalculator.h
    EmissionMatrix.h
    ForwardCalculator.h
    Matrix.h
    NUCLEOTIDE.h
    PairHMM.h
    STATE.h
    TransitionMatrix.h
    ViterbiCalculator.h
)

set(HMM_SRCS
    BackwardCalculator.cpp
    EmissionMatrix.cpp
    ForwardCalculator.cpp
    NUCLEOTIDE.cpp
    PairHMM.cpp
    STATE.cpp
    TransitionMatrix.cpp
    ViterbiCalculator.cpp
)

add_library(PairHMM ${HMM_HDRS} ${HMM_SRCS})

add_executable(HiddenMarkovModel main.cpp)

target_link_libraries(HiddenMarkovModel PairHMM)

