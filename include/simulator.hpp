#pragma once
#include <filesystem>
#include <vector>
#include "fmt/core.h"
#include <mpi.h>

struct counter {
    unsigned count = 0;
    double time = 0;
    std::string methodName;
    unsigned bytesMoved;

    void setCounter(std::string _methodName, unsigned _bytesMoved) {
        methodName = _methodName;
        bytesMoved = _bytesMoved;
    }

    double getBandwidth() {
        return (bytesMoved/time/1000) * count;
    }

    void printCounter() {
        fmt::print("{:<25} {:<25} {:<25} {:<25.5f} \n", methodName, time, count, getBandwidth());
    }
};


class Simulator {
    static inline bool printing = false;
    counter counterSolveU;
    counter counterSolveV;
    counter counterSolveP;
    counter counterBoundaryU;
    counter counterBoundaryV;
    counter counterBoundaryP;
public:
    static void setPrinting(bool toPrint);

    typedef unsigned SizeType;
    typedef double FloatType;

    // these should be a private section (until constructor), but for testing and benchmarking we keep it public
    SizeType grid;
    const FloatType dx, dy, dt;
    std::vector<FloatType> u, un, v, vn, p, pn, m;

    // helper functions for constructor
    void initU();
    void initV();
    void initP();
    void initCounters();

    // helper functions for run
    void solveUMomentum(const FloatType Re);
    void applyBoundaryU();

    void solveVMomentum(const FloatType Re);
    void applyBoundaryV();

    void solveContinuityEquationP(const FloatType delta);
    void applyBoundaryP();

    FloatType calculateError();

    void iterateU();
    void iterateV();
    void iterateP();

    void deallocate();

public:
    Simulator(SizeType gridP);
    void run(const FloatType delta, const FloatType Re, unsigned maxSteps = std::numeric_limits<unsigned>::max());
    ~Simulator();
};
