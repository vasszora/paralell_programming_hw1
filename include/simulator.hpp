#pragma once
#include <filesystem>
#include <vector>

class Simulator {
public:
    typedef unsigned SizeType;

    // these should be a private section (until constructor), but for testing and benchmarking we keep it public
    const SizeType      grid;
    const double        dx, dy, dt;
    std::vector<double> u, un, v, vn, p, pn, m;

    // helper functions for constructor
    void initU();
    void initV();
    void initP();

    // helper functions for run
    void solveUMomentum(const double Re);
    void applyBoundaryU();

    void solveVMomentum(const double Re);
    void applyBoundaryV();

    void solveContinuityEquationP(const double delta);
    void applyBoundaryP();

    double calculateError();

    void iterateU();
    void iterateV();
    void iterateP();

    void deallocate();

public:
    Simulator(SizeType gridP);
    void run(const double delta, const double Re);
    ~Simulator();
};