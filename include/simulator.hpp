#pragma once
#include <vector>
#include <filesystem>

class Simulator {
public:
    typedef unsigned SizeType;

private:
    const SizeType grid;
    const double dx, dy, dt;
    std::vector<double> u, un, uc, v, vn, vc, p, pn, pc, m;

    // this should be a private section, but for testing and benchmarking we keep it public
public:
    // helper functions for constructor
    void initU();
    void initV();
    void initP();

    // helper functions for run
    void solveUMomentum(double Re);
    void applyBoundaryU();

    void solveVMomentum(double Re);
    void applyBoundaryV();

    void solveContinuityEquationP(double delta);
    void applyBoundaryP();

    double calculateError();

    void iterateU();
    void iterateV();
    void iterateP();

    void calculatingCentralArrays();

public:
    Simulator(SizeType gridP);
    void run(double delta, double Re);
    void printUVP(std::filesystem::path filePath);
    void printUCentral(std::filesystem::path filePath);
};