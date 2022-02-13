#pragma once
#include <vector>

class Simulator {
    const unsigned grid;
    const double dx, dy, dt;
    std::vector<double> u, un, uc, v, vn, vc, p, pn, pc, m;

public:
    Simulator(unsigned gridP);
    void run(double delta, double error, double Re);
};