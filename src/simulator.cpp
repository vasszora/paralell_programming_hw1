#include "simulator.hpp"
#include <cmath>
#include <iostream>
#include <exception>

Simulator::Simulator(unsigned gridP)
    : grid([](unsigned g) {
          if (g <= 1) {
              throw std::runtime_error("Grid is smaller or equal to 1.0, give larger number");
          }
          return g;
      }(gridP)),
      dx(1.0 / static_cast<double>(grid - 1)),
      dy(1.0 / static_cast<double>(grid - 1)),
      dt(0.001 / std::pow(static_cast<double>(grid / 128 * 2), 2.0)),
      u(grid * (grid + 1)),
      un(grid * (grid + 1)),
      uc(grid * grid),
      v((grid + 1) * grid),
      vn((grid + 1) * grid),
      vc(grid * grid),
      p((grid + 1) * (grid + 1)),
      pn((grid + 1) * (grid + 1)),
      pc(grid * grid),
      m((grid + 1) * (grid + 1)) {}

void Simulator::run(double delta, double error, double Re) {
    std::cout << "Running simulation with delta: " << delta << ", error: " << error
              << ", Re: " << Re << "\n";
}