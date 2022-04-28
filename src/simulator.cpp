#include "simulator.hpp"
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <chrono>
#include <mpi.h>

void Simulator::setPrinting(bool toPrint) { printing = toPrint && my_rank == 0; }

void Simulator::initU() {
    for (SizeType i = 0; i <= (grid - 1); i++) {
        if (get_coords()[0] == 0) u[(i) * (grid + 1) + grid] = 1.0;
        if (get_coords()[0] == 0) u[(i) * (grid + 1) + grid - 1] = 1.0;
        for (SizeType j = 0; j < (grid - 1); j++) {
            u[(i) * (grid + 1) + j] = 0.0;
        }
    }
}

void Simulator::initV() {
    for (SizeType i = 0; i <= (grid); i++) {
        for (SizeType j = 0; j <= (grid); j++) {
            v[(i)*(grid + 1) + j] = 0.0;
        }
    }
}

void Simulator::initP() {
    for (SizeType i = 0; i <= (grid); i++) {
        for (SizeType j = 0; j <= (grid); j++) {
            p[(i) * (grid + 1) + j] = 1.0;
        }
    }
}

void Simulator::solveUMomentum(const FloatType Re) {
    auto t1 = std::chrono::high_resolution_clock::now();
    exchangeHalo(grid, grid, &u[0]);
    exchangeHalo(grid, grid, &v[0]);
    exchangeHalo(grid, grid, &p[0]);
    #pragma omp parallel for
    for (SizeType i = 1; i <= (grid - 2); i++) {
        for (SizeType j = 1; j <= (grid - 1); j++) {
            un[(i) * (grid + 1) + j] = u[(i) * (grid + 1) + j]
                - dt
                    * ((u[(i + 1) * (grid + 1) + j] * u[(i + 1) * (grid + 1) + j] - u[(i - 1) * (grid + 1) + j] * u[(i - 1) * (grid + 1) + j]) / 2.0 / dx
                    + 0.25 * ((u[(i) * (grid + 1) + j] + u[(i) * (grid + 1) + j + 1]) * (v[(i)*(grid + 1) + j] + v[(i + 1) * (grid + 1) + j])
                            - (u[(i) * (grid + 1) + j] + u[(i) * (grid + 1) + j - 1]) * (v[(i + 1) * (grid + 1) + j - 1] + v[(i)*(grid + 1) + j - 1])) / dy)
                    - dt / dx * (p[(i + 1) * (grid + 1) + j] - p[(i) * (grid + 1) + j]) + dt * 1.0 / Re
                    * ((u[(i + 1) * (grid + 1) + j] - 2.0 * u[(i) * (grid + 1) + j] + u[(i - 1) * (grid + 1) + j]) / dx / dx
                     + (u[(i) * (grid + 1) + j + 1] - 2.0 * u[(i) * (grid + 1) + j] + u[(i) * (grid + 1) + j - 1]) / dy / dy);
        }//writing u, reading u, v, p neighbours -> 3 ghost exchanges before loop
        
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    counterSolveU.time += static_cast<FloatType>(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count());
    counterSolveU.count++;
}

void Simulator::applyBoundaryU() {//no reading -> no halo exchange
    auto t1 = std::chrono::high_resolution_clock::now();
    for (SizeType j = 1; j <= (grid - 1); j++) {
        un[(0) * (grid + 1) + j] = 0.0;
        un[(grid - 1) * (grid + 1) + j] = 0.0;
    }

    for (SizeType i = 0; i <= (grid - 1); i++) {
        un[(i) * (grid + 1) + 0] = -un[(i) * (grid + 1) + 1];
        un[(i) * (grid + 1) + grid] = 2 - un[(i) * (grid + 1) + grid - 1];
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    counterBoundaryU.time += static_cast<FloatType>(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count());
    counterBoundaryU.count++;
}

void Simulator::solveVMomentum(const FloatType Re) {
    auto t1 = std::chrono::high_resolution_clock::now();
    exchangeHalo(grid, grid, &u[0]);
    exchangeHalo(grid, grid, &v[0]);
    exchangeHalo(grid, grid, &p[0]);
    #pragma omp parallel for collapse(2)
    for (SizeType i = 1; i <= (grid - 1); i++) {
        for (SizeType j = 1; j <= (grid - 2); j++) {
            vn[(i)*(grid + 1) + j] = v[(i)*(grid + 1) + j]
                - dt * (0.25 * ((u[(i) * (grid + 1) + j] + u[(i) * (grid + 1) + j + 1]) * (v[(i)*(grid + 1) + j] + v[(i + 1) * (grid + 1) + j])
                              - (u[(i - 1) * (grid + 1) + j] + u[(i - 1) * (grid + 1) + j + 1]) * (v[(i)*(grid + 1) + j] + v[(i - 1) * (grid + 1) + j])) / dx
                              + (v[(i)*(grid + 1) + j + 1] * v[(i)*(grid + 1) + j + 1] - v[(i)*(grid + 1) + j - 1] * v[(i)*(grid + 1) + j - 1]) / 2.0 / dy)
                              - dt / dy * (p[(i) * (grid + 1) + j + 1] - p[(i) * (grid + 1) + j]) + dt * 1.0 / Re
                              * ((v[(i + 1) * (grid + 1) + j] - 2.0 * v[(i)*(grid + 1) + j] + v[(i - 1) * (grid + 1) + j]) / dx / dx
                              + (v[(i)*(grid + 1) + j + 1] - 2.0 * v[(i)*(grid + 1) + j] + v[(i)*(grid + 1) + j - 1]) / dy / dy);
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    counterSolveV.time += static_cast<FloatType>(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count());
    counterSolveV.count++;
}

void Simulator::applyBoundaryV() {
    auto t1 = std::chrono::high_resolution_clock::now();
    for (SizeType j = 1; j <= (grid - 2); j++) {
        vn[(0) * (grid + 1) + j] = -vn[(1) * (grid + 1) + j];
        vn[(grid)*(grid + 1) + j] = -vn[(grid - 1) * (grid + 1) + j];
    }

    for (SizeType i = 0; i <= (grid); i++) {
        vn[(i)*(grid + 1) + 0] = 0.0;
        vn[(i)*(grid + 1) + grid - 1] = 0.0;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    counterBoundaryV.time += static_cast<FloatType>(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count());
    counterBoundaryV.count++;
}

void Simulator::solveContinuityEquationP(const FloatType delta) {
    auto t1 = std::chrono::high_resolution_clock::now();
    exchangeHalo(grid, grid, &un[0]);
    exchangeHalo(grid, grid, &vn[0]);
    exchangeHalo(grid, grid, &p[0]);
    #pragma omp parallel for
    for (SizeType i = 1; i <= (grid - 1); i++) {
        for (SizeType j = 1; j <= (grid - 1); j++) {
            pn[(i) * (grid + 1) + j] = p[(i) * (grid + 1) + j]
                - dt * delta * ((un[(i) * (grid + 1) + j] - un[(i - 1) * (grid + 1) + j]) / dx + (vn[(i)*(grid + 1) + j] - vn[(i)*(grid + 1) + j - 1]) / dy);
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    counterSolveP.time += static_cast<FloatType>(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count());
    counterSolveP.count++;
}

void Simulator::applyBoundaryP() {
    auto t1 = std::chrono::high_resolution_clock::now();
    for (SizeType i = 1; i <= (grid - 1); i++) {
        pn[(i) * (grid + 1) + 0] = pn[(i) * (grid + 1) + 1];
        pn[(i) * (grid + 1) + grid] = pn[(i) * (grid + 1) + grid - 1];
    }

    for (SizeType j = 0; j <= (grid); j++) {
        pn[(0) * (grid + 1) + j] = pn[(1) * (grid + 1) + j];
        pn[(grid) * (grid + 1) + j] = pn[(grid - 1) * (grid + 1) + j];
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    counterBoundaryP.time += static_cast<FloatType>(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count());
    counterBoundaryP.count++;
}

Simulator::FloatType Simulator::calculateError() {//write m, read un, vn -> halo exchange
    FloatType error = 0.0;
    exchangeHalo(grid, grid, &un[0]);
    exchangeHalo(grid, grid, &vn[0]);
    #pragma omp parallel for reduction (+:error)
    for (SizeType i = 1; i <= (grid - 1); i++) {
        for (SizeType j = 1; j <= (grid - 1); j++) {
            m[(i) * (grid + 1) + j] =
                ((un[(i) * (grid + 1) + j] - un[(i - 1) * (grid + 1) + j]) / dx + (vn[(i)*(grid + 1) + j] - vn[(i)*(grid + 1) + j - 1]) / dy);
            error += fabs(m[(i) * (grid + 1) + j]);
        }
    }
    double error_g;
    MPI_Allreduce(&error, &error_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    error = error_g;
    return error;
}

void Simulator::iterateU() {
    std::swap(u, un);
    // for (SizeType i = 0; i <= (grid - 1); i++) {
    //     for (SizeType j = 0; j <= (grid); j++) {
    //         u[(i) * (grid + 1) + j] = un[(i) * (grid + 1) + j];
    //     }
    // }
}

void Simulator::iterateV() {
    std::swap(v, vn);
    // for (SizeType i = 0; i <= (grid); i++) {
    //     for (SizeType j = 0; j <= (grid - 1); j++) {
    //         v[(i)*(grid + 1) + j] = vn[(i)*(grid + 1) + j];
    //     }
    // }
}

void Simulator::iterateP() {
    std::swap(p, pn);
    // for (SizeType i = 0; i <= (grid); i++) {
    //     for (SizeType j = 0; j <= (grid); j++) {
    //         p[(i) * (grid + 1) + j] = pn[(i) * (grid + 1) + j];
    //     }
    // }
}

void Simulator::deallocate() {
    // it doesn't do anything until we use vectors
    // because that deallocates automatically
    // but if we have to use a more raw data structure later it is needed
    // and when the the Tests overwrites some member those might won't deallocate
}

Simulator::Simulator(SizeType gridP)
    : grid([](auto g) {
          if (g <= 1) {
              throw std::runtime_error("Grid is smaller or equal to 1.0, give larger number");
          }
          return g;
      }(gridP)),
      dx(1.0 / static_cast<FloatType>(grid - 1)),
      dy(1.0 / static_cast<FloatType>(grid - 1)),
      dt(0.001 / std::pow(grid / 128.0 * 2.0, 2.0)) {

    u.resize((grid + 2) * (grid + 2));
    un.resize((grid + 2) * (grid + 2));
    v.resize((grid + 2) * (grid + 2));
    vn.resize((grid + 2) * (grid + 2));
    p.resize((grid + 2) * (grid + 2));
    pn.resize((grid + 2) * (grid + 2));
    m.resize((grid + 2) * (grid + 2));

    initU();
    initV();
    initP();
    initCounters();
}

void Simulator::initCounters() {
    unsigned byteMoved = static_cast<unsigned>((grid - 1) * (grid - 1) * sizeof(FloatType));
    unsigned byteMoved2 = static_cast<unsigned>((grid - 1) * (grid - 2) * sizeof(FloatType));
    counterSolveU.setCounter("solveUMomentum", byteMoved2 * 4);
    counterSolveV.setCounter("solveVMomentum", byteMoved2 * 4);
    counterSolveP.setCounter("solveContinuityEquationP", byteMoved * 4);
    counterBoundaryU.setCounter("applyBoundaryU", static_cast<unsigned>((grid - 1) * sizeof(FloatType) * 5));
    counterBoundaryV.setCounter("applyBoundaryV", static_cast<unsigned>(((grid - 2) * 4 + grid * 2) * sizeof(FloatType)));
    counterBoundaryP.setCounter("applyBoundaryP", static_cast<unsigned>(((grid - 1) * 4 + grid * 4) * sizeof(FloatType)));
}

void Simulator::run(const FloatType delta, const FloatType Re, unsigned maxSteps) {
    if (printing) {
        fmt::print("Running simulation with delta: {}, Re: {}\n", delta, Re);
    }
    auto error = std::numeric_limits<FloatType>::max();
    unsigned step = 1;
    while (error > 0.00000001 && step <= maxSteps) {
        solveUMomentum(Re);
        applyBoundaryU();

        solveVMomentum(Re);
        applyBoundaryV();

        solveContinuityEquationP(delta);
        applyBoundaryP();

        error = calculateError();

        if (printing && (step % 100 == 1)) {
            fmt::print("Error is {} for the step {}\n", error, step);
        }

        iterateU();
        iterateV();
        iterateP();
        ++step;
    }
    //Print bandwidth table
    if (printing) {
        fmt::print("{:<25} {:<25} {:<25} {:<25}\n", "Method", "Time (microseconds)", "Count", "Bandwidth (GB/s)");
        counterSolveU.printCounter();
        counterSolveV.printCounter();
        counterSolveP.printCounter();
        counterBoundaryU.printCounter();
        counterBoundaryV.printCounter();
        counterBoundaryP.printCounter();
    }
}

Simulator::~Simulator() { deallocate(); }
