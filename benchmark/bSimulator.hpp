#include "simulator.hpp"
#include "benchmark/benchmark.h"

void benchConstructor(benchmark::State& state) {
    auto grid = static_cast<unsigned>(state.range(0));
    for(auto _ : state) {
        Simulator s{grid};
    }
}