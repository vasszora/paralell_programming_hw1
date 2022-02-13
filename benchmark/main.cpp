#include "benchmark/benchmark.h"
#include "bSimulator.hpp"

// be careful don't run any code which prints

BENCHMARK(benchConstructor)->RangeMultiplier(2)->Range(32, 2<<10);

BENCHMARK_MAIN();