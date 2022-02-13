#include "simulator.hpp"
#include "gtest/gtest.h"
#include <exception>

TEST(SimulatorTests, ConstructorSafe) {
    EXPECT_NO_THROW(Simulator s{ 1024 });
    EXPECT_THROW(Simulator s{ 0 }, std::runtime_error);
    EXPECT_THROW(Simulator s{ 1 }, std::runtime_error);
}