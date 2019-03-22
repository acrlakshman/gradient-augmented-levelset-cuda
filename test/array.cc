/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include <cpu/array.h>
#include <cpu/vec_n.h>

#include <gtest/gtest.h>

using namespace std;

TEST(CPU, ARRAY) {
  GALS::CPU::Grid<double, 3> grid(10);
  GALS::CPU::Array<GALS::CPU::Grid<double, 3>, double> levelset(grid);

  EXPECT_TRUE(levelset.size() == 108);

  GALS::CPU::Array<GALS::CPU::Grid<double, 3>, GALS::CPU::VecN<double, 2>> velocity(grid);
  EXPECT_TRUE(velocity.size() == 108);
  EXPECT_TRUE(velocity[0].size() == 2);

  for (int i = 0; i < velocity[1].size(); ++i) {
    velocity[1][i] = static_cast<double>(i + 0.5);
  }
  EXPECT_TRUE(velocity[1][1] == 1.5);
}
