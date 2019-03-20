/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include <cpu/array.h>

#include <gtest/gtest.h>

using namespace std;

TEST(CPU, ARRAY) {
  GALS::CPU::Grid<double, 3> grid(10);
  GALS::CPU::Array<GALS::CPU::Grid<double, 3>> levelset(grid);

  EXPECT_TRUE(levelset.size() == 108);
}
