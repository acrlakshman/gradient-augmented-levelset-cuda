/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include <cpu/grid.h>
#include <gtest/gtest.h>

using namespace std;

TEST(CPU, GRID)
{
  GALS::CPU::Grid<double, 3> grid(10);

  EXPECT_TRUE(grid.dimension() == 3);
  EXPECT_TRUE(grid.getNumCells()[0] == 10);
  EXPECT_TRUE(grid.getNumCells()[1] == 1);
  EXPECT_TRUE(grid.getNumCells()[2] == 1);
  EXPECT_TRUE(grid.getPadding() == 1);

  grid.setPadding(2);
  EXPECT_TRUE(grid.getPadding() == 2);
}
