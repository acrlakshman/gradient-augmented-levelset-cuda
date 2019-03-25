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

#include <math.h>

TEST(CPU, GRID_1D)
{
  int num_cells = 10;
  GALS::CPU::Grid<double, 1> grid(num_cells);

  EXPECT_TRUE(grid.dimension() == 1);
  EXPECT_TRUE(grid.getNumCells()[0] == 10);
  EXPECT_TRUE(grid.getNumCells()[1] == 1);
  EXPECT_TRUE(grid.getNumCells()[2] == 1);

  grid.setPadding(2);
  EXPECT_TRUE(grid.getPadding() == 2);

  grid.setPadding(1);
  grid.generate(-1, 1, -1, 1, -1, 1);

  double dx = 2. / num_cells;
  double x_min = -1 + (dx * 0.5);

  for (int i = 0; i < num_cells; ++i) {
    double x_expect = x_min + i * dx;
    EXPECT_TRUE(fabs(grid(i, 0, 0)[0] - x_expect) < 1e-10);
  }
}

TEST(CPU, GRID_2D)
{
  int n_x = 10, n_y = 10;

  GALS::CPU::Grid<double, 2> grid(n_x, n_y);

  EXPECT_TRUE(grid.dimension() == 2);

  grid.generate(-1, 1, -1, 1, -1, 1);

  double dx = 2. / n_x, dy = 2. / n_y;
  double x_min = -1 + (dx * 0.5);
  double y_min = -1 + (dy * 0.5);

  for (int i = 0; i < n_x; ++i)
    for (int j = 0; j < n_y; ++j) {
      double x_expect = x_min + i * dx;
      double y_expect = y_min + j * dy;
      EXPECT_TRUE(fabs(grid(i, j, 0)[0] - x_expect) < 1e-10);
      EXPECT_TRUE(fabs(grid(i, j, 0)[1] - y_expect) < 1e-10);
      EXPECT_TRUE(fabs(grid.x(i, j, 0)[0] - x_expect) < 1e-10);
      EXPECT_TRUE(fabs(grid.x(i, j, 0)[1] - y_expect) < 1e-10);
    }

  EXPECT_TRUE(grid.size() == 144);
  EXPECT_TRUE(grid.getMask()[0] == 1);
  EXPECT_TRUE(grid.getMask()[1] == 1);
  EXPECT_TRUE(grid.getMask()[2] == 0);

  EXPECT_TRUE(grid.getIndex(0, 0, 0) == 13);

  EXPECT_TRUE(fabs(grid.dX()[0] - dx) < 1e-10);
  EXPECT_TRUE(fabs(grid.dX()[1] - dy) < 1e-10);
  EXPECT_TRUE(fabs(grid.dX()[2] - 2.) < 1e-10);
}
