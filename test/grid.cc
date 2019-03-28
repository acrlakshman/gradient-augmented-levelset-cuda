///////////////////////////////////////////////////////////////////////////////
// Copyright 2019 Lakshman Anumolu, Raunak Bardia.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
// may be used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#include <cpu/grid.h>

#include <gtest/gtest.h>

#include <math.h>

TEST(CPU, GRID_1D)
{
  int num_cells = 10;
  GALS::CPU::Grid<double, 1> grid_1(num_cells, 1, 1);
  GALS::CPU::Grid<double, 1> grid_2(num_cells, 1);
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
    EXPECT_TRUE(fabs(grid.x(i, 0, 0)[0] - x_expect) < 1e-10);
  }

  EXPECT_TRUE(grid.size() == 12);
  EXPECT_TRUE(grid.getMask()[0] == 1);
  EXPECT_TRUE(grid.getMask()[1] == 0);
  EXPECT_TRUE(grid.getMask()[2] == 0);

  EXPECT_TRUE(grid.getIndex(0, 0, 0) == 1);

  EXPECT_TRUE(fabs(grid.dX()[0] - dx) < 1e-10);

  grid.writeToFile();
}

TEST(CPU, GRID_2D)
{
  int n_x = 10, n_y = 10;

  GALS::CPU::Grid<double, 2> grid_1(n_x, n_y, 1);
  GALS::CPU::Grid<double, 2> grid_2(n_x);
  GALS::CPU::Grid<double, 2> grid(n_x, n_y);

  EXPECT_TRUE(grid.dimension() == 2);
  EXPECT_TRUE(grid.getNumCells()[0] == 10);
  EXPECT_TRUE(grid.getNumCells()[1] == 10);
  EXPECT_TRUE(grid.getNumCells()[2] == 1);

  grid.setPadding(2);
  EXPECT_TRUE(grid.getPadding() == 2);

  grid.setPadding(1);
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

  grid.writeToFile();
}

TEST(CPU, GRID_3D)
{
  int num_cells = 10;
  GALS::CPU::Grid<double, 3> grid_1(num_cells, 1, 1);
  GALS::CPU::Grid<double, 3> grid_2(num_cells, 1);
  GALS::CPU::Grid<double, 3> grid(num_cells, num_cells, num_cells);

  EXPECT_TRUE(grid.dimension() == 3);
  EXPECT_TRUE(grid.getNumCells()[0] == 10);
  EXPECT_TRUE(grid.getNumCells()[1] == 10);
  EXPECT_TRUE(grid.getNumCells()[2] == 10);

  grid.setPadding(2);
  EXPECT_TRUE(grid.getPadding() == 2);

  grid.setPadding(1);
  grid.generate(-1, 1, -1, 1, -1, 1);

  double dx = 2. / num_cells;
  double x_min = -1 + (dx * 0.5);

  for (int i = 0; i < num_cells; ++i) {
    double x_expect = x_min + i * dx;
    EXPECT_TRUE(fabs(grid(i, 0, 0)[0] - x_expect) < 1e-10);
    EXPECT_TRUE(fabs(grid.x(i, 0, 0)[0] - x_expect) < 1e-10);
  }

  EXPECT_TRUE(grid.size() == 1728);
  EXPECT_TRUE(grid.getMask()[0] == 1);
  EXPECT_TRUE(grid.getMask()[1] == 1);
  EXPECT_TRUE(grid.getMask()[2] == 1);

  EXPECT_TRUE(grid.getIndex(0, 0, 0) == 157);

  EXPECT_TRUE(fabs(grid.dX()[0] - dx) < 1e-10);

  grid.writeToFile();
}
