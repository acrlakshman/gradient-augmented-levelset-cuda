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

#include <cpu/array.h>
#include <cpu/gradient.h>
#include <cpu/utilities.h>
#include <cpu/vec3.h>
#include <cpu/vec_n.h>

#include <gtest/gtest.h>

#include <math.h>
#include <iostream>

TEST(CPU, GRADIENT_SECOND_ORDER_CENTRAL_DOUBLE_1D)
{
  // scalar array on 1D grid.
  GALS::CPU::Grid<double, 1> grid(10, 1, 1);
  grid.generate(-1, 1, -1, 1, -1, 1);
  GALS::CPU::Array<GALS::CPU::Grid<double, 1>, double> levelset(grid);

  const auto mask = grid.getMask();
  const int pad = grid.getPadding();
  const auto num_cells = grid.numCells();

  int i_min = -pad * mask[0];
  int j_min = -pad * mask[1];
  int k_min = -pad * mask[2];
  int i_max = num_cells[0] + pad * mask[0];
  int j_max = num_cells[1] + pad * mask[1];
  int k_max = num_cells[2] + pad * mask[2];

  for (int i = i_min; i < i_max; ++i)
    for (int j = j_min; j < j_max; ++j)
      for (int k = k_min; k < k_max; ++k) {
        levelset(i, j, k) = i;
        // std::cout << "levelset(" << i << ", " << j << ", " << k << ") = " << levelset(i, j, k) << std::endl;
      }

  GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::Vec3<double>> grad_levelset(grid);

  GALS::CPU::Gradient<double, GALS::CPU::Grid<double, 1>> gradient;
  gradient.compute(levelset, grad_levelset);

  // for (int i = i_min; i < i_max; ++i)
  // for (int j = j_min; j < j_max; ++j)
  // for (int k = k_min; k < k_max; ++k) {
  // std::cout << "grad_levelset(" << i << ", " << j << ", " << k << ") = {" << grad_levelset(i, j, k)[0] << ", "
  //<< grad_levelset(i, j, k)[1] << ", " << grad_levelset(i, j, k)[2] << "}" << std::endl;
  //}

  for (int i = 0; i < num_cells[0]; ++i)
    for (int j = 0; j < num_cells[1]; ++j)
      for (int k = 0; k < num_cells[2]; ++k) EXPECT_TRUE(GALS::is_equal(5., grad_levelset(i, j, k)[0]));
}

TEST(CPU, GRADIENT_SECOND_ORDER_CENTRAL_DOUBLE_2D)
{
  // scalar array on 2D grid.
  GALS::CPU::Grid<double, 2> grid(10, 10, 1);
  grid.generate(-1, 1, -1, 1, -1, 1);
  GALS::CPU::Array<GALS::CPU::Grid<double, 2>, double> levelset(grid);

  const auto mask = grid.getMask();
  const int pad = grid.getPadding();
  const auto num_cells = grid.numCells();

  int i_min = -pad * mask[0];
  int j_min = -pad * mask[1];
  int k_min = -pad * mask[2];
  int i_max = num_cells[0] + pad * mask[0];
  int j_max = num_cells[1] + pad * mask[1];
  int k_max = num_cells[2] + pad * mask[2];

  for (int i = i_min; i < i_max; ++i)
    for (int j = j_min; j < j_max; ++j)
      for (int k = k_min; k < k_max; ++k) {
        levelset(i, j, k) = i;
      }

  GALS::CPU::Array<GALS::CPU::Grid<double, 2>, GALS::CPU::Vec3<double>> grad_levelset(grid);

  GALS::CPU::Gradient<double, GALS::CPU::Grid<double, 2>> gradient;
  gradient.compute(levelset, grad_levelset);

  for (int i = 0; i < num_cells[0]; ++i)
    for (int j = 0; j < num_cells[1]; ++j)
      for (int k = 0; k < num_cells[2]; ++k) EXPECT_TRUE(GALS::is_equal(5., grad_levelset(i, j, k)[0]));
}

TEST(CPU, GRADIENT_SECOND_ORDER_CENTRAL_DOUBLE_3D)
{
  // scalar array on 2D grid.
  GALS::CPU::Grid<double, 3> grid(10, 10, 10);
  grid.generate(-1, 1, -1, 1, -1, 1);
  GALS::CPU::Array<GALS::CPU::Grid<double, 3>, double> levelset(grid);

  const auto mask = grid.getMask();
  const int pad = grid.getPadding();
  const auto num_cells = grid.numCells();

  int i_min = -pad * mask[0];
  int j_min = -pad * mask[1];
  int k_min = -pad * mask[2];
  int i_max = num_cells[0] + pad * mask[0];
  int j_max = num_cells[1] + pad * mask[1];
  int k_max = num_cells[2] + pad * mask[2];

  for (int i = i_min; i < i_max; ++i)
    for (int j = j_min; j < j_max; ++j)
      for (int k = k_min; k < k_max; ++k) {
        levelset(i, j, k) = i;
      }

  GALS::CPU::Array<GALS::CPU::Grid<double, 3>, GALS::CPU::Vec3<double>> grad_levelset(grid);

  GALS::CPU::Gradient<double, GALS::CPU::Grid<double, 3>> gradient;
  gradient.compute(levelset, grad_levelset);

  for (int i = 0; i < num_cells[0]; ++i)
    for (int j = 0; j < num_cells[1]; ++j)
      for (int k = 0; k < num_cells[2]; ++k) EXPECT_TRUE(GALS::is_equal(5., grad_levelset(i, j, k)[0]));
}
