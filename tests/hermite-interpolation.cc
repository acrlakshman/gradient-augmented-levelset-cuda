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

#include "gals/cpu/interpolate.h"
#include "gals/utilities/array.h"
#include "gals/utilities/utilities.h"

#include <gtest/gtest.h>

#include <math.h>
#include <iostream>

const double xo = 0.;
const double ro = 0.5;

static double oned_levelset(double x, double xo, double ro) { return (x - xo) * (x - xo) - ro * ro; }

static double test_oned(const int nx)
{
  // scalar array on 1D grid.
  GALS::CPU::Grid<double, 1> grid(nx, 1, 1);
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
        levelset(i, j, k) = oned_levelset(grid(i, j, k)[0], xo, ro);
        // std::cout << "grid(" << i << ", " << j << ", " << k << ") = " << grid(i, j, k) << "; levelset(" << i << ", "
        //<< j << ", " << k << ") = " << levelset(i, j, k) << std::endl;
      }

  // Create interpolating points.
  GALS::CPU::Grid<double, 1> grid_interp(23, 1, 1);
  grid_interp.generate(-0.87, 0.87, -1, 1, -1, 1);
  GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::Vec3<double>> x_interp(grid_interp);
  for (int i = 0; i < grid_interp.numCells()[0]; ++i)
    for (int j = 0; j < grid_interp.numCells()[1]; ++j)
      for (int k = 0; k < grid_interp.numCells()[2]; ++k) {
        x_interp(i, j, k)[0] = grid_interp(i, j, k)[0];

        // std::cout << "x_interp(" << i << ", " << j << ", " << k << ") = " << x_interp(i, j, k) << std::endl;
      }

  GALS::CPU::Array<GALS::CPU::Grid<double, 1>, double> levelset_interpolated(grid_interp);
  GALS::CPU::Interpolate<
      double, GALS::CPU::Grid<double, 1>,
      GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 1>>>::compute(x_interp, levelset,
                                                                                 levelset_interpolated);

  // Compute error.
  double l1err = 0.;

  for (int i = 0; i < x_interp.grid().numCells()[0]; ++i)
    for (int j = 0; j < x_interp.grid().numCells()[1]; ++j)
      for (int k = 0; k < x_interp.grid().numCells()[2]; ++k) {
        auto levelset_interp = levelset_interpolated(i, j, k);
        auto levelset_exact = oned_levelset(x_interp(i, j, k)[0], xo, ro);
        l1err += fabs(levelset_interp - levelset_exact);
      }
  l1err /= grid.totalCells();

  return l1err;
}

TEST(CPU, INTERPOLATION_HERMITE_DOUBLE_1D)
{
  std::vector<int> nx = {10, 20, 40, 80, 160};
  std::vector<double> l1err(nx.size());
  std::vector<double> rate(nx.size() - 1);

  for (int i = 0; i < nx.size(); ++i) {
    // std::cout << "------------ nx = " << nx[i] << "------------" << std::endl;
    l1err[i] = test_oned(nx[i]);
  }

  // Debug: output error.
  // for (int i = 0; i < nx.size(); ++i) {
  // std::cout << "l1err[" << i << "]: " << l1err[i] << std::endl;
  //}

  // Compute rate of convergence.
  for (int i = 0; i < nx.size() - 1; ++i) {
    rate[i] = log(l1err[i + 1] / l1err[i]) / log(static_cast<double>(nx[i]) / static_cast<double>(nx[i + 1]));
    std::cout << "rate (" << nx[i] << " -> " << nx[i + 1] << "): " << rate[i] << std::endl;
  }

  // temperary object, for code coverage.
  GALS::CPU::Interpolate<double, GALS::CPU::Grid<double, 1>> interpolate_tmp;

  GALS::CPU::Interpolate<double, GALS::CPU::Grid<double, 1>,
                         GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 1>>>
      hermite_interpolant_1d;

  GALS::INTERPOLATION::ControlPoints<double> control_points =
      GALS::INTERPOLATION::get_control_points<double>(0, 0, 0, 0, 0, false);
}

TEST(CPU, INTERPOLATION_HERMITE_DOUBLE_2D)
{
  const int dim = 2;
  // scalar array on 1D grid.
  GALS::CPU::Grid<double, dim> grid(10, 10, 1);
  grid.generate(-1, 1, -1, 1, -1, 1);
  GALS::CPU::Array<GALS::CPU::Grid<double, dim>, double> levelset(grid);

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

  GALS::CPU::Array<GALS::CPU::Grid<double, dim>, double> levelset_interpolated(grid);
  GALS::CPU::Array<GALS::CPU::Grid<double, dim>, GALS::CPU::Vec3<double>> x_interp(grid);

  GALS::CPU::Interpolate<double, GALS::CPU::Grid<double, dim>> interpolate_tmp;

  GALS::CPU::Interpolate<
      double, GALS::CPU::Grid<double, dim>,
      GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, dim>>>::compute(x_interp, levelset,
                                                                                   levelset_interpolated);
}

TEST(CPU, INTERPOLATION_HERMITE_DOUBLE_3D)
{
  const int dim = 3;
  // scalar array on 1D grid.
  GALS::CPU::Grid<double, dim> grid(10, 10, 10);
  grid.generate(-1, 1, -1, 1, -1, 1);
  GALS::CPU::Array<GALS::CPU::Grid<double, dim>, double> levelset(grid);

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

  GALS::CPU::Array<GALS::CPU::Grid<double, dim>, double> levelset_interpolated(grid);
  GALS::CPU::Array<GALS::CPU::Grid<double, dim>, GALS::CPU::Vec3<double>> x_interp(grid);

  GALS::CPU::Interpolate<double, GALS::CPU::Grid<double, dim>> interpolate_tmp;

  GALS::CPU::Interpolate<
      double, GALS::CPU::Grid<double, dim>,
      GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, dim>>>::compute(x_interp, levelset,
                                                                                   levelset_interpolated);
}
