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
#include "gals/cpu/levelset.h"
#include "gals/utilities/array.h"
#include "gals/utilities/utilities.h"

#include <gtest/gtest.h>

#include <math.h>
#include <iostream>

const double xo = 0.;
const double ro = 0.5;

static double oned_levelset(double x, double xo, double ro) { return exp((x - xo) * (x - xo) - ro * ro); }
static double oned_levelset_grad(double x, double xo, double ro) { return 2. * (x - xo) * oned_levelset(x, xo, ro); }

struct Errors {
  double l1_phi;
  double l1_psi[3];
  Errors() { l1_phi = 0., l1_psi[0] = 0., l1_psi[1] = 0., l1_psi[2] = 0.; }
};

static Errors test_oned(const int nx)
{
  // scalar array on 1D grid.
  GALS::CPU::Grid<double, 1> grid(nx, 1, 1);
  grid.generate(-1, 1, -1, 1, -1, 1);
  GALS::CPU::Levelset<GALS::CPU::Grid<double, 1>, double> levelset(grid);

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
        levelset.phiTm1()(i, j, k) = oned_levelset(grid(i, j, k)[0], xo, ro);
        levelset.phi()(i, j, k) = levelset.phiTm1()(i, j, k);

        levelset.psiTm1()(i, j, k)[0] = oned_levelset_grad(grid(i, j, k)[0], xo, ro);
        levelset.psi()(i, j, k) = levelset.psiTm1()(i, j, k);

        // std::cout << "grid(" << i << ", " << j << ", " << k << ") = " << grid(i, j, k) << "; levelset.phiTm1()(" << i
        //<< ", " << j << ", " << k << ") = " << levelset.phiTm1()(i, j, k) << std::endl;
      }

  // Create interpolating points.
  GALS::CPU::Grid<double, 1> grid_interp(11, 1, 1);
  grid_interp.generate(-0.9, 0.9, -1, 1, -1, 1);
  GALS::CPU::Array<GALS::CPU::Grid<double, 1>, GALS::CPU::Vec3<double>> x_interp(grid_interp);
  for (int i = 0; i < grid_interp.numCells()[0]; ++i)
    for (int j = 0; j < grid_interp.numCells()[1]; ++j)
      for (int k = 0; k < grid_interp.numCells()[2]; ++k) {
        x_interp(i, j, k)[0] = grid_interp(i, j, k)[0];

        // std::cout << "x_interp(" << i << ", " << j << ", " << k << ") = " << x_interp(i, j, k) << std::endl;
      }

  GALS::CPU::Interpolate<double, GALS::CPU::Grid<double, 1>,
                         GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 1>>>::compute(x_interp, levelset);

  // Compute error.
  Errors errs;

  for (int i = 0; i < x_interp.grid().numCells()[0]; ++i)
    for (int j = 0; j < x_interp.grid().numCells()[1]; ++j)
      for (int k = 0; k < x_interp.grid().numCells()[2]; ++k) {
        auto levelset_interp = levelset.phiInterpTm1()(i, j, k);
        auto levelset_exact = oned_levelset(x_interp(i, j, k)[0], xo, ro);

        errs.l1_phi += fabs(levelset_interp - levelset_exact);

        auto psi_interp = levelset.psiInterpTm1()(i, j, k)[0];
        auto psi_exact = oned_levelset_grad(x_interp(i, j, k)[0], xo, ro);

        errs.l1_psi[0] += fabs(psi_interp - psi_exact);

        // std::cout << "x_interp(" << i << ", " << j << ", " << k << ") = " << x_interp(i, j, k)
        //<< "; levelset.phiInterpTm1()(" << i << ", " << j << ", " << k
        //<< ") = " << levelset.phiInterpTm1()(i, j, k) << "; levelset_exact"
        //<< " = " << levelset_exact << std::endl;
      }
  errs.l1_phi /= x_interp.grid().totalCells();
  errs.l1_psi[0] /= x_interp.grid().totalCells();

  return errs;
}

TEST(CPU, INTERPOLATION_HERMITE_DOUBLE_1D)
{
  std::vector<int> nx = {11, 21, 41, 81};

  // These reference values are obtained from matlab.
  std::vector<double> l1_phi_ref = {71.7922473960529e-006, 4.19183043139219e-006, 263.113371645751e-009,
                                    14.4915820757515e-009};
  std::vector<double> l1_psi_ref = {542.528257269213e-006, 133.121515337428e-006, 21.4701778071460e-006,
                                    2.50561543448733e-006};

  std::vector<double> l1_phi(nx.size());
  std::vector<double> l1_psi(nx.size());
  std::vector<double> rate(nx.size() - 1);

  for (int i = 0; i < nx.size(); ++i) {
    auto errs = test_oned(nx[i]);
    l1_phi[i] = errs.l1_phi;
    l1_psi[i] = errs.l1_psi[0];

    EXPECT_TRUE(fabs(l1_phi[i] - l1_phi_ref[i]) < 1e-12);
    EXPECT_TRUE(fabs(l1_psi[i] - l1_psi_ref[i]) < 1e-12);
  }

  // Debug: output error.
  // for (int i = 0; i < nx.size(); ++i) {
  // std::cout << "l1_phi[" << i << "]: " << l1_phi[i] << std::endl;
  //}
  // std::cout << "--------" << std::endl;
  // for (int i = 0; i < nx.size(); ++i) {
  // std::cout << "l1_psi[" << i << "]: " << l1_psi[i] << std::endl;
  //}

  // Compute rate of convergence.
  // for (int i = 0; i < nx.size() - 1; ++i) {
  // rate[i] = log(l1_phi[i + 1] / l1_phi[i]) / log(static_cast<double>(nx[i]) / static_cast<double>(nx[i + 1]));
  // std::cout << "rate (" << nx[i] << " -> " << nx[i + 1] << "): " << rate[i] << std::endl;
  //}

  // temperary object, for code coverage.
  // GALS::CPU::Interpolate<double, GALS::CPU::Grid<double, 1>> interpolate_tmp;

  // GALS::CPU::Interpolate<double, GALS::CPU::Grid<double, 1>,
  // GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 1>>>
  // hermite_interpolant_1d;

  // GALS::INTERPOLATION::ControlPoints<double> control_points =
  // GALS::INTERPOLATION::get_control_points<double>(0, 0, 0, 0, 0, false);
}

/********************* Testing 2D ***********************/

static Errors test_twod(const int nx, const int ny)
{
  auto twod_exp_circle = [](double x, double y, double xo, double yo, double ro) {
    return exp((x - xo) * (x - xo) + (y - yo) * (y - yo) - ro * ro);
  };
  auto twod_exp_circle_grad_x = [](double x, double y, double xo, double yo, double ro) {
    return 2 * (x - xo) * exp((x - xo) * (x - xo) + (y - yo) * (y - yo) - ro * ro);
  };
  auto twod_exp_circle_grad_y = [](double x, double y, double xo, double yo, double ro) {
    return 2 * (y - yo) * exp((x - xo) * (x - xo) + (y - yo) * (y - yo) - ro * ro);
  };
  auto twod_exp_circle_grad_xy = [](double x, double y, double xo, double yo, double ro) {
    return 4 * (x - xo) * (y - yo) * exp((x - xo) * (x - xo) + (y - yo) * (y - yo) - ro * ro);
  };

  // scalar array on 2D grid.
  GALS::CPU::Grid<double, 2> grid(nx, ny, 1);
  grid.generate(-1, 1, -1, 1, -1, 1);
  GALS::CPU::Levelset<GALS::CPU::Grid<double, 2>, double> levelset(grid);

  const auto mask = grid.getMask();
  const int pad = grid.getPadding();
  const auto num_cells = grid.numCells();

  int i_min = -pad * mask[0];
  int j_min = -pad * mask[1];
  int k_min = -pad * mask[2];
  int i_max = num_cells[0] + pad * mask[0];
  int j_max = num_cells[1] + pad * mask[1];
  int k_max = num_cells[2] + pad * mask[2];

  double xo = 0., yo = 0., ro = 0.5;
  for (int i = i_min; i < i_max; ++i)
    for (int j = j_min; j < j_max; ++j)
      for (int k = k_min; k < k_max; ++k) {
        levelset.phiTm1()(i, j, k) = twod_exp_circle(grid(i, j, k)[0], grid(i, j, k)[1], xo, yo, ro);
        levelset.phi()(i, j, k) = levelset.phiTm1()(i, j, k);

        levelset.psiTm1()(i, j, k)[0] = twod_exp_circle_grad_x(grid(i, j, k)[0], grid(i, j, k)[1], xo, yo, ro);
        levelset.psi()(i, j, k)[0] = levelset.psiTm1()(i, j, k)[0];

        levelset.psiTm1()(i, j, k)[1] = twod_exp_circle_grad_y(grid(i, j, k)[0], grid(i, j, k)[1], xo, yo, ro);
        levelset.psi()(i, j, k)[1] = levelset.psiTm1()(i, j, k)[1];

        levelset.phiMixedDerivativesTm1()(i, j, k)[0] =
            twod_exp_circle_grad_xy(grid(i, j, k)[0], grid(i, j, k)[1], xo, yo, ro);
        levelset.phiMixedDerivatives()(i, j, k)[0] = levelset.phiMixedDerivativesTm1()(i, j, k)[0];

        // std::cout << "grid(" << i << ", " << j << ", " << k << ") = " << grid(i, j, k) << "; levelset.phiTm1()(" << i
        //<< ", " << j << ", " << k << ") = " << levelset.phiTm1()(i, j, k) << std::endl;
      }

  // Create interpolating points.
  GALS::CPU::Grid<double, 2> grid_interp(11, 11, 1);
  grid_interp.generate(-0.9, 0.9, -0.9, 0.9, -1, 1);
  GALS::CPU::Array<GALS::CPU::Grid<double, 2>, GALS::CPU::Vec3<double>> x_interp(grid_interp);
  for (int i = 0; i < grid_interp.numCells()[0]; ++i)
    for (int j = 0; j < grid_interp.numCells()[1]; ++j)
      for (int k = 0; k < grid_interp.numCells()[2]; ++k) {
        x_interp(i, j, k)[0] = grid_interp(i, j, k)[0];
        x_interp(i, j, k)[1] = grid_interp(i, j, k)[1];
        // std::cout << "x_interp(" << i << ", " << j << ", " << k << ") = " << x_interp(i, j, k) << std::endl;
      }

  GALS::CPU::Interpolate<double, GALS::CPU::Grid<double, 2>,
                         GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 2>>>::compute(x_interp, levelset);

  // Compute error.
  Errors errs;

  for (int i = 0; i < x_interp.grid().numCells()[0]; ++i)
    for (int j = 0; j < x_interp.grid().numCells()[1]; ++j)
      for (int k = 0; k < x_interp.grid().numCells()[2]; ++k) {
        auto levelset_interp = levelset.phiInterpTm1()(i, j, k);
        auto levelset_exact = twod_exp_circle(x_interp(i, j, k)[0], x_interp(i, j, k)[1], xo, yo, ro);

        errs.l1_phi += fabs(levelset_interp - levelset_exact);

        auto psix_interp = levelset.psiInterpTm1()(i, j, k)[0];
        auto psix_exact = twod_exp_circle_grad_x(x_interp(i, j, k)[0], x_interp(i, j, k)[1], xo, yo, ro);

        auto psiy_interp = levelset.psiInterpTm1()(i, j, k)[1];
        auto psiy_exact = twod_exp_circle_grad_y(x_interp(i, j, k)[0], x_interp(i, j, k)[1], xo, yo, ro);

        errs.l1_psi[0] += fabs(psix_interp - psix_exact);
        errs.l1_psi[1] += fabs(psiy_interp - psiy_exact);

        // std::cout << "x_interp(" << i << ", " << j << ", " << k << ") = " << x_interp(i, j, k) << std::endl
        //<< "\tlevelset.phiInterpTm1()(" << i << ", " << j << ", " << k
        //<< ") = " << levelset.phiInterpTm1()(i, j, k) << "; levelset_exact"
        //<< " = " << levelset_exact << std::endl
        //<< "\tlevelset.psiInterpTm1()(" << i << ", " << j << ", " << k
        //<< ") = " << levelset.psiInterpTm1()(i, j, k) << std::endl
        //<< "\tpsi_exact"
        //<< " = " << psix_exact << ", " << psiy_exact << std::endl;
      }

  errs.l1_phi /= x_interp.grid().totalCells();
  errs.l1_psi[0] /= x_interp.grid().totalCells();
  errs.l1_psi[1] /= x_interp.grid().totalCells();

  return errs;
}

TEST(CPU, INTERPOLATION_HERMITE_DOUBLE_2D)
{
  std::vector<int> nx = {11, 21, 41, 81, 161};
  std::vector<int> ny = {11, 21, 41, 81, 161};

  // These reference values are obtained from matlab.
  std::vector<double> l1_phi_ref = {193.196907657978e-006, 11.2808113414184e-006, 708.076798341652e-009,
                                    38.9989844917323e-009, 2.66600931322065e-009};
  std::vector<double> l1_psix_ref = {715.246551139715e-006, 176.408187501144e-006, 28.8739250246707e-006,
                                     3.37235914274825e-006, 479.180656372292e-009};
  std::vector<double> l1_psiy_ref = {715.246551140005e-006, 176.408187502064e-006, 28.8739250235467e-006,
                                     3.37235914080582e-006, 479.180654183501e-009};

  std::vector<double> l1_phi(nx.size());
  std::vector<double> l1_psix(nx.size());
  std::vector<double> l1_psiy(nx.size());

  for (int i = 0; i < nx.size(); ++i) {
    auto errs = test_twod(nx[i], ny[i]);
    l1_phi[i] = errs.l1_phi;
    l1_psix[i] = errs.l1_psi[0];
    l1_psiy[i] = errs.l1_psi[1];

    EXPECT_TRUE(fabs(l1_phi[i] - l1_phi_ref[i]) < 1e-12);
    EXPECT_TRUE(fabs(l1_psix[i] - l1_psix_ref[i]) < 1e-12);
    EXPECT_TRUE(fabs(l1_psiy[i] - l1_psiy_ref[i]) < 1e-12);
  }

  // Debug
  // std::cout << "Hermite interpolation 2D unit test" << std::endl;
  // for (int i = 0; i < nx.size(); ++i) std::cout << "\tl1_phi[" << i << "]: " << l1_phi[i] << std::endl;
  // for (int i = 0; i < nx.size(); ++i) std::cout << "\tl1_psix[" << i << "]: " << l1_psix[i] << std::endl;
  // for (int i = 0; i < nx.size(); ++i) std::cout << "\tl1_psiy[" << i << "]: " << l1_psiy[i] << std::endl;
}

/********************* Testing 3D ***********************/

static Errors test_threed(const int nx, const int ny, const int nz)
{
  auto threed_exp_circle = [](double x, double y, double z, double xo, double yo, double zo, double ro) {
    return exp((x - xo) * (x - xo) + (y - yo) * (y - yo) + (z - zo) * (z - zo) - ro * ro);
  };
  auto threed_exp_circle_grad_x = [](double x, double y, double z, double xo, double yo, double zo, double ro) {
    return 2 * (x - xo) * exp((x - xo) * (x - xo) + (y - yo) * (y - yo) + (z - zo) * (z - zo) - ro * ro);
  };
  auto threed_exp_circle_grad_y = [](double x, double y, double z, double xo, double yo, double zo, double ro) {
    return 2 * (y - yo) * exp((x - xo) * (x - xo) + (y - yo) * (y - yo) + (z - zo) * (z - zo) - ro * ro);
  };
  auto threed_exp_circle_grad_z = [](double x, double y, double z, double xo, double yo, double zo, double ro) {
    return 2 * (z - zo) * exp((x - xo) * (x - xo) + (y - yo) * (y - yo) + (z - zo) * (z - zo) - ro * ro);
  };
  auto threed_exp_circle_grad_xy = [](double x, double y, double z, double xo, double yo, double zo, double ro) {
    return 4 * (x - xo) * (y - yo) * exp((x - xo) * (x - xo) + (y - yo) * (y - yo) + (z - zo) * (z - zo) - ro * ro);
  };
  auto threed_exp_circle_grad_yz = [](double x, double y, double z, double xo, double yo, double zo, double ro) {
    return 4 * (y - yo) * (z - zo) * exp((x - xo) * (x - xo) + (y - yo) * (y - yo) + (z - zo) * (z - zo) - ro * ro);
  };
  auto threed_exp_circle_grad_zx = [](double x, double y, double z, double xo, double yo, double zo, double ro) {
    return 4 * (x - xo) * (z - zo) * exp((x - xo) * (x - xo) + (y - yo) * (y - yo) + (z - zo) * (z - zo) - ro * ro);
  };
  auto threed_exp_circle_grad_xyz = [](double x, double y, double z, double xo, double yo, double zo, double ro) {
    return 8 * (x - xo) * (y - yo) * (z - zo) *
           exp((x - xo) * (x - xo) + (y - yo) * (y - yo) + (z - zo) * (z - zo) - ro * ro);
  };

  // scalar array on 2D grid.
  GALS::CPU::Grid<double, 3> grid(nx, ny, nz);
  grid.generate(-1, 1, -1, 1, -1, 1);
  GALS::CPU::Levelset<GALS::CPU::Grid<double, 3>, double> levelset(grid);

  const auto mask = grid.getMask();
  const int pad = grid.getPadding();
  const auto num_cells = grid.numCells();

  int i_min = -pad * mask[0];
  int j_min = -pad * mask[1];
  int k_min = -pad * mask[2];
  int i_max = num_cells[0] + pad * mask[0];
  int j_max = num_cells[1] + pad * mask[1];
  int k_max = num_cells[2] + pad * mask[2];

  double xo = 0., yo = 0., zo = 0., ro = 0.5;
  for (int i = i_min; i < i_max; ++i)
    for (int j = j_min; j < j_max; ++j)
      for (int k = k_min; k < k_max; ++k) {
        levelset.phiTm1()(i, j, k) =
            threed_exp_circle(grid(i, j, k)[0], grid(i, j, k)[1], grid(i, j, k)[2], xo, yo, zo, ro);
        levelset.phi()(i, j, k) = levelset.phiTm1()(i, j, k);

        levelset.psiTm1()(i, j, k)[0] =
            threed_exp_circle_grad_x(grid(i, j, k)[0], grid(i, j, k)[1], grid(i, j, k)[2], xo, yo, zo, ro);
        levelset.psi()(i, j, k)[0] = levelset.psiTm1()(i, j, k)[0];

        levelset.psiTm1()(i, j, k)[1] =
            threed_exp_circle_grad_y(grid(i, j, k)[0], grid(i, j, k)[1], grid(i, j, k)[2], xo, yo, zo, ro);
        levelset.psi()(i, j, k)[1] = levelset.psiTm1()(i, j, k)[1];

        levelset.psiTm1()(i, j, k)[2] =
            threed_exp_circle_grad_z(grid(i, j, k)[0], grid(i, j, k)[1], grid(i, j, k)[2], xo, yo, zo, ro);
        levelset.psi()(i, j, k)[2] = levelset.psiTm1()(i, j, k)[2];

        levelset.phiMixedDerivativesTm1()(i, j, k)[0] =
            threed_exp_circle_grad_xy(grid(i, j, k)[0], grid(i, j, k)[1], grid(i, j, k)[2], xo, yo, zo, ro);
        levelset.phiMixedDerivatives()(i, j, k)[0] = levelset.phiMixedDerivativesTm1()(i, j, k)[0];

        levelset.phiMixedDerivativesTm1()(i, j, k)[1] =
            threed_exp_circle_grad_yz(grid(i, j, k)[0], grid(i, j, k)[1], grid(i, j, k)[2], xo, yo, zo, ro);
        levelset.phiMixedDerivatives()(i, j, k)[1] = levelset.phiMixedDerivativesTm1()(i, j, k)[1];

        levelset.phiMixedDerivativesTm1()(i, j, k)[2] =
            threed_exp_circle_grad_zx(grid(i, j, k)[0], grid(i, j, k)[1], grid(i, j, k)[2], xo, yo, zo, ro);
        levelset.phiMixedDerivatives()(i, j, k)[2] = levelset.phiMixedDerivativesTm1()(i, j, k)[2];

        levelset.phiMixedDerivativesTm1()(i, j, k)[3] =
            threed_exp_circle_grad_xyz(grid(i, j, k)[0], grid(i, j, k)[1], grid(i, j, k)[2], xo, yo, zo, ro);
        levelset.phiMixedDerivatives()(i, j, k)[3] = levelset.phiMixedDerivativesTm1()(i, j, k)[3];

        // std::cout << "grid(" << i << ", " << j << ", " << k << ") = " << grid(i, j, k) << "; levelset.phiTm1()(" << i
        //<< ", " << j << ", " << k << ") = " << levelset.phiTm1()(i, j, k) << std::endl;
      }

  // Create interpolating points.
  GALS::CPU::Grid<double, 3> grid_interp(11, 11, 11);
  grid_interp.generate(-0.9, 0.9, -0.9, 0.9, -0.9, 0.9);
  GALS::CPU::Array<GALS::CPU::Grid<double, 3>, GALS::CPU::Vec3<double>> x_interp(grid_interp);
  for (int i = 0; i < grid_interp.numCells()[0]; ++i)
    for (int j = 0; j < grid_interp.numCells()[1]; ++j)
      for (int k = 0; k < grid_interp.numCells()[2]; ++k) {
        x_interp(i, j, k)[0] = grid_interp(i, j, k)[0];
        x_interp(i, j, k)[1] = grid_interp(i, j, k)[1];
        x_interp(i, j, k)[2] = grid_interp(i, j, k)[2];
        // std::cout << "x_interp(" << i << ", " << j << ", " << k << ") = " << x_interp(i, j, k) << std::endl;
      }

  GALS::CPU::Interpolate<double, GALS::CPU::Grid<double, 3>,
                         GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 3>>>::compute(x_interp, levelset);

  // Compute error.
  Errors errs;

  for (int i = 0; i < x_interp.grid().numCells()[0]; ++i)
    for (int j = 0; j < x_interp.grid().numCells()[1]; ++j)
      for (int k = 0; k < x_interp.grid().numCells()[2]; ++k) {
        auto levelset_interp = levelset.phiInterpTm1()(i, j, k);
        auto levelset_exact =
            threed_exp_circle(x_interp(i, j, k)[0], x_interp(i, j, k)[1], x_interp(i, j, k)[2], xo, yo, zo, ro);

        errs.l1_phi += fabs(levelset_interp - levelset_exact);

        auto psix_interp = levelset.psiInterpTm1()(i, j, k)[0];
        auto psix_exact =
            threed_exp_circle_grad_x(x_interp(i, j, k)[0], x_interp(i, j, k)[1], x_interp(i, j, k)[2], xo, yo, zo, ro);

        auto psiy_interp = levelset.psiInterpTm1()(i, j, k)[1];
        auto psiy_exact =
            threed_exp_circle_grad_y(x_interp(i, j, k)[0], x_interp(i, j, k)[1], x_interp(i, j, k)[2], xo, yo, zo, ro);

        auto psiz_interp = levelset.psiInterpTm1()(i, j, k)[2];
        auto psiz_exact =
            threed_exp_circle_grad_z(x_interp(i, j, k)[0], x_interp(i, j, k)[1], x_interp(i, j, k)[2], xo, yo, zo, ro);

        errs.l1_psi[0] += fabs(psix_interp - psix_exact);
        errs.l1_psi[1] += fabs(psiy_interp - psiy_exact);
        errs.l1_psi[2] += fabs(psiz_interp - psiz_exact);

        // std::cout << "x_interp(" << i << ", " << j << ", " << k << ") = " << x_interp(i, j, k) << std::endl
        //<< "\tlevelset.phiInterpTm1()(" << i << ", " << j << ", " << k
        //<< ") = " << levelset.phiInterpTm1()(i, j, k) << "; levelset_exact"
        //<< " = " << levelset_exact << std::endl
        //<< "\tlevelset.psiInterpTm1()(" << i << ", " << j << ", " << k
        //<< ") = " << levelset.psiInterpTm1()(i, j, k) << std::endl
        //<< "\tpsi_exact"
        //<< " = " << psix_exact << ", " << psiy_exact << std::endl;
      }

  errs.l1_phi /= x_interp.grid().totalCells();
  errs.l1_psi[0] /= x_interp.grid().totalCells();
  errs.l1_psi[1] /= x_interp.grid().totalCells();
  errs.l1_psi[2] /= x_interp.grid().totalCells();

  return errs;
}

TEST(CPU, INTERPOLATION_HERMITE_DOUBLE_3D)
{
  std::vector<int> nx = {11, 21, 41, 81, 161};
  std::vector<int> ny = {11, 21, 41, 81, 161};
  std::vector<int> nz = {11, 21, 41, 81, 161};

  // These reference values are obtained from matlab.
  // TODO (obtain reference values from matlab).
  // std::vector<double> l1_phi_ref = {193.196907657978e-006, 11.2808113414184e-006, 708.076798341652e-009,
  // 38.9989844917323e-009, 2.66600931322065e-009};
  // std::vector<double> l1_psix_ref = {715.246551139715e-006, 176.408187501144e-006, 28.8739250246707e-006,
  // 3.37235914274825e-006, 479.180656372292e-009};
  // std::vector<double> l1_psiy_ref = {715.246551140005e-006, 176.408187502064e-006, 28.8739250235467e-006,
  // 3.37235914080582e-006, 479.180654183501e-009};

  std::vector<double> l1_phi(nx.size());
  std::vector<double> l1_psix(nx.size());
  std::vector<double> l1_psiy(nx.size());
  std::vector<double> l1_psiz(nx.size());

  for (int i = 0; i < nx.size(); ++i) {
    auto errs = test_threed(nx[i], ny[i], nz[i]);
    l1_phi[i] = errs.l1_phi;
    l1_psix[i] = errs.l1_psi[0];
    l1_psiy[i] = errs.l1_psi[1];
    l1_psiz[i] = errs.l1_psi[2];

    // EXPECT_TRUE(fabs(l1_phi[i] - l1_phi_ref[i]) < 1e-12);
    // EXPECT_TRUE(fabs(l1_psix[i] - l1_psix_ref[i]) < 1e-12);
    // EXPECT_TRUE(fabs(l1_psiy[i] - l1_psiy_ref[i]) < 1e-12);
  }

  // Debug
  std::cout << "Hermite interpolation 3D unit test" << std::endl;
  for (int i = 0; i < nx.size(); ++i) std::cout << "\tl1_phi[" << i << "]: " << l1_phi[i] << std::endl;
  for (int i = 0; i < nx.size(); ++i) std::cout << "\tl1_psix[" << i << "]: " << l1_psix[i] << std::endl;
  for (int i = 0; i < nx.size(); ++i) std::cout << "\tl1_psiy[" << i << "]: " << l1_psiy[i] << std::endl;
  for (int i = 0; i < nx.size(); ++i) std::cout << "\tl1_psiz[" << i << "]: " << l1_psiz[i] << std::endl;
}
