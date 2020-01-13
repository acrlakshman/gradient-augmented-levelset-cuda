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

#include "gals/cpu/interpolation/hermite.h"

#include <math.h>

template <typename T>
GALS::CPU::InterpolatedFields<GALS::CPU::Vec3<T>> GALS::INTERPOLATION::Hermite<T, GALS::CPU::Grid<T, 2>>::interpolate(
    const GALS::CPU::Grid<typename GALS::CPU::Grid<T, 2>::value_type, GALS::CPU::Grid<T, 2>::dim> &grid,
    const GALS::CPU::Vec3<int> &node_id, const typename GALS::CPU::Grid<T, 2>::position_type &x_interp,
    const GALS::CPU::Levelset<GALS::CPU::Grid<T, 2>, T> &levelset, const bool use_gradient_limiting)
{
  GALS::CPU::InterpolatedFields<GALS::CPU::Vec3<T>> hermite_fields;

  using T_GRID = GALS::CPU::Grid<T, 2>;
  const int dim = T_GRID::dim;
  const int axis_x = 0, axis_y = 1;
  const auto &axis_vectors = GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim>::axis_vectors;
  const GALS::CPU::Vec3<int> base_i_j = grid.baseNodeId(x_interp);
  const GALS::CPU::Vec3<int> base_ip1_j = GALS::CPU::Vec3<int>(base_i_j[0] + axis_vectors(0, 0), base_i_j[1], 0);
  const GALS::CPU::Vec3<int> base_i_jp1 = GALS::CPU::Vec3<int>(base_i_j[0], base_i_j[1] + axis_vectors(1, 1), 0);
  const GALS::CPU::Vec3<int> base_ip1_jp1 =
      GALS::CPU::Vec3<int>(base_i_j[0] + axis_vectors(0, 0), base_i_j[1] + axis_vectors(1, 1), 0);
  // std::cout << "x_interp = " << x_interp << std::endl;
  // std::cout << "\t\t\tbase i j = " << base_i_j << std::endl;
  // std::cout << "here" << std::endl;

  hermite_fields.phi_interpolated = levelset.phiPrev()(node_id);
  hermite_fields.psi_interpolated = levelset.psiPrev()(node_id);

  if (!grid.isIndexInDomain(base_i_j) || !grid.isIndexInDomain(base_ip1_j) || !grid.isIndexInDomain(base_i_jp1) ||
      !grid.isIndexInDomain(base_ip1_jp1))
    return hermite_fields;

  const typename T_GRID::position_type x_base = grid(base_i_j);
  const auto &dx = grid.dX();
  const auto &one_over_dx = grid.oneOverDX();

  GALS::CPU::Vec3<T> eta = (x_interp - x_base) * one_over_dx;
  // std::cout << "here" << std::endl;

  const ControlPoints<T> &control_points_bottom = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_i_j), levelset.psiPrev()(base_i_j)[axis_x], levelset.phiPrev()(base_ip1_j),
      levelset.psiPrev()(base_ip1_j)[axis_x], dx[axis_x], use_gradient_limiting);
  const ControlPoints<T> &control_points_top = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_i_jp1), levelset.psiPrev()(base_i_jp1)[axis_x], levelset.phiPrev()(base_ip1_jp1),
      levelset.psiPrev()(base_ip1_jp1)[axis_x], dx[axis_x], use_gradient_limiting);
  const ControlPoints<T> &control_points_left = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_i_j), levelset.psiPrev()(base_i_j)[axis_y], levelset.phiPrev()(base_i_jp1),
      levelset.psiPrev()(base_i_jp1)[axis_y], dx[axis_y], use_gradient_limiting);
  const ControlPoints<T> &control_points_right = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_ip1_j), levelset.psiPrev()(base_ip1_j)[axis_y], levelset.phiPrev()(base_ip1_jp1),
      levelset.psiPrev()(base_ip1_jp1)[axis_y], dx[axis_y], use_gradient_limiting);
  // std::cout << "here" << std::endl;

  T c_1 =
      levelset.psiPrev()(base_i_j)[axis_y] + dx[axis_x] * one_third * levelset.phiMixedDerivativesPrev()(base_i_j)[0];
  // std::cout << "c_1" << std::endl;
  T c_2 = levelset.psiPrev()(base_ip1_j)[axis_x] +
          dx[axis_y] * one_third * levelset.phiMixedDerivativesPrev()(base_ip1_j)[0];
  // std::cout << "c_2" << std::endl;
  T c_2_prime = levelset.psiPrev()(base_ip1_jp1)[axis_x] -
                dx[axis_y] * one_third * levelset.phiMixedDerivativesPrev()(base_ip1_jp1)[0];
  // std::cout << "c_2p" << std::endl;
  T c_3 = levelset.psiPrev()(base_i_jp1)[axis_y] +
          dx[axis_x] * one_third * levelset.phiMixedDerivativesPrev()(base_i_jp1)[0];
  // std::cout << "c_3" << std::endl;
  // std::cout << "here" << std::endl;

  T bx[] = {B0(eta[0]), B1(eta[0]), B2(eta[0]), B3(eta[0])};
  T bx_prime[] = {B0_Prime(eta[0]), B1_Prime(eta[0]), B2_Prime(eta[0]), B3_Prime(eta[0])};
  T by[] = {B0(eta[1]), B1(eta[1]), B2(eta[1]), B3(eta[1])};
  T by_prime[] = {B0_Prime(eta[1]), B1_Prime(eta[1]), B2_Prime(eta[1]), B3_Prime(eta[1])};
  // std::cout << "here" << std::endl;

  T control_points_all[] = {control_points_left.c_30,
                            control_points_left.c_21,
                            control_points_left.c_12,
                            control_points_left.c_03,
                            control_points_bottom.c_21,
                            control_points_bottom.c_21 + dx[axis_y] * one_third * c_1,
                            control_points_top.c_21 - dx[axis_y] * one_third * c_3,
                            control_points_top.c_21,
                            control_points_bottom.c_12,
                            control_points_right.c_21 - dx[axis_x] * one_third * c_2,
                            control_points_right.c_12 - dx[axis_x] * one_third * c_2_prime,
                            control_points_top.c_12,
                            control_points_right.c_30,
                            control_points_right.c_21,
                            control_points_right.c_12,
                            control_points_right.c_03};
  // std::cout << "here" << std::endl;

  for (int i = 0, k = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      hermite_fields.phi_interpolated += bx[i] * by[j] * control_points_all[k];
      ++k;
    }
  // std::cout << "computed phi interp" << std::endl;

  for (int i = 0, k = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      hermite_fields.psi_interpolated[axis_x] += bx_prime[i] * by[j] * control_points_all[k];
      ++k;
    }
  hermite_fields.psi_interpolated[axis_x] *= one_over_dx[axis_x];
  // std::cout << "computed psi interp x" << std::endl;

  for (int i = 0, k = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      hermite_fields.psi_interpolated[axis_y] += bx[i] * by_prime[j] * control_points_all[k];
      ++k;
    }
  hermite_fields.psi_interpolated[axis_y] *= one_over_dx[axis_y];
  // std::cout << "computed psi interp y" << std::endl;

  return hermite_fields;
}

template <typename T>
void GALS::INTERPOLATION::Hermite<T, GALS::CPU::Grid<T, 2>>::compute(
    const GALS::CPU::Array<GALS::CPU::Grid<T, 2>, typename GALS::CPU::Grid<T, 2>::position_type> &x_interp,
    GALS::CPU::Levelset<GALS::CPU::Grid<T, 2>, T> &levelset)
{
  typedef GALS::CPU::Grid<T, 2> T_GRID;

  const GALS::CPU::Vec3<int> num_cells_interp = x_interp.numCells();
  const T_GRID &grid = levelset.grid();
  const GALS::CPU::Vec3<typename T_GRID::value_type> dx = grid.dX();
  const auto &axis_vectors = GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim>::axis_vectors;

  std::cout << "Num cells: " << num_cells_interp << std::endl;
  for (int i = 0; i < num_cells_interp[0]; ++i)
    for (int j = 0; j < num_cells_interp[1]; ++j)
      for (int k = 0; k < num_cells_interp[2]; ++k) {
        GALS::CPU::Vec3<int> node_id(i, j, k);
        const auto &hermite_fields = this->interpolate(grid, node_id, x_interp(i, j, k), levelset);

        levelset.phiInterpPrev()(i, j, k) = hermite_fields.phi_interpolated;
        levelset.psiInterpPrev()(i, j, k) = hermite_fields.psi_interpolated;

        // std::cout << "\tgrid.X(" << i << "," << j << "," << k << "): " << grid(i,j,k) << std::endl;
        // std::cout << "\t\tlevelset.phi(" << i << "," << j << "," << k << "): " << levelset.phiPrev()(i,j,k) <<
        // std::endl; std::cout << "\t\tlevelset.psi(" << i << "," << j << "," << k << "): " << levelset.psiPrev()(i,j,k)
        // << std::endl; std::cout << "\t\tx_interp(" << i << "," << j << "," << k << "): " << x_interp(i,j,k) <<
        // std::endl; std::cout << "\t\tphi_interp(" << i << "," << j << "," << k << "): " <<
        // levelset.phiInterpPrev()(i,j,k) << std::endl; std::cout << "\t\tpsi_interp(" << i << "," << j << "," << k <<
        // "): " << levelset.psiInterpPrev()(i,j,k) << std::endl;
      }
}

template class GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 2>>;
