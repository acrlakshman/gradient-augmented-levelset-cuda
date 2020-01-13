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
GALS::CPU::InterpolatedFields<GALS::CPU::Vec3<T>> GALS::INTERPOLATION::Hermite<T, GALS::CPU::Grid<T, 1>>::interpolate(
    const GALS::CPU::Grid<typename GALS::CPU::Grid<T, 1>::value_type, GALS::CPU::Grid<T, 1>::dim> &grid,
    const GALS::CPU::Vec3<int> &node_id, const typename GALS::CPU::Grid<T, 1>::position_type &x_interp,
    const GALS::CPU::Levelset<GALS::CPU::Grid<T, 1>, T> &levelset, const bool use_gradient_limiting)
{
  GALS::CPU::InterpolatedFields<GALS::CPU::Vec3<T>> hermite_fields;

  using T_GRID = GALS::CPU::Grid<T, 1>;
  const int dim = T_GRID::dim;
  const int axis = 0;
  const auto &axis_vectors = GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim>::axis_vectors;
  const GALS::CPU::Vec3<int> base_node_id = grid.baseNodeId(x_interp);
  const GALS::CPU::Vec3<int> base_node_id_p1 = GALS::CPU::Vec3<int>(base_node_id[0] + axis_vectors(0, 0), 0, 0);

  const typename T_GRID::position_type x_base = grid(base_node_id);
  const auto &dx = grid.dX();
  const auto &one_over_dx = grid.oneOverDX();

  T eta = (x_interp[0] - x_base[0]) * one_over_dx[0];

  const ControlPoints<T> &control_points = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_node_id), levelset.psiPrev()(base_node_id)[axis], levelset.phiPrev()(base_node_id_p1),
      levelset.psiPrev()(base_node_id_p1)[axis], dx[axis], use_gradient_limiting);

  hermite_fields.phi_interpolated = control_points.c_30 * B0(eta) + control_points.c_21 * B1(eta) +
                                    control_points.c_12 * B2(eta) + control_points.c_03 * B3(eta);
  hermite_fields.psi_interpolated[axis] = (control_points.c_30 * B0_Prime(eta) + control_points.c_21 * B1_Prime(eta) +
                                           control_points.c_12 * B2_Prime(eta) + control_points.c_03 * B3_Prime(eta)) *
                                          one_over_dx[axis];

  // Debug
  // std::cout << std::scientific;
  // std::cout << "eta = " << eta << "; xi = " << x_interp[0] << std::endl
  //<< "\t; phi_b = " << levelset.phiPrev()(base_node_id)
  //<< "\t; phi_bp1 = " << levelset.phiPrev()(base_node_id_p1) << std::endl
  //<< "\t; psi_b = " << levelset.psiPrev()(base_node_id)
  //<< "\t; psi_bp1 = " << levelset.psiPrev()(base_node_id_p1) << std::endl
  //<< "\t; phi_i = " << hermite_fields.phi_interpolated
  //<< "\t; psi_i = " << hermite_fields.psi_interpolated[axis] << std::endl;
  // std::cout << "\tfirst: c30 = " << control_points.c_30 << "; B0(eta) = " << B0(eta) << std::endl;
  // std::cout << "\tsecond: c21 = " << control_points.c_21 << "; B1(eta) = " << B1(eta) << std::endl;
  // std::cout << "\tthird: c12 = " << control_points.c_12 << "; B2(eta) = " << B2(eta) << std::endl;
  // std::cout << "\tfourth: c03 = " << control_points.c_03 << "; B3(eta) = " << B3(eta) << std::endl;

  return hermite_fields;
}

template <typename T>
void GALS::INTERPOLATION::Hermite<T, GALS::CPU::Grid<T, 1>>::compute(
    const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, typename GALS::CPU::Grid<T, 1>::position_type> &x_interp,
    GALS::CPU::Levelset<GALS::CPU::Grid<T, 1>, T> &levelset)
{
  typedef GALS::CPU::Grid<T, 1> T_GRID;

  const GALS::CPU::Vec3<int> num_cells_interp = x_interp.numCells();
  const T_GRID &grid = levelset.grid();
  const GALS::CPU::Vec3<typename T_GRID::value_type> dx = grid.dX();
  const auto &axis_vectors = GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim>::axis_vectors;

  for (int i = 0; i < num_cells_interp[0]; ++i)
    for (int j = 0; j < num_cells_interp[1]; ++j)
      for (int k = 0; k < num_cells_interp[2]; ++k) {
        GALS::CPU::Vec3<int> node_id(i, j, k);
        const auto &hermite_fields = this->interpolate(grid, node_id, x_interp(i, j, k), levelset);

        levelset.phiInterpPrev()(i, j, k) = hermite_fields.phi_interpolated;
        levelset.psiInterpPrev()(i, j, k) = hermite_fields.psi_interpolated;
      }
}

template class GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 1>>;
