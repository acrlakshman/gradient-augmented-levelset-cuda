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
T GALS::INTERPOLATION::Hermite<T, GALS::CPU::Grid<T, 1>>::interpolate(
    const GALS::CPU::Grid<typename GALS::CPU::Grid<T, 1>::value_type, GALS::CPU::Grid<T, 1>::dim> &grid,
    const typename GALS::CPU::Grid<T, 1>::position_type &x_interp,
    const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, T> &alpha)
{
  typedef GALS::CPU::Grid<T, 1> T_GRID;

  const int dim = T_GRID::dim;
  const auto &axis_vectors = GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim>::axis_vectors;
  const GALS::CPU::Vec3<int> base_node_id = grid.baseNodeId(x_interp);
  const GALS::CPU::Vec3<int> base_node_id_p1 = GALS::CPU::Vec3<int>(base_node_id[0] + axis_vectors(0, 0), 0, 0);

  // TODO (lakshman): Update
  const typename T_GRID::position_type x_base = grid(base_node_id);
  const auto &one_over_dx = grid.oneOverDX();

  /////
  // const int axis = 0;
  // T eta = (x_interp[axis] - x_base[axis]) * one_over_dx[axis];
  // const struct ControlPoints<T> &control_points = get_control_points(
  // alpha(base_node_id), phi_gradient_ghost(reference_index)(axis),
  // phi_ghost(reference_index + T_INDEX::Axis_Vector(axis)),
  // phi_gradient_ghost(reference_index + T_INDEX::Axis_Vector(axis))(axis), axis, use_gradient_limiting);

  // hermite_fields.h_phi_location = control_points.c_30 * B0(eta) + control_points.c_21 * B1(eta) +
  // control_points.c_12 * B2(eta) + control_points.c_03 * B3(eta);
  // hermite_fields.h_phi_gradient_location(axis) =
  //(control_points.c_30 * B0_Prime(eta) + control_points.c_21 * B1_Prime(eta) + control_points.c_12 * B2_Prime(eta) +
  // control_points.c_03 * B3_Prime(eta)) *
  // grid.One_Over_DX()(axis);
  //...
  T eta = (x_interp[0] - x_base[0]) * one_over_dx[0];

  T alpha_interpolated = (1. - eta) * alpha(base_node_id) + eta * alpha(base_node_id_p1);

  // TODO (lakshman): DELETE
  // std::cout << "base_node_id = " << base_node_id << "; eta = " << eta
  //<< "; alpha_interpolation = " << alpha_interpolated << std::endl;
  // T xo = 0., ro = 0.5;
  // T alpha_exact = (x_interp[0] - xo) * (x_interp[0] - xo) - ro * ro;
  // std::cout << "error = " << fabs(alpha_interpolated - alpha_exact) << std::endl;

  return alpha_interpolated;
}

template <typename T>
void GALS::INTERPOLATION::Hermite<T, GALS::CPU::Grid<T, 1>>::compute(
    const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, typename GALS::CPU::Grid<T, 1>::position_type> &x_interp,
    const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, T> &alpha,
    GALS::CPU::Array<GALS::CPU::Grid<T, 1>, T> &alpha_interpolated)
{
  typedef GALS::CPU::Grid<T, 1> T_GRID;

  const GALS::CPU::Vec3<int> num_cells = alpha.numCells();
  const T_GRID &grid = alpha.grid();
  const GALS::CPU::Vec3<typename T_GRID::value_type> dx = grid.dX();
  const auto &axis_vectors = GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim>::axis_vectors;

  for (int i = 0; i < x_interp.grid().numCells()[0]; ++i)
    for (int j = 0; j < x_interp.grid().numCells()[1]; ++j)
      for (int k = 0; k < x_interp.grid().numCells()[2]; ++k)
        alpha_interpolated(i, j, k) = interpolate(grid, x_interp(i, j, k), alpha);
}

template class GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 1>>;
