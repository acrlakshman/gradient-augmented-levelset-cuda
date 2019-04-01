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

#include "linear.h"

#include <iostream>

#include "../vec3.h"

template <typename T, typename T_GRID>
GALS::INTERPOLATION::Linear<T, T_GRID>::Linear()
{
}

template <typename T, typename T_GRID>
GALS::INTERPOLATION::Linear<T, T_GRID>::~Linear()
{
}

template <typename T, typename T_GRID>
T GALS::INTERPOLATION::Linear<T, T_GRID>::linearInterpolation(
    const GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim> &grid,
    const typename T_GRID::position_type &x_interp, const GALS::CPU::Array<T_GRID, T> &alpha)
{
  T alpha_interpolated;

  const GALS::CPU::Vec3<int> base_node_id = grid.baseNodeId(x_interp);

  for (int d = 0; d < T_GRID::dim; ++d) {
    // TODO (lakshman), complete this.
  }

  return alpha_interpolated;
}

template <typename T, typename T_GRID>
void GALS::INTERPOLATION::Linear<T, T_GRID>::compute(
    const GALS::CPU::Array<T_GRID, typename T_GRID::position_type> &x_interp, const GALS::CPU::Array<T_GRID, T> &alpha,
    GALS::CPU::Array<T_GRID, T> &alpha_interpolated)
{
  // const Vec3<int> num_cells = alpha.numCells();
  // const T_GRID &grid = alpha.grid();
  // const Vec3<typename T_GRID::value_type> dx = grid.dX();
  // const auto &axis_vectors = GALS::INTERPOLATION::Grid<typename T_GRID::value_type, T_GRID::dim>::axis_vectors;

  // for (int i = 0; i < num_cells[0]; ++i)
  // for (int j = 0; j < num_cells[1]; ++j)
  // for (int k = 0; k < num_cells[2]; ++k) {
  // for (int axis = 0; axis < T_GRID::dim; ++axis) {
  // typename T_GRID::value_type one_by_dx = static_cast<typename T_GRID::value_type>(1.) / dx[axis];

  // grad_alpha(i, j, k)[axis] =
  //(alpha(i + axis_vectors(axis, 0), j + axis_vectors(axis, 1), k + axis_vectors(axis, 2)) -
  // alpha(i - axis_vectors(axis, 0), j - axis_vectors(axis, 1), k - axis_vectors(axis, 2))) *
  // one_by_dx * static_cast<T>(0.5);
  //}
  //}
}

template class GALS::INTERPOLATION::Linear<double, GALS::CPU::Grid<double, 1>>;
template class GALS::INTERPOLATION::Linear<double, GALS::CPU::Grid<double, 2>>;
template class GALS::INTERPOLATION::Linear<double, GALS::CPU::Grid<double, 3>>;
