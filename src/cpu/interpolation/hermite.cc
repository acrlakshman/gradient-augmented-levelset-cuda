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

#include "hermite.h"

#include <iostream>

#include "../vec3.h"

template <typename T, typename T_GRID>
GALS::INTERPOLATION::Hermite<T, T_GRID>::Hermite()
{
}

template <typename T, typename T_GRID>
GALS::INTERPOLATION::Hermite<T, T_GRID>::~Hermite()
{
}

template <typename T, typename T_GRID>
T GALS::INTERPOLATION::Hermite<T, T_GRID>::interpolate(
    const GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim> &grid,
    const typename T_GRID::position_type &x_interp, const GALS::CPU::Array<T_GRID, T> &alpha)
{
  std::cout << "generalized" << std::endl;
  T alpha_interpolated;

  const GALS::CPU::Vec3<int> base_node_id = grid.baseNodeId(x_interp);

  const typename T_GRID::position_type x_base = grid(base_node_id);
  const auto &one_over_dx = grid.oneOverDX();

  GALS::CPU::Vec3<T> eta;
  for (int d = 0; d < T_GRID::dim; ++d) eta[d] = (x_interp[d] - x_base[d]) * one_over_dx[d];

  for (int d = 0; d < T_GRID::dim; ++d) {
    // TODO (lakshman), complete this.
  }

  return alpha_interpolated;
}

template <typename T, typename T_GRID>
void GALS::INTERPOLATION::Hermite<T, T_GRID>::compute(
    const GALS::CPU::Array<T_GRID, typename T_GRID::position_type> &x_interp, const GALS::CPU::Array<T_GRID, T> &alpha,
    GALS::CPU::Array<T_GRID, T> &alpha_interpolated)
{
  std::cout << "compute: generalized" << std::endl;
}

template class GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 1>>;
template class GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 2>>;
template class GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 3>>;
