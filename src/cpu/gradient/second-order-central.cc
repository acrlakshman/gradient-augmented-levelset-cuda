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

#include "gals/cpu/gradient/second-order-central.h"

template <typename T, typename T_GRID>
GALS::CPU::SecondOrderCentral<T, T_GRID>::SecondOrderCentral()
{
}

template <typename T, typename T_GRID>
GALS::CPU::SecondOrderCentral<T, T_GRID>::~SecondOrderCentral()
{
}

template <typename T, typename T_GRID>
void GALS::CPU::SecondOrderCentral<T, T_GRID>::compute(const Array<T_GRID, T> &alpha,
                                                       Array<T_GRID, Vec3<T>> &grad_alpha)
{
  const Vec3<int> num_cells = alpha.numCells();
  const T_GRID &grid = alpha.grid();
  const Vec3<typename T_GRID::value_type> dx = grid.dX();
  const auto &axis_vectors = GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim>::axis_vectors;

  for (int i = 0; i < num_cells[0]; ++i)
    for (int j = 0; j < num_cells[1]; ++j)
      for (int k = 0; k < num_cells[2]; ++k) {
        for (int axis = 0; axis < T_GRID::dim; ++axis) {
          typename T_GRID::value_type one_by_dx = static_cast<typename T_GRID::value_type>(1.) / dx[axis];

          grad_alpha(i, j, k)[axis] =
              (alpha(i + axis_vectors(axis, 0), j + axis_vectors(axis, 1), k + axis_vectors(axis, 2)) -
               alpha(i - axis_vectors(axis, 0), j - axis_vectors(axis, 1), k - axis_vectors(axis, 2))) *
              one_by_dx * static_cast<T>(0.5);
        }
      }
}

template <typename T, typename T_GRID>
void GALS::CPU::SecondOrderCentral<T, T_GRID>::compute(const Array<T_GRID, Vec3<T>> &alpha,
                                                       Array<T_GRID, Mat3<T>> &grad_alpha)
{
  const Vec3<int> num_cells = alpha.numCells();
  const T_GRID &grid = alpha.grid();
  const Vec3<typename T_GRID::value_type> dx = grid.dX();
  const auto &axis_vectors = GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim>::axis_vectors;

  for (int i = 0; i < num_cells[0]; ++i)
    for (int j = 0; j < num_cells[1]; ++j)
      for (int k = 0; k < num_cells[2]; ++k) {
        for (int axis = 0; axis < T_GRID::dim; ++axis) {
          for (int cmpt = 0; cmpt < T_GRID::dim; ++cmpt) {
            typename T_GRID::value_type one_by_dx = static_cast<typename T_GRID::value_type>(1.) / dx[cmpt];

            grad_alpha(i, j, k)(axis, cmpt) =
                (alpha(i + axis_vectors(axis, 0), j + axis_vectors(axis, 1), k + axis_vectors(axis, 2))[axis] -
                 alpha(i - axis_vectors(axis, 0), j - axis_vectors(axis, 1), k - axis_vectors(axis, 2))[axis]) *
                one_by_dx * static_cast<T>(0.5);
          }
        }
      }
}

template class GALS::CPU::SecondOrderCentral<double, GALS::CPU::Grid<double, 1>>;
template class GALS::CPU::SecondOrderCentral<double, GALS::CPU::Grid<double, 2>>;
template class GALS::CPU::SecondOrderCentral<double, GALS::CPU::Grid<double, 3>>;
