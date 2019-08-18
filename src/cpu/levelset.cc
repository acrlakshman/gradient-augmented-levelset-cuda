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

#include "gals/cpu/levelset.h"
#include "gals/cpu/gradient.h"

#include <iostream>

template <typename T_GRID, typename T>
GALS::CPU::Levelset<T_GRID, T>::Levelset(const T_GRID& grid)
    : m_grid(grid),
      m_phi(grid),
      m_psi(grid),
      m_phi_mixed_derivatives(grid),
      m_phi_prev(grid),
      m_psi_prev(grid),
      m_phi_mixed_derivatives_prev(grid),
      m_phi_interp_prev(grid),
      m_psi_interp_prev(grid)
{
}

template <typename T_GRID, typename T>
GALS::CPU::Levelset<T_GRID, T>::~Levelset()
{
}

template <typename T_GRID, typename T>
void GALS::CPU::Levelset<T_GRID, T>::computeMixedDerivatives(
    const Array<T_GRID, Vec3<T>>& psi, Array<T_GRID, VecN<T, T_GRID::num_mixed_derivatives>>& phi_mixed_derivatives)
{
  // For now gradient computations are hardcoded in this function.
  const Vec3<int> num_cells = psi.numCells();
  const T_GRID& grid = psi.grid();
  const Vec3<typename T_GRID::value_type> dx = grid.dX();
  const auto one_over_dx = grid.oneOverDX();
  const auto& axis_vectors = GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim>::axis_vectors;

  if constexpr (T_GRID::dim == 2) {
    int axis_x = 0, axis_y = 1;

    // Computing \frac{\partial}{\partial y} \left( \frac{\partial}{\partial x} \right).
    for (int i = 0; i < num_cells[0]; ++i)
      for (int j = 0; j < num_cells[1]; ++j)
        for (int k = 0; k < num_cells[2]; ++k) {
          phi_mixed_derivatives(i, j, k)[0] =
              (psi(i + axis_vectors(axis_y, 0), j + axis_vectors(axis_y, 1), k + axis_vectors(axis_y, 2))[axis_x] -
               psi(i - axis_vectors(axis_y, 0), j - axis_vectors(axis_y, 1), k - axis_vectors(axis_y, 2))[axis_x]) *
              one_over_dx[axis_y] * static_cast<T>(0.5);
        }
  } else if constexpr (T_GRID::dim == 3) {
    int axis_x = 0, axis_y = 1, axis_z = 2;

    // Computing \frac{\partial}{\partial x} \left( \frac{\partial}{\partial y} \right).
    for (int i = 0; i < num_cells[0]; ++i)
      for (int j = 0; j < num_cells[1]; ++j)
        for (int k = 0; k < num_cells[2]; ++k) {
          phi_mixed_derivatives(i, j, k)[0] =
              (psi(i + axis_vectors(axis_x, 0), j + axis_vectors(axis_x, 1), k + axis_vectors(axis_x, 2))[axis_y] -
               psi(i - axis_vectors(axis_x, 0), j - axis_vectors(axis_x, 1), k - axis_vectors(axis_x, 2))[axis_y]) *
              one_over_dx[axis_x] * static_cast<T>(0.5);
        }

    // Computing \frac{\partial}{\partial y} \left( \frac{\partial}{\partial z} \right).
    for (int i = 0; i < num_cells[0]; ++i)
      for (int j = 0; j < num_cells[1]; ++j)
        for (int k = 0; k < num_cells[2]; ++k) {
          phi_mixed_derivatives(i, j, k)[1] =
              (psi(i + axis_vectors(axis_y, 0), j + axis_vectors(axis_y, 1), k + axis_vectors(axis_y, 2))[axis_z] -
               psi(i - axis_vectors(axis_y, 0), j - axis_vectors(axis_y, 1), k - axis_vectors(axis_y, 2))[axis_z]) *
              one_over_dx[axis_y] * static_cast<T>(0.5);
        }

    // Computing \frac{\partial}{\partial z} \left( \frac{\partial}{\partial x} \right).
    for (int i = 0; i < num_cells[0]; ++i)
      for (int j = 0; j < num_cells[1]; ++j)
        for (int k = 0; k < num_cells[2]; ++k) {
          phi_mixed_derivatives(i, j, k)[2] =
              (psi(i + axis_vectors(axis_z, 0), j + axis_vectors(axis_z, 1), k + axis_vectors(axis_z, 2))[axis_x] -
               psi(i - axis_vectors(axis_z, 0), j - axis_vectors(axis_z, 1), k - axis_vectors(axis_z, 2))[axis_x]) *
              one_over_dx[axis_z] * static_cast<T>(0.5);
        }

    // Computing \frac{\partial}{\partial y} \left( \frac{\partial}{\partial z} \right).
    for (int i = 0; i < num_cells[0]; ++i)
      for (int j = 0; j < num_cells[1]; ++j)
        for (int k = 0; k < num_cells[2]; ++k) {
          // phi_mixed_derivatives(i, j, k)[1] =
          //(psi(i + axis_vectors(axis_y, 0), j + axis_vectors(axis_y, 1), k + axis_vectors(axis_y, 2))[axis_z] -
          // psi(i - axis_vectors(axis_y, 0), j - axis_vectors(axis_y, 1), k - axis_vectors(axis_y, 2))[axis_z]) *
          // one_over_dx[axis_y] * static_cast<T>(0.5);
        }
  }
}

template <typename T_GRID, typename T>
void GALS::CPU::Levelset<T_GRID, T>::print()
{
  std::cout << "inside levelset print" << std::endl;
}

template class GALS::CPU::Levelset<GALS::CPU::Grid<double, 1>, double>;
template class GALS::CPU::Levelset<GALS::CPU::Grid<double, 2>, double>;
template class GALS::CPU::Levelset<GALS::CPU::Grid<double, 3>, double>;
