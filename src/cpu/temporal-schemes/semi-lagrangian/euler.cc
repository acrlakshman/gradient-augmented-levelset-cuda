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

#include "gals/cpu/temporal-schemes/semi-lagrangian/euler.h"
#include "gals/cpu/gradient.h"
#include "gals/cpu/interpolate.h"

#include <math.h>

template <typename T, typename T_GRID>
GALS::TEMPORAL_SCHEMES::SEMI_LAGRANGIAN::Euler<T, T_GRID>::Euler()
{
}

template <typename T, typename T_GRID>
GALS::TEMPORAL_SCHEMES::SEMI_LAGRANGIAN::Euler<T, T_GRID>::~Euler()
{
}

template <typename T, typename T_GRID>
void GALS::TEMPORAL_SCHEMES::SEMI_LAGRANGIAN::Euler<T, T_GRID>::compute(
    const T dt, const GALS::CPU::LevelsetVelocity<T_GRID, T> &levelset_velocity,
    GALS::CPU::Levelset<T_GRID, T> &levelset)
{
  const auto &grid = levelset.grid();
  const auto &velocity = levelset_velocity.velocity();
  const auto &velocity_grad = levelset_velocity.velocityGradient();
  const GALS::CPU::Vec3<int> num_cells = grid.numCells();
  auto &phi = levelset.phi();
  auto &psi = levelset.psi();
  // std::cout << "phi = " << phi << std::endl;

  // Compute x_root `root of the characteristic`.
  GALS::CPU::Array<T_GRID, GALS::CPU::Vec3<T>> x_root(grid);

  for (int i = 0; i < num_cells[0]; ++i)
    for (int j = 0; j < num_cells[1]; ++j)
      for (int k = 0; k < num_cells[2]; ++k) {
        x_root(i, j, k) = grid(i, j, k) - velocity(i, j, k) * dt;
        // std::cout << "grid(" << i << "," << j << "," << k << "): " << grid(i, j, k) << std::endl;
        // std::cout << "\tvelocity(" << i << "," << j << "," << k << "): " << velocity(i,j,k) << std::endl;
        // std::cout << "\tdt = " << dt << std::endl;
        // std::cout << "\tx_root(" << i << "," << j << "," << k << "): " << x_root(i,j,k) << std::endl;
      }
  // std::cout << "x_root = " << x_root << std::endl;

  // Compute phi_interp_prev, psi_interp_prev at x_root.
  GALS::CPU::Interpolate<T, T_GRID, GALS::INTERPOLATION::Hermite<T, T_GRID>>::compute(x_root, levelset);
  // std::cout << "after interp" << std::endl;

  // Compute x_root_grad.
  GALS::CPU::Array<T_GRID, GALS::CPU::Mat3<T>> x_root_grad(grid);

  const auto &identity_mat = T_GRID::identity_mat;

  for (int i = 0; i < num_cells[0]; ++i)
    for (int j = 0; j < num_cells[1]; ++j)
      for (int k = 0; k < num_cells[2]; ++k) x_root_grad(i, j, k) = identity_mat - velocity_grad(i, j, k) * dt;
  // std::cout << "x_root_grad = " << x_root_grad << std::endl;

  const auto &phi_interp_prev = levelset.phiInterpPrev();
  const auto &psi_interp_prev = levelset.psiInterpPrev();
  // std::cout << "phi_interp_prev = " << phi_interp_prev << std::endl;
  // std::cout << "psi_interp_prev = " << psi_interp_prev << std::endl;

  // Update phi and psi.
  for (int i = 0; i < num_cells[0]; ++i)
    for (int j = 0; j < num_cells[1]; ++j)
      for (int k = 0; k < num_cells[2]; ++k) {
        phi(i, j, k) = phi_interp_prev(i, j, k);
        psi(i, j, k) = x_root_grad(i, j, k).dot(psi_interp_prev(i, j, k));
      }
  // std::cout << "phi = " << phi << std::endl;
  // std::cout << "psi = " << psi << std::endl;
}

template class GALS::TEMPORAL_SCHEMES::SEMI_LAGRANGIAN::Euler<double, GALS::CPU::Grid<double, 1>>;
template class GALS::TEMPORAL_SCHEMES::SEMI_LAGRANGIAN::Euler<double, GALS::CPU::Grid<double, 2>>;
template class GALS::TEMPORAL_SCHEMES::SEMI_LAGRANGIAN::Euler<double, GALS::CPU::Grid<double, 3>>;
