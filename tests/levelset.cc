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
#include "gals/utilities/array.h"
#include "gals/utilities/utilities.h"
#include "gals/utilities/vec3.h"
#include "gals/utilities/vec_n.h"

#include <gtest/gtest.h>

#include <math.h>
#include <iostream>

/* * * * * *  TEST #0  * * * * * */
TEST(CPU, INTERPOLATED_FIELDS) { GALS::CPU::InterpolatedFields<GALS::CPU::Vec3<int>> interpolated_fields; }

/* * * * * *  TEST #1  * * * * * */
TEST(CPU, LEVELSET_1D)
{
  typedef GALS::CPU::Grid<double, 1> T_GRID;

  // initializing 1-D test grid.
  GALS::CPU::Grid<double, 1> grid(10, 1, 1);

  // grid generation
  grid.generate(-1, 1, -1, 1, -1, 1);

  // accessing grid details
  const auto mask = grid.getMask();
  const int pad = grid.getPadding();
  const auto num_cells = grid.numCells();

  // initializing 1-D scalar array
  GALS::CPU::Levelset<T_GRID, double> levelset(grid);
  levelset.print();

  // Test getters.
  const auto& grid_levelset = levelset.grid();
  const auto& phi = levelset.phi();
  const auto& psi = levelset.psi();
  const auto& phi_mixed_derivatives = levelset.phiMixedDerivatives();
  const auto& phi_tm1 = levelset.phiTm1();
  const auto& psi_tm1 = levelset.psiTm1();
  const auto& phi_mixed_derivatives_tm1 = levelset.phiMixedDerivativesTm1();
  const auto& phi_interp_tm1 = levelset.phiInterpTm1();
  const auto& psi_interp_tm1 = levelset.psiInterpTm1();
  auto& phi_1 = levelset.phi();
  auto& psi_1 = levelset.psi();
  auto& phi_mixed_derivatives_1 = levelset.phiMixedDerivatives();
  auto& phi_tm1_1 = levelset.phiTm1();
  auto& psi_tm1_1 = levelset.psiTm1();
  auto& phi_mixed_derivatives_tm1_1 = levelset.phiMixedDerivativesTm1();
  auto& phi_interp_tm1_1 = levelset.phiInterpTm1();
  auto& psi_interp_tm1_1 = levelset.psiInterpTm1();
}

/* * * * * *  TEST #2  * * * * * */
TEST(CPU, LEVELSET_2D)
{
  typedef GALS::CPU::Grid<double, 2> T_GRID;

  // initializing 1-D test grid.
  GALS::CPU::Grid<double, 2> grid(10, 10, 1);

  // grid generation
  grid.generate(-1, 1, -1, 1, -1, 1);

  // accessing grid details
  const auto mask = grid.getMask();
  const int pad = grid.getPadding();
  const auto num_cells = grid.numCells();

  // initializing 1-D scalar array
  GALS::CPU::Levelset<T_GRID, double> levelset(grid);
  levelset.print();

  // Test getters.
  const auto& grid_levelset = levelset.grid();
  const auto& phi = levelset.phi();
  const auto& psi = levelset.psi();
  const auto& phi_mixed_derivatives = levelset.phiMixedDerivatives();
  const auto& phi_tm1 = levelset.phiTm1();
  const auto& psi_tm1 = levelset.psiTm1();
  const auto& phi_mixed_derivatives_tm1 = levelset.phiMixedDerivativesTm1();
  const auto& phi_interp_tm1 = levelset.phiInterpTm1();
  const auto& psi_interp_tm1 = levelset.psiInterpTm1();
  auto& phi_1 = levelset.phi();
  auto& psi_1 = levelset.psi();
  auto& phi_mixed_derivatives_1 = levelset.phiMixedDerivatives();
  auto& phi_tm1_1 = levelset.phiTm1();
  auto& psi_tm1_1 = levelset.psiTm1();
  auto& phi_mixed_derivatives_tm1_1 = levelset.phiMixedDerivativesTm1();
  auto& phi_interp_tm1_1 = levelset.phiInterpTm1();
  auto& psi_interp_tm1_1 = levelset.psiInterpTm1();
}
