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

#include "gals/cpu/temporal.h"
#include "gals/utilities/array.h"
#include "gals/utilities/utilities.h"

#include <gtest/gtest.h>

#include <math.h>
#include <iostream>

TEST(CPU, TEMPORAL_SCHEME_EULER)
{
  using T = double;
  using T_GRID = GALS::CPU::Grid<T, 1>;
  // scalar array on 1D grid.
  T_GRID grid(10, 1, 1);
  grid.generate(-1, 1, -1, 1, -1, 1);

  const T dt = 1.;
  GALS::CPU::Array<T_GRID, T> alpha(grid);
  GALS::CPU::Array<T_GRID, T> convection(grid);
  GALS::CPU::Array<T_GRID, T> alpha_new(grid);
  GALS::CPU::Levelset<T_GRID, T> levelset(grid);

  levelset.phi() = alpha_new;
  levelset.phiTm1() = alpha;

  // TODO: Complete test case.
  GALS::CPU::Temporal<T, T_GRID, GALS::TEMPORAL_SCHEMES::Euler<T, T_GRID>>::compute(dt, convection, levelset);

  // For test converage.
  GALS::CPU::Temporal<T, T_GRID> temporal_scheme;
}
