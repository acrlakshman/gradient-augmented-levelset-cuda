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

#include "gals/cpu/temporal-schemes/euler.h"

#include <math.h>

template <typename T, typename T_GRID>
GALS::TEMPORAL_SCHEMES::Euler<T, T_GRID>::Euler()
{
}

template <typename T, typename T_GRID>
GALS::TEMPORAL_SCHEMES::Euler<T, T_GRID>::~Euler()
{
}

template <typename T, typename T_GRID>
void GALS::TEMPORAL_SCHEMES::Euler<T, T_GRID>::compute(const T dt,
                                                       const GALS::CPU::LevelsetVelocity<T_GRID, T> &levelset_velocity,
                                                       GALS::CPU::Levelset<T_GRID, T> &levelset)
{
  const auto &phi_prev = levelset.phiPrev();
  auto &phi = levelset.phi();
  const GALS::CPU::Vec3<int> num_cells = phi.numCells();
  const auto &velocity = levelset_velocity.velocity();

  for (int i = 0; i < num_cells[0]; ++i)
    for (int j = 0; j < num_cells[1]; ++j)
      for (int k = 0; k < num_cells[2]; ++k) {
        phi(i, j, k) = phi_prev(i, j, k) - dt * velocity(i, j, k)[0];  // FIXME
      }
}

template class GALS::TEMPORAL_SCHEMES::Euler<double, GALS::CPU::Grid<double, 1>>;
