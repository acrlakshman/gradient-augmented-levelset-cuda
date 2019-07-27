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

#include "gals/analytical-fields/levelset.h"
#include "gals/utilities/utilities.h"

template <typename T_GRID, typename T>
GALS::ANALYTICAL_FIELDS::Levelset<T_GRID, T>::Levelset(const T_GRID& grid, const GALS::INPUT_FIELDS::Levelset& inputs)
    : m_grid(grid), m_inputs(inputs)
{
}

template <typename T_GRID, typename T>
GALS::ANALYTICAL_FIELDS::Levelset<T_GRID, T>::~Levelset()
{
}

template <typename T_GRID, typename T>
void GALS::ANALYTICAL_FIELDS::Levelset<T_GRID, T>::compute(
    const GALS::CPU::Array<T_GRID, GALS::CPU::Vec3<T>>& positions, GALS::CPU::Levelset<T_GRID, T>& levelset)
{
  const GALS::CPU::Vec3<int> num_cells = positions.numCells();
  const auto& levelset_name = m_inputs.name;

  switch (levelset_name_map[levelset_name]) {
    case LevelsetFieldNames::CIRCLE: {
      auto& phi = levelset.phi();

      GALS::CPU::Vec3<T> distance_vec;

      for (int i = 0; i < num_cells[0]; ++i)
        for (int j = 0; j < num_cells[1]; ++j)
          for (int k = 0; k < num_cells[2]; ++k) {
            distance_vec = positions(i, j, k) - m_inputs.center;
            phi(i, j, k) = distance_vec.mag() - m_inputs.radius;
          }

      break;
    }

    default:
      break;
  }
}

template class GALS::ANALYTICAL_FIELDS::Levelset<GALS::CPU::Grid<double, 1>, double>;
template class GALS::ANALYTICAL_FIELDS::Levelset<GALS::CPU::Grid<double, 2>, double>;
template class GALS::ANALYTICAL_FIELDS::Levelset<GALS::CPU::Grid<double, 3>, double>;
