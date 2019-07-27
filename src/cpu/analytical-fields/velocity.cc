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

#include "gals/analytical-fields/velocity.h"
#include "gals/utilities/utilities.h"

template <typename T_GRID, typename T>
GALS::ANALYTICAL_FIELDS::Velocity<T_GRID, T>::Velocity(const T_GRID& grid, const GALS::INPUT_FIELDS::Velocity& inputs)
    : m_grid(grid), m_inputs(inputs)
{
}

template <typename T_GRID, typename T>
GALS::ANALYTICAL_FIELDS::Velocity<T_GRID, T>::~Velocity()
{
}

template <typename T_GRID, typename T>
void GALS::ANALYTICAL_FIELDS::Velocity<T_GRID, T>::compute(
    const GALS::CPU::Array<T_GRID, GALS::CPU::Vec3<T>>& positions,
    GALS::CPU::LevelsetVelocity<T_GRID, T>& levelset_velocity)
{
  const GALS::CPU::Vec3<int> num_cells = positions.numCells();
  const auto& velocity_name = m_inputs.name;

  switch (velocity_name_map[velocity_name]) {
    case VelocityFieldNames::CIRCULAR: {
      auto& velocity = levelset_velocity.velocity();

      if (T_GRID::dim == 1)
        GALS_FUNCTION_NOT_IMPLEMENTED(
            "GALS::ANALYTICAL_FIELDS::Velocity::compute: CIRCULAR velocity field for 1D grid.");

      T coeff = pi() / static_cast<T>(3.15);
      GALS::CPU::Vec3<T> velocity_node;

      for (int i = 0; i < num_cells[0]; ++i)
        for (int j = 0; j < num_cells[1]; ++j)
          for (int k = 0; k < num_cells[2]; ++k) {
            for (int cmpt = 0; cmpt < T_GRID::dim; ++cmpt) {
              if (cmpt == 0) velocity_node[cmpt] = coeff * (m_inputs.center[cmpt] - positions(i, j, k)[cmpt + 1]);
              if (cmpt == 1) velocity_node[cmpt] = coeff * (positions(i, j, k)[cmpt - 1] - m_inputs.center[cmpt]);
            }

            velocity(i, j, k) = velocity_node;
          }

      // Compute velocity gradient using analytical expressions.
      if (!std::strcmp(m_inputs.gradient_scheme, "ANALYTICAL")) {
        auto& velocity_gradient = levelset_velocity.velocityGradient();

        for (int i = 0; i < num_cells[0]; ++i)
          for (int j = 0; j < num_cells[1]; ++j)
            for (int k = 0; k < num_cells[2]; ++k) {
              for (int cmpt = 0; cmpt < T_GRID::dim; ++cmp) {
                for (int derv = 0; derv < T_GRID::dim; ++derv) {
                  if (cmpt == 0 && derv == 0) velocity_gradient[cmpt][derv] = static_cast<T>(0);
                  if (cmpt == 0 && derv == 1) velocity_gradient[cmpt][derv] = -coeff;
                  if (cmpt == 1 && derv == 0) velocity_gradient[cmpt][derv] = coeff;
                  if (cmpt == 1 && derv == 1) velocity_gradient[cmpt][derv] = static_cast<T>(0);
                }
              }
            }
      }

      break;
    }

    default:
      break;
  }
}

template class GALS::ANALYTICAL_FIELDS::Velocity<GALS::CPU::Grid<double, 1>, double>;
template class GALS::ANALYTICAL_FIELDS::Velocity<GALS::CPU::Grid<double, 2>, double>;
template class GALS::ANALYTICAL_FIELDS::Velocity<GALS::CPU::Grid<double, 3>, double>;
