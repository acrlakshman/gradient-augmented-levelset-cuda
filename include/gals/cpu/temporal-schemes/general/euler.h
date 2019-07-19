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

#pragma once

#include "gals/cpu/levelset.h"
#include "gals/utilities/array.h"
#include "gals/utilities/grid.h"

namespace GALS
{
namespace TEMPORAL_SCHEMES
{
namespace GENERAL
{
/*! \class Euler
 *
 * Euler scheme for temporal update.
 */
template <typename T, typename T_GRID>
class Euler
{
 public:
  using value_type = T;
  using grid_type = T_GRID;

  /*! Default constructor
   */
  Euler() {}

  /*! Destructor
   */
  ~Euler() {}

  /*! Advect in time.
   *
   * \f$\frac{\phi^{n+1} - \phi^n}{dt} + U\cdot\nabla\phi = 0.\f$
   * Ghost cells are not updated during this step.
   *
   * \param dt time step.
   * \param velocity velocity term.
   * \param levelset levelset function that needs to be advected.
   */
  void compute(const T dt, const CPU::Array<T_GRID, CPU::Vec3<T>> &velocity, GALS::CPU::Levelset<T_GRID, T> &levelset);

  /*! Advect in time.
   *
   * \f$\frac{\phi^{n+1} - \phi^n}{dt} + U\cdot\nable\phi = 0.\f$
   * Ghost cells are not updated during this step.
   *
   * \param dt time step.
   * \param velocity velocity term.
   * \param levelset levelset function that needs to be advected.
   */
  void operator()(const T dt, const CPU::Array<T_GRID, CPU::Vec3<T>> &velocity,
                  GALS::CPU::Levelset<T_GRID, T> &levelset)
  {
    this->compute(dt, velocity, levelset);
  }
};

}  // namespace GENERAL
}  // namespace TEMPORAL_SCHEMES
}  // namespace GALS
