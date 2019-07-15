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
#include "gals/cpu/temporal-schemes/euler.h"
#include "gals/utilities/array.h"
#include "gals/utilities/grid.h"

namespace GALS
{
namespace CPU
{
/*! \class Temporal
 *
 * Class to perform temporal update. Default temporal scheme is set to GALS::TEMPORAL_SCHEMES::Euler<...>.
 */
template <typename T, typename T_GRID, typename TEMPORAL_SCHEME = GALS::TEMPORAL_SCHEMES::Euler<T, T_GRID>>
class Temporal
{
 public:
  using value_type = T;

  /*! Default constructor
   */
  Temporal();

  /*! Destructor
   */
  ~Temporal();

  /*! Advect in time.
   *
   * NOTE: Ghost cells are not updated during this step.
   *
   * \param dt time step.
   * \param convection convection term.
   * \param levelset levelset function that needs to be advected.
   */
  static void compute(const T dt, const Array<T_GRID, T> &convection, GALS::CPU::Levelset<T_GRID, T> &levelset);
};

}  // namespace CPU
}  // namespace GALS
