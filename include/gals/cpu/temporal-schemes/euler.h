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

#include "gals/utilities/array.h"
#include "gals/utilities/grid.h"

namespace GALS
{
namespace TEMPORAL_SCHEMES
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

  /*! Default constructor
   */
  Euler();

  /*! Destructor
   */
  ~Euler();

  /*! Advect in time.
   *
   * \frac{\alpha^{n+1} - \alpha^n}{dt} + convection = 0.
   * Ghost cells are not updated during this step.
   *
   * \param dt time step.
   * \param alpha variable to advect in time.
   * \param convection convection term.
   * \param alpha_new values after advection.
   */
  void compute(const T dt, const CPU::Array<T_GRID, T> &alpha, const CPU::Array<T_GRID, T> &convection,
               CPU::Array<T_GRID, T> &alpha_new);

  /*! Advect in time.
   *
   * \frac{\alpha^{n+1} - \alpha^n}{dt} + convection = 0.
   * Ghost cells are not updated during this step.
   *
   * \param dt time step.
   * \param alpha variable to advect in time.
   * \param convection convection term.
   * \param alpha_new values after advection.
   */
  void operator()(const T dt, const CPU::Array<T_GRID, T> &alpha, const CPU::Array<T_GRID, T> &convection,
                  CPU::Array<T_GRID, T> &alpha_new)
  {
    this->compute(dt, alpha, convection, alpha_new);
  }
};

}  // namespace TEMPORAL_SCHEMES
}  // namespace GALS
