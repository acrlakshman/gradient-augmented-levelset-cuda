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

#include "gals/cpu/interpolation/hermite.h"
#include "gals/cpu/interpolation/linear.h"
#include "gals/cpu/levelset.h"
#include "gals/utilities/array.h"
#include "gals/utilities/grid.h"
#include "gals/utilities/mat3.h"
#include "gals/utilities/vec3.h"

namespace GALS
{
namespace CPU
{
/*! \class Interpolate
 *
 * Class to perform interpolation. Default interpolation scheme is set to GALS::INTERPOLATION::Linear<...>.
 */
template <typename T, typename T_GRID, typename INTERPOLATION_SCHEME = GALS::INTERPOLATION::Linear<T, T_GRID>>
class Interpolate
{
 public:
  using value_type = T;

  /*! Default constructor
   */
  Interpolate();

  /*! Destructor
   */
  ~Interpolate();

  /*! Compute interpolation of a scalar.
   *
   * \param x_interp interpolation points.
   * \param alpha variable to interpolate.
   * \param alpha_interpolated interpolated values are written to this variable.
   */
  static void compute(const Array<T_GRID, typename T_GRID::position_type> &x_interp, const Array<T_GRID, T> &alpha,
                      Array<T_GRID, T> &alpha_interpolated);

  /*! Compute interpolation of a levelset field.
   *
   * Depending on the interpolation scheme, phi and/or psi will be computed and the corresponding member variables of
   * Levelset class are updated. Variables of `_prev` are used to updated `_interp_prev` variables.
   *
   * \param x_interp interpolation points.
   * \param levelset variable whose members will be used and udpated during interpolation.
   */
  static void compute(const Array<T_GRID, typename T_GRID::position_type> &x_interp,
                      GALS::CPU::Levelset<T_GRID, double> &levelset);
};

}  // namespace CPU
}  // namespace GALS
