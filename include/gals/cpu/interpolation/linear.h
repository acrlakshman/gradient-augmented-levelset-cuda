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
#include "gals/utilities/mat3.h"
#include "gals/utilities/vec3.h"

namespace GALS
{
namespace INTERPOLATION
{
/*! \class Linear
 *
 * Piece-wise linear interpolation.
 */
template <typename T, typename T_GRID>
class Linear;

/*! \class Linear
 *
 * Piece-wise linear interpolation for 1D.
 */
// Template specialized for 1D
template <typename T>
class Linear<T, GALS::CPU::Grid<T, 1>>
{
 public:
  using value_type = T;
  using T_GRID = GALS::CPU::Grid<T, 1>;

  /*! Default constructor
   */
  Linear(){};

  /*! Destructor
   */
  ~Linear(){};

  /*! Piece-wise linear interpolation.
   *
   * \param grid reference to grid.
   * \param x_interp position of a point where interpolation need to be performed.
   * \param alpha variable to interpolate.
   */
  T linearInterpolation(const GALS::CPU::Grid<typename GALS::CPU::Grid<T, 1>::value_type, T_GRID::dim> &grid,
                        const typename GALS::CPU::Grid<T, 1>::position_type &x_interp,
                        const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, T> &alpha);

  /*! Interpolate scalar field.
   *
   * \param x_interp interpolation points.
   * \param alpha variable to interpolate.
   * \param alpha_interpolated interpolated values are written to this variable.
   */
  void compute(const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, typename GALS::CPU::Grid<T, 1>::position_type> &x_interp,
               const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, T> &alpha,
               GALS::CPU::Array<GALS::CPU::Grid<T, 1>, T> &alpha_interpolated);

  /*! Overload operator to compute linear interpolation of a scalar field.
   *
   * \param x_interp interpolation points.
   * \param alpha variable to interpolate.
   * \param alpha_interpolated interpolated values are written to this variable.
   */
  void operator()(
      const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, typename GALS::CPU::Grid<T, 1>::position_type> &x_interp,
      const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, T> &alpha,
      GALS::CPU::Array<GALS::CPU::Grid<T, 1>, T> &alpha_interpolated)
  {
    compute(x_interp, alpha, alpha_interpolated);
  }

  /*! Interpolate scalar field.
   *
   * \param x_interp interpolation points.
   * \param levelset levelset field to interpolate.
   */
  void compute(const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, typename GALS::CPU::Grid<T, 1>::position_type> &x_interp,
               CPU::Levelset<GALS::CPU::Grid<T, 1>, T> &levelset);

  /*! Overload operator to compute linear interpolation of a scalar field.
   *
   * \param x_interp interpolation points.
   * \param levelset levelset field to interpolate.
   */
  void operator()(
      const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, typename GALS::CPU::Grid<T, 1>::position_type> &x_interp,
      CPU::Levelset<GALS::CPU::Grid<T, 1>, T> &levelset)
  {
    compute(x_interp, levelset);
  }
};

}  // namespace INTERPOLATION
}  // namespace GALS
