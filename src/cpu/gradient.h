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

#include "./gradient/second-order-central.h"
#include "./gradient/third-order.h"
#include "array.h"
#include "grid.h"
#include "mat3.h"
#include "vec3.h"

namespace GALS
{
namespace CPU
{
/*! \class Gradient
 *
 * Class to create 3 component elements at a computational cell. For e.x. velocity, gradients, etc.
 */
template <typename T, typename T_GRID, typename GRADIENT_SCHEME = ThirdOrder<T, T_GRID>>
class Gradient
{
 public:
  typedef T value_type;

  /*! Default constructor
   */
  Gradient();

  /*! Destructor
   */
  ~Gradient();

  /*! Compute gradient of scalar.
   *
   * For a scalar the gradient is stored in Vec3, for e.x.:
   *  \f$\left(\nabla \alpha\right)_i = \{\partial_x \alpha, \partial_y \alpha, \partial_z \alpha\}_i.\f$
   *
   * \param alpha scalar variable for which gradient will be computed.
   * \param grad_alpha gradient of alpha is written to this.
   */
  static void compute(const Array<T_GRID, T> &alpha, Array<T_GRID, Vec3<T>> &grad_alpha);

  /*! Compute gradient of vector.
   *
   * For a Vec3 the gradient is stored in Mat3, for e.x.:
   *  Let \f$\boldsymbol{u} = \{u, v, w\} \f$ then,
   *  \f$\left(\nabla \boldsymbol{u}\right)_i = \begin{bmatrix}\partial_x u & \partial_y u & \partial_z u \\ \partial_x
   * v & \partial_y v & \partial_z v \\ \partial_x w & \partial_y w & \partial_z w \end{bmatrix}_i.\f$
   *
   * \param alpha vec3 variable for which gradient will be computed.
   * \param grad_alpha gradient of alpha is written to this.
   */
  static void compute(const Array<T_GRID, Vec3<T>> &alpha, Array<T_GRID, Mat3<T>> &grad_alpha);
};

}  // namespace CPU
}  // namespace GALS
