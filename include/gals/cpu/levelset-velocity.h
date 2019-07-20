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

#include <string>
#include <vector>

#include "gals/utilities/array.h"
#include "gals/utilities/grid.h"

namespace GALS
{
namespace CPU
{
/*! \class LevelsetVelocity
 *
 * Class to create LevelsetVelocity.
 */
template <typename T_GRID, typename T = double>
class LevelsetVelocity
{
 public:
  /*! Constructor with grid.
   *
   * \param grid grid.
   */
  LevelsetVelocity(const T_GRID& grid);

  /*! Remove default constructor.
   */
  LevelsetVelocity() = delete;

  /*! Destructor
   */
  ~LevelsetVelocity();

  /*! Return grid.
   *
   * \return grid.
   */
  const T_GRID& grid() const { return m_grid; }

  /*! Return velocity field.
   *
   * \return velocity field.
   */
  Array<T_GRID, Vec3<T>>& velocity() { return m_velocity; }

  /*! Return velocity field.
   *
   * \return velocity field.
   */
  const Array<T_GRID, Vec3<T>>& velocity() const { return m_velocity; }

  /*! Return gradient of velocity field.
   *
   * \return gradient of velocity field.
   */
  Array<T_GRID, Mat3<T>>& velocityGradient() { return m_velocity_grad; }

  /*! Return gradient of velocity field.
   *
   * \return gradient of velocity field.
   */
  const Array<T_GRID, Mat3<T>>& velocityGradient() const { return m_velocity_grad; }

 private:
  const T_GRID& m_grid;
  Array<T_GRID, Vec3<T>> m_velocity;       //! Velocity field.
  Array<T_GRID, Mat3<T>> m_velocity_grad;  //! Gradient of velocity field.
};

}  // namespace CPU
}  // namespace GALS
