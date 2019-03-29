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

#include <vector>

namespace GALS
{
namespace CPU
{
/*! \class Vec3
 *
 * Class to create 3 component elements at a computational cell. For e.x. velocity, gradients, etc.
 */
template <typename T>
class Vec3
{
 public:
  typedef T value_type;

  /*! Default constructor
   */
  Vec3();

  /*! Destructor
   */
  ~Vec3();

  /*! Returns number of elements.
   *
   * \return number of elements.
   */
  const int size() const;

  /*! Overloaded subscript operator that returns a reference.
   *
   * \param idx zero based index of element.
   *
   * \return element at index (idx).
   */
  const T operator[](const int idx) const;

  /*! Overloaded subscript operator that returns a reference.
   *
   * \param idx zero based index of element.
   *
   * \return element at index (idx).
   */
  T &operator[](const int idx);

  /*! Overload assignment operator.
   *
   * \param vec variable whose values will be assigned.
   */
  void operator=(const Vec3<T> &vec);

  /*! Equality operator.
   *
   * \param vec variable to compare against.
   *
   * \return true if equal, false otherwise.
   */
  bool operator==(const Vec3<T> &vec);

 private:
  T m_data[3];
};

}  // namespace CPU
}  // namespace GALS
