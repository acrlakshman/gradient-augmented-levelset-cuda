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

#include <iostream>
#include <vector>

namespace GALS
{
namespace CPU
{
/*! \class VecN
 *
 * Class to create varying size elements at a computational cell. For e.x. gradients, etc.
 */
template <typename T, int SIZE = 3>
class VecN
{
 public:
  using value_type = T;

  /*! Default constructor
   */
  VecN();

  /*! Destructor
   */
  ~VecN();

  /*! Returns number of elements.
   *
   * \return number of elements.
   */
  const int size() const;

  /*! Overloaded subscript operator.
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

  /*! Output operator overload.
   *
   * \param out output stream.
   * \param vec VecN object to output stream.
   *
   * \return reference to output stream.
   */
  friend std::ostream &operator<<(std::ostream &out, const VecN<T, SIZE> &vec)
  {
    for (int i = 0; i < vec.size(); ++i) out << vec[i] << "\t";

    return out;
  }

 private:
  std::vector<T> m_data;
};

}  // namespace CPU
}  // namespace GALS
