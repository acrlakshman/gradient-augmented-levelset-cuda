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

namespace GALS
{
namespace CPU
{
/*! \class Mat3
 *
 * Class to create 9 component variable at a computational cell. For e.x. gradient of velocity.
 */
template <typename T>
class Mat3
{
 public:
  using value_type = T;

  static constexpr int SIZE = 9;

  /*! Constructor.
   *
   * \param a00 row 0, col 0
   * \param a01 row 0, col 1
   * \param a02 row 0, col 2
   * \param a10 row 1, col 0
   * \param a11 row 1, col 1
   * \param a12 row 1, col 2
   * \param a20 row 2, col 0
   * \param a21 row 2, col 1
   * \param a22 row 2, col 2
   */
  Mat3(T a00, T a01, T a02, T a10, T a11, T a12, T a20, T a21, T a22);

  /*! Default constructor.
   */
  Mat3();

  /*! Destructor
   */
  ~Mat3();

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

  /*! Row column based element access.
   *
   * \param i_row row index (zero-based).
   * \param i_col column index (zero-based).
   *
   * \return element at index (i_row, i_col).
   */
  const T operator()(const int i_row, const int i_col) const;

  /*! Row column based element access.
   *
   * \param i_row row index (zero-based).
   * \param i_col column index (zero-based).
   *
   * \return element at index (i_row, i_col).
   */
  T &operator()(const int i_row, const int i_col);

  /*! Overload assignment operator.
   *
   * \param mat variable whose values will be assigned.
   */
  void operator=(const Mat3<T> &mat);

  /*! Equality operator.
   *
   * \param mat variable to compare against.
   *
   * \return true if equal, false otherwise.
   */
  bool operator==(const Mat3<T> &mat) const;

  /*! Output operator overload.
   *
   * \param out output stream.
   * \param mat Mat3 object to output stream.
   *
   * \return reference to output stream.
   */
  friend std::ostream &operator<<(std::ostream &out, const Mat3<T> &mat)
  {
    for (int i = 0; i < mat.size(); ++i) out << mat[i] << "\t";

    return out;
  }

 private:
  T m_data[SIZE];
};

}  // namespace CPU
}  // namespace GALS
