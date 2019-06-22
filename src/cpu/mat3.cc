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

#include "gals/utilities/mat3.h"

#include "gals/utilities/utilities.h"

template <typename T>
GALS::CPU::Mat3<T>::Mat3(T a00, T a01, T a02, T a10, T a11, T a12, T a20, T a21, T a22)
{
  m_data[0] = a00;
  m_data[1] = a01;
  m_data[2] = a02;
  m_data[3] = a10;
  m_data[4] = a11;
  m_data[5] = a12;
  m_data[6] = a20;
  m_data[7] = a21;
  m_data[8] = a22;
}

template <typename T>
GALS::CPU::Mat3<T>::Mat3()
    : Mat3(static_cast<T>(0), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
           static_cast<T>(0), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0))
{
}

template <typename T>
GALS::CPU::Mat3<T>::~Mat3()
{
}

template <typename T>
const int GALS::CPU::Mat3<T>::size() const
{
  return SIZE;
}

template <typename T>
const T GALS::CPU::Mat3<T>::operator[](const int idx) const
{
  return m_data[idx];
}

template <typename T>
T& GALS::CPU::Mat3<T>::operator[](const int idx)
{
  return m_data[idx];
}

template <typename T>
const T GALS::CPU::Mat3<T>::operator()(const int i_row, const int i_col) const
{
  return m_data[i_row * 3 + i_col];
}

template <typename T>
T& GALS::CPU::Mat3<T>::operator()(const int i_row, const int i_col)
{
  return m_data[i_row * 3 + i_col];
}

template <typename T>
void GALS::CPU::Mat3<T>::operator=(const Mat3<T>& mat)
{
  for (int i = 0; i < SIZE; ++i) m_data[i] = mat[i];
}

template <typename T>
bool GALS::CPU::Mat3<T>::operator==(const Mat3<T>& mat) const
{
  return (GALS::is_equal(m_data[0], mat[0]) && GALS::is_equal(m_data[1], mat[1]) && GALS::is_equal(m_data[2], mat[2]) &&
          GALS::is_equal(m_data[3], mat[3]) && GALS::is_equal(m_data[4], mat[4]) && GALS::is_equal(m_data[5], mat[5]) &&
          GALS::is_equal(m_data[6], mat[6]) && GALS::is_equal(m_data[7], mat[7]) && GALS::is_equal(m_data[8], mat[8]));
}

template class GALS::CPU::Mat3<int>;
template class GALS::CPU::Mat3<double>;
