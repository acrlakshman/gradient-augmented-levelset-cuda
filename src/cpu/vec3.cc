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

#include "gals/utilities/vec3.h"

#include "gals/utilities/utilities.h"

template <typename T>
GALS::CPU::Vec3<T>::Vec3(const T a, const T b, const T c)
{
  m_data[0] = a;
  m_data[1] = b;
  m_data[2] = c;
}

template <typename T>
GALS::CPU::Vec3<T>::Vec3() : Vec3(static_cast<T>(0), static_cast<T>(0), static_cast<T>(0))
{
}

template <typename T>
GALS::CPU::Vec3<T>::~Vec3()
{
}

template <typename T>
const int GALS::CPU::Vec3<T>::size() const
{
  return SIZE;
}

template <typename T>
const T GALS::CPU::Vec3<T>::operator[](const int idx) const
{
  return m_data[idx];
}

template <typename T>
T& GALS::CPU::Vec3<T>::operator[](const int idx)
{
  return m_data[idx];
}

template <typename T>
void GALS::CPU::Vec3<T>::operator=(const Vec3<T>& vec)
{
  for (int i = 0; i < SIZE; ++i) m_data[i] = vec[i];
}

template <typename T>
bool GALS::CPU::Vec3<T>::operator==(const Vec3<T>& vec) const
{
  return (GALS::is_equal(m_data[0], vec[0]) && GALS::is_equal(m_data[1], vec[1]) && GALS::is_equal(m_data[2], vec[2]));
}

template <typename T>
const GALS::CPU::Vec3<T> GALS::CPU::Vec3<T>::operator-(const Vec3<T>& vec) const
{
  return Vec3<T>(m_data[0] - vec[0], m_data[1] - vec[1], m_data[2] - vec[2]);
}

template <typename T>
const GALS::CPU::Vec3<T> GALS::CPU::Vec3<T>::operator*(const Vec3<T>& vec) const
{
  return Vec3<T>(m_data[0] * vec[0], m_data[1] * vec[1], m_data[2] * vec[2]);
}

template <typename T>
const GALS::CPU::Vec3<T> GALS::CPU::Vec3<T>::operator*(const T var) const
{
  return Vec3<T>(m_data[0] * var, m_data[1] * var, m_data[2] * var);
}

template <typename T>
const GALS::CPU::Vec3<T> GALS::CPU::Vec3<T>::operator/(const Vec3<T>& vec) const
{
  return Vec3<T>(m_data[0] / vec[0], m_data[1] / vec[1], m_data[2] / vec[2]);
}

template class GALS::CPU::Vec3<int>;
template class GALS::CPU::Vec3<double>;
