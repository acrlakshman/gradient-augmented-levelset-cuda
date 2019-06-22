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

#include "gals/utilities/vec_n.h"

template <typename T, int SIZE>
GALS::CPU::VecN<T, SIZE>::VecN()
{
  m_data.resize(SIZE);

  for (int i = 0; i < SIZE; ++i) m_data[i] = static_cast<T>(0);
}

template <typename T, int SIZE>
GALS::CPU::VecN<T, SIZE>::~VecN()
{
  m_data.clear();
  m_data.shrink_to_fit();
}

template <typename T, int SIZE>
const int GALS::CPU::VecN<T, SIZE>::size() const
{
  return m_data.size();
}

template <typename T, int SIZE>
const T GALS::CPU::VecN<T, SIZE>::operator[](const int idx) const
{
  return m_data[idx];
}

template <typename T, int SIZE>
T& GALS::CPU::VecN<T, SIZE>::operator[](const int idx)
{
  return m_data[idx];
}

template class GALS::CPU::VecN<int, 1>;
template class GALS::CPU::VecN<int, 2>;
template class GALS::CPU::VecN<int, 3>;
template class GALS::CPU::VecN<double, 1>;
template class GALS::CPU::VecN<double, 2>;
template class GALS::CPU::VecN<double, 3>;
