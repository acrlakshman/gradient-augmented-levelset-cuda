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

#include "gals/cpu/interpolation/linear.h"

// Template specialized for 2D
template <typename T>
T GALS::INTERPOLATION::Linear<T, GALS::CPU::Grid<T, 2>>::linearInterpolation(
    const GALS::CPU::Grid<typename GALS::CPU::Grid<T, 2>::value_type, GALS::CPU::Grid<T, 2>::dim> &grid,
    const typename GALS::CPU::Grid<T, 2>::position_type &x_interp,
    const GALS::CPU::Array<GALS::CPU::Grid<T, 2>, T> &alpha)
{
  std::cout << "specialized for 2D" << std::endl;
  T alpha_interpolated = 0;

  return alpha_interpolated;
}

template <typename T>
void GALS::INTERPOLATION::Linear<T, GALS::CPU::Grid<T, 2>>::compute(
    const GALS::CPU::Array<GALS::CPU::Grid<T, 2>, typename GALS::CPU::Grid<T, 2>::position_type> &x_interp,
    const GALS::CPU::Array<GALS::CPU::Grid<T, 2>, T> &alpha,
    GALS::CPU::Array<GALS::CPU::Grid<T, 2>, T> &alpha_interpolated)
{
  std::cout << "compute: specialized for 2D" << std::endl;
}

template class GALS::INTERPOLATION::Linear<double, GALS::CPU::Grid<double, 2>>;
