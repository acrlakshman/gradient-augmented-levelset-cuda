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

#include "gals/cpu/gradient.h"

template <typename T, typename T_GRID, typename GRADIENT_SCHEME>
GALS::CPU::Gradient<T, T_GRID, GRADIENT_SCHEME>::Gradient()
{
}

template <typename T, typename T_GRID, typename GRADIENT_SCHEME>
GALS::CPU::Gradient<T, T_GRID, GRADIENT_SCHEME>::~Gradient()
{
}

template <typename T, typename T_GRID, typename GRADIENT_SCHEME>
void GALS::CPU::Gradient<T, T_GRID, GRADIENT_SCHEME>::compute(const Array<T_GRID, T> &alpha,
                                                              Array<T_GRID, Vec3<T>> &grad_alpha)
{
  GRADIENT_SCHEME()(alpha, grad_alpha);
}

template <typename T, typename T_GRID, typename GRADIENT_SCHEME>
void GALS::CPU::Gradient<T, T_GRID, GRADIENT_SCHEME>::compute(const Array<T_GRID, Vec3<T>> &alpha,
                                                              Array<T_GRID, Mat3<T>> &grad_alpha)
{
  GRADIENT_SCHEME()(alpha, grad_alpha);
}

template class GALS::CPU::Gradient<double, GALS::CPU::Grid<double, 1>,
                                   GALS::CPU::SecondOrderCentral<double, GALS::CPU::Grid<double, 1>>>;

template class GALS::CPU::Gradient<double, GALS::CPU::Grid<double, 2>,
                                   GALS::CPU::SecondOrderCentral<double, GALS::CPU::Grid<double, 2>>>;

template class GALS::CPU::Gradient<double, GALS::CPU::Grid<double, 3>,
                                   GALS::CPU::SecondOrderCentral<double, GALS::CPU::Grid<double, 3>>>;

template class GALS::CPU::Gradient<double, GALS::CPU::Grid<double, 1>,
                                   GALS::CPU::ThirdOrder<double, GALS::CPU::Grid<double, 1>>>;

template class GALS::CPU::Gradient<double, GALS::CPU::Grid<double, 2>,
                                   GALS::CPU::ThirdOrder<double, GALS::CPU::Grid<double, 2>>>;

template class GALS::CPU::Gradient<double, GALS::CPU::Grid<double, 3>,
                                   GALS::CPU::ThirdOrder<double, GALS::CPU::Grid<double, 3>>>;
