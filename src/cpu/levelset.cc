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

#include "gals/cpu/levelset.h"

#include <iostream>

template <typename T_GRID, typename T>
GALS::CPU::Levelset<T_GRID, T>::Levelset(const T_GRID& grid)
    : m_grid(grid),
      m_phi(grid),
      m_psi(grid),
      m_phi_mixed_derivatives(grid),
      m_phi_tm1(grid),
      m_psi_tm1(grid),
      m_phi_mixed_derivatives_tm1(grid),
      m_phi_interp_tm1(grid),
      m_psi_interp_tm1(grid)
{
}

template <typename T_GRID, typename T>
GALS::CPU::Levelset<T_GRID, T>::~Levelset()
{
}

template <typename T_GRID, typename T>
void GALS::CPU::Levelset<T_GRID, T>::print()
{
  std::cout << "inside levelset print" << std::endl;
}

template class GALS::CPU::Levelset<GALS::CPU::Grid<double, 1>, double>;
template class GALS::CPU::Levelset<GALS::CPU::Grid<double, 2>, double>;
// template class GALS::CPU::Levelset<GALS::CPU::Grid<double, 3>, double>;
