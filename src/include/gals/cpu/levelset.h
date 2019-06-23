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
#include "gals/utilities/vec_n.h"

namespace GALS
{
namespace CPU
{
/*! \struct InterpolatedFields
 *
 * Struct to store field interpolated cell quantities.
 */
template <class T_VECTOR>
struct InterpolatedFields {
  typedef typename T_VECTOR::value_type T;
  T h_phi_location;
  T_VECTOR h_phi_gradient_location;

  InterpolatedFields()
  {
    h_phi_location = T();
    h_phi_gradient_location = T_VECTOR();
  };

  ~InterpolatedFields(){};
};

/*! \class Levelset
 *
 * Class to create Levelset.
 */
template <typename T_GRID, typename T = double>
class Levelset
{
 public:
  /*! Constructor with grid.
   *
   * \param grid grid.
   */
  Levelset(const T_GRID& grid);

  /*! Remove default constructor.
   */
  Levelset() = delete;

  /*! Destructor
   */
  ~Levelset();

  /*! Return grid.
   *
   * \return grid.
   */
  const T_GRID& grid() const { return m_grid; }

  /*! Return phi.
   *
   * \return phi.
   */
  Array<T_GRID, T>& phi() { return m_phi; }

  /*! Return psi.
   *
   * \return psi.
   */
  Array<T_GRID, Vec3<T>>& psi() { return m_psi; }

  /*! Return phi_mixed_derivatives.
   *
   * \return phi_mixed_derivatives.
   */
  Array<T_GRID, VecN<T, T_GRID::num_mixed_derivatives>>& phiMixedDerivatives() { return m_phi_mixed_derivatives; }

  /*! Return phi_tm1.
   *
   * \return phi_tm1.
   */
  Array<T_GRID, T>& phiTm1() { return m_phi_tm1; }

  /*! Return psi_tm1.
   *
   * \return psi_tm1.
   */
  Array<T_GRID, Vec3<T>>& psiTm1() { return m_psi_tm1; }

  /*! Return phi_interp_tm1.
   *
   * \return phi_interp_tm1.
   */
  Array<T_GRID, T>& phiInterpTm1() { return m_phi_interp_tm1; }

  /*! Return psi_interp_tm1.
   *
   * \return psi_interp_tm1.
   */
  Array<T_GRID, Vec3<T>>& psiInterpTm1() { return m_psi_interp_tm1; }

  //! Print levelset values.
  void print();

 private:
  const T_GRID& m_grid;
  Array<T_GRID, T> m_phi;                                                         //! levelset field.
  Array<T_GRID, Vec3<T>> m_psi;                                                   //! gradient of levelset field.
  Array<T_GRID, VecN<T, T_GRID::num_mixed_derivatives>> m_phi_mixed_derivatives;  //! Mixed derivatives of phi.
  Array<T_GRID, T> m_phi_tm1;               //! levelset field at previous time step.
  Array<T_GRID, Vec3<T>> m_psi_tm1;         //! gradient of levelset field at previous time step.
  Array<T_GRID, T> m_phi_interp_tm1;        //! interpolated levelset field at previous time step.
  Array<T_GRID, Vec3<T>> m_psi_interp_tm1;  //! interpolated gradient of levelset field at previous time step.
};

}  // namespace CPU
}  // namespace GALS
