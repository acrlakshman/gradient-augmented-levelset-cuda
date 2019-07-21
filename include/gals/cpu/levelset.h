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
  using T = typename T_VECTOR::value_type;

  T phi_interpolated;
  T_VECTOR psi_interpolated;

  InterpolatedFields()
  {
    phi_interpolated = T();
    psi_interpolated = T_VECTOR();
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

  /*! Return phi.
   *
   * \return phi.
   */
  const Array<T_GRID, T>& phi() const { return m_phi; }

  /*! Return psi.
   *
   * \return psi.
   */
  Array<T_GRID, Vec3<T>>& psi() { return m_psi; }

  /*! Return psi.
   *
   * \return psi.
   */
  const Array<T_GRID, Vec3<T>>& psi() const { return m_psi; }

  /*! Return phi_mixed_derivatives.
   *
   * \return phi_mixed_derivatives.
   */
  Array<T_GRID, VecN<T, T_GRID::num_mixed_derivatives>>& phiMixedDerivatives() { return m_phi_mixed_derivatives; }

  /*! Return phi_mixed_derivatives.
   *
   * \return phi_mixed_derivatives.
   */
  const Array<T_GRID, VecN<T, T_GRID::num_mixed_derivatives>>& phiMixedDerivatives() const
  {
    return m_phi_mixed_derivatives;
  }

  /*! Return phi_prev.
   *
   * \return phi_prev.
   */
  Array<T_GRID, T>& phiPrev() { return m_phi_prev; }

  /*! Return phi_prev.
   *
   * \return phi_prev.
   */
  const Array<T_GRID, T>& phiPrev() const { return m_phi_prev; }

  /*! Return psi_prev.
   *
   * \return psi_prev.
   */
  Array<T_GRID, Vec3<T>>& psiPrev() { return m_psi_prev; }

  /*! Return psi_prev.
   *
   * \return psi_prev.
   */
  const Array<T_GRID, Vec3<T>>& psiPrev() const { return m_psi_prev; }

  /*! Return phi_mixed_derivatives_prev.
   *
   * \return phi_mixed_derivatives_prev.
   */
  Array<T_GRID, VecN<T, T_GRID::num_mixed_derivatives>>& phiMixedDerivativesPrev()
  {
    return m_phi_mixed_derivatives_prev;
  }

  /*! Return phi_mixed_derivatives_prev.
   *
   * \return phi_mixed_derivatives_prev.
   */
  const Array<T_GRID, VecN<T, T_GRID::num_mixed_derivatives>>& phiMixedDerivativesPrev() const
  {
    return m_phi_mixed_derivatives_prev;
  }

  /*! Return phi_interp_prev.
   *
   * \return phi_interp_prev.
   */
  Array<T_GRID, T>& phiInterpPrev() { return m_phi_interp_prev; }

  /*! Return phi_interp_prev.
   *
   * \return phi_interp_prev.
   */
  const Array<T_GRID, T>& phiInterpPrev() const { return m_phi_interp_prev; }

  /*! Return psi_interp_prev.
   *
   * \return psi_interp_prev.
   */
  Array<T_GRID, Vec3<T>>& psiInterpPrev() { return m_psi_interp_prev; }

  /*! Return psi_interp_prev.
   *
   * \return psi_interp_prev.
   */
  const Array<T_GRID, Vec3<T>>& psiInterpPrev() const { return m_psi_interp_prev; }

  //! Print levelset values.
  void print();

 private:
  const T_GRID& m_grid;
  Array<T_GRID, T> m_phi;                                                         //! levelset field.
  Array<T_GRID, Vec3<T>> m_psi;                                                   //! gradient of levelset field.
  Array<T_GRID, VecN<T, T_GRID::num_mixed_derivatives>> m_phi_mixed_derivatives;  //! Mixed derivatives of phi.
  Array<T_GRID, T> m_phi_prev;        //! levelset field at previous time step.
  Array<T_GRID, Vec3<T>> m_psi_prev;  //! gradient of levelset field at previous time step.
  Array<T_GRID, VecN<T, T_GRID::num_mixed_derivatives>>
      m_phi_mixed_derivatives_prev;          //! Mixed derivatives of phi at previous time step.
  Array<T_GRID, T> m_phi_interp_prev;        //! interpolated levelset field at previous time step.
  Array<T_GRID, Vec3<T>> m_psi_interp_prev;  //! interpolated gradient of levelset field at previous time step.
};

}  // namespace CPU
}  // namespace GALS
