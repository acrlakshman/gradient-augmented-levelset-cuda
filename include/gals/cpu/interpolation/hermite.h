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

#include "gals/cpu/levelset.h"
#include "gals/utilities/array.h"
#include "gals/utilities/grid.h"
#include "gals/utilities/mat3.h"
#include "gals/utilities/utilities.h"
#include "gals/utilities/vec3.h"

namespace GALS
{
namespace INTERPOLATION
{
/*! \struct ControlPoints
 *
 * Control points for Bernstein form.
 */
template <typename T>
struct ControlPoints {
  T c_30, c_21, c_12, c_03;
};

// TODO (lakshman): Cleanup below functions.
//! cubic Bernstein polynomials
template <typename T = double>
T B0(T eta)
{
  return cube((T)1.0 - eta);
}

template <typename T = double>
T B1(T eta)
{
  return (T)3.0 * eta * sqr((T)1.0 - eta);
}

template <typename T = double>
T B2(T eta)
{
  return (T)3.0 * sqr(eta) * ((T)1.0 - eta);
}

template <typename T = double>
T B3(T eta)
{
  return cube(eta);
}

template <typename T = double>
T B0_Prime(T eta)
{
  return -(T)3.0 * sqr((T)1.0 - eta);
}

template <typename T = double>
T B1_Prime(T eta)
{
  return (T)3 * ((T)1 - eta) * ((T)1 - (T)3 * eta);
}

template <typename T = double>
T B2_Prime(T eta)
{
  return (T)3 * eta * ((T)2 - (T)3 * eta);
}

template <typename T = double>
T B3_Prime(T eta)
{
  return (T)3 * sqr(eta);
}

/*! Return control points as a struct.
 *
 * \param phi_0
 * \param phi_x_0
 * \param phi_1
 * \param phi_x_1
 * \param dx
 * \param use_gradient_limiting Implementation of gradient limiting algorithm by Bockmann et al. JCP 258(2014)
 *
 * \return ControlPoints<T>
 */
template <typename T>
static struct ControlPoints<T> get_control_points(const T &phi_0, const T &phi_x_0, const T &phi_1, const T &phi_x_1,
                                                  const T &dx, bool use_gradient_limiting = false) {
  struct ControlPoints<T> control_points;
  control_points.c_30 = phi_0;
  control_points.c_21 = phi_0 + (phi_x_0 * dx * (T)one_third);
  control_points.c_12 = phi_1 - (phi_x_1 * dx * (T)one_third);
  control_points.c_03 = phi_1;

  if (use_gradient_limiting) {
    // Implementation of gradient limiting algorithm by Bockmann et al. JCP 258(2014)
    T D1 = (T)3.0 * (control_points.c_21 - control_points.c_30);
    T D2 = (T)3.0 * (control_points.c_03 - control_points.c_12);
    T D3 = (T)1.5 * (control_points.c_12 - control_points.c_30);
    T D4 = (T)1.5 * (control_points.c_03 - control_points.c_21);
    T D5 = control_points.c_03 - control_points.c_30;

    if ((D1 - D3) * (D1 - D5) < (T)0)
      control_points.c_12 = (T)2 * control_points.c_21 - control_points.c_30;
    else if ((D2 - D4) * (D2 - D5) < (T)0)
      control_points.c_21 = (T)2 * control_points.c_12 - control_points.c_03;
    else if ((D1 - D5) * (D2 - D5) <= (T)0) {
      control_points.c_21 = (two_thirds * control_points.c_30) + (one_third * control_points.c_03);
      control_points.c_12 = (two_thirds * control_points.c_03) + (one_third * control_points.c_30);
    }
  }

  return control_points;
}

/*! \class Hermite
 *
 * Piece-wise hermite interpolation.
 */
template <typename T, typename T_GRID>
class Hermite;

/*! \class Hermite
 *
 * Piece-wise hermite interpolation for 1D.
 */
// Template specialized for 1D
template <typename T>
class Hermite<T, GALS::CPU::Grid<T, 1>>
{
 public:
  using value_type = T;
  using T_GRID = GALS::CPU::Grid<T, 1>;

  /*! Default constructor
   */
  Hermite(){};

  /*! Destructor
   */
  ~Hermite(){};

  /*! Piece-wise hermite interpolation.
   *
   * \param grid reference to grid.
   * \param x_interp position of a point where interpolation need to be performed.
   * \param levelset levelset field to interpolate.
   * \param use_gradient_limiting true/false
   */
  CPU::InterpolatedFields<CPU::Vec3<T>> interpolate(
      const GALS::CPU::Grid<typename GALS::CPU::Grid<T, 1>::value_type, GALS::CPU::Grid<T, 1>::dim> &grid,
      const typename GALS::CPU::Grid<T, 1>::position_type &x_interp,
      const CPU::Levelset<GALS::CPU::Grid<T, 1>, T> &levelset, const bool use_gradient_limiting = false);

  /*! Interpolate scalar field.
   *
   * \param x_interp interpolation points.
   * \param levelset levelset field to interpolate.
   */
  void compute(const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, typename GALS::CPU::Grid<T, 1>::position_type> &x_interp,
               CPU::Levelset<GALS::CPU::Grid<T, 1>, T> &levelset);

  /*! Overload operator to compute hermite interpolation of a scalar field.
   *
   * \param x_interp interpolation points.
   * \param levelset levelset field to interpolate.
   */
  void operator()(
      const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, typename GALS::CPU::Grid<T, 1>::position_type> &x_interp,
      CPU::Levelset<GALS::CPU::Grid<T, 1>, T> &levelset)
  {
    compute(x_interp, levelset);
  }

  /*! Overload operator to compute hermite interpolation of a scalar field.
   *
   * This is a placeholder function to avoid compilation error.
   *
   * \param x_interp interpolation points.
   * \param alpha variable to interpolate.
   * \param alpha_interpolated interpolated values are written to this variable.
   */
  void operator()(
      const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, typename GALS::CPU::Grid<T, 1>::position_type> &x_interp,
      const GALS::CPU::Array<GALS::CPU::Grid<T, 1>, T> &alpha,
      GALS::CPU::Array<GALS::CPU::Grid<T, 1>, T> &alpha_interpolated)
  {
    GALS_FUNCTION_NOT_IMPLEMENTED();
  }
};

/*! \class Hermite
 *
 * Piece-wise hermite interpolation for 2D.
 */
// Template specialized for 2D
template <typename T>
class Hermite<T, GALS::CPU::Grid<T, 2>>
{
 public:
  using value_type = T;
  using T_GRID = GALS::CPU::Grid<T, 2>;

  /*! Default constructor
   */
  Hermite(){};

  /*! Destructor
   */
  ~Hermite(){};

  /*! Piece-wise hermite interpolation.
   *
   * \param grid reference to grid.
   * \param x_interp position of a point where interpolation need to be performed.
   * \param levelset levelset field to interpolate.
   * \param use_gradient_limiting true/false
   */
  CPU::InterpolatedFields<CPU::Vec3<T>> interpolate(
      const GALS::CPU::Grid<typename GALS::CPU::Grid<T, 2>::value_type, GALS::CPU::Grid<T, 2>::dim> &grid,
      const typename GALS::CPU::Grid<T, 2>::position_type &x_interp,
      const CPU::Levelset<GALS::CPU::Grid<T, 2>, T> &levelset, const bool use_gradient_limiting = false);

  /*! Interpolate scalar field.
   *
   * \param x_interp interpolation points.
   * \param levelset levelset field to interpolate.
   */
  void compute(const GALS::CPU::Array<GALS::CPU::Grid<T, 2>, typename GALS::CPU::Grid<T, 2>::position_type> &x_interp,
               CPU::Levelset<GALS::CPU::Grid<T, 2>, T> &levelset);

  /*! Overload operator to compute hermite interpolation of a scalar field.
   *
   * \param x_interp interpolation points.
   * \param levelset levelset field to interpolate.
   */
  void operator()(
      const GALS::CPU::Array<GALS::CPU::Grid<T, 2>, typename GALS::CPU::Grid<T, 2>::position_type> &x_interp,
      CPU::Levelset<GALS::CPU::Grid<T, 2>, T> &levelset)
  {
    compute(x_interp, levelset);
  }

  /*! Overload operator to compute hermite interpolation of a scalar field.
   *
   * This is a placeholder function to avoid compilation error.
   *
   * \param x_interp interpolation points.
   * \param alpha variable to interpolate.
   * \param alpha_interpolated interpolated values are written to this variable.
   */
  void operator()(
      const GALS::CPU::Array<GALS::CPU::Grid<T, 2>, typename GALS::CPU::Grid<T, 2>::position_type> &x_interp,
      const GALS::CPU::Array<GALS::CPU::Grid<T, 2>, T> &alpha,
      GALS::CPU::Array<GALS::CPU::Grid<T, 2>, T> &alpha_interpolated)
  {
    GALS_FUNCTION_NOT_IMPLEMENTED();
  }
};

/*! \class Hermite
 *
 * Piece-wise hermite interpolation for 3D.
 */
// Template specialized for 3D
template <typename T>
class Hermite<T, GALS::CPU::Grid<T, 3>>
{
 public:
  using value_type = T;
  using T_GRID = GALS::CPU::Grid<T, 3>;

  /*! Default constructor
   */
  Hermite(){};

  /*! Destructor
   */
  ~Hermite(){};

  /*! Piece-wise hermite interpolation.
   *
   * \param grid reference to grid.
   * \param x_interp position of a point where interpolation need to be performed.
   * \param levelset levelset field to interpolate.
   * \param use_gradient_limiting true/false
   */
  CPU::InterpolatedFields<CPU::Vec3<T>> interpolate(
      const GALS::CPU::Grid<typename GALS::CPU::Grid<T, 3>::value_type, GALS::CPU::Grid<T, 3>::dim> &grid,
      const typename GALS::CPU::Grid<T, 3>::position_type &x_interp,
      const CPU::Levelset<GALS::CPU::Grid<T, 3>, T> &levelset, const bool use_gradient_limiting = false);

  /*! Interpolate scalar field.
   *
   * \param x_interp interpolation points.
   * \param levelset levelset field to interpolate.
   */
  void compute(const GALS::CPU::Array<GALS::CPU::Grid<T, 3>, typename GALS::CPU::Grid<T, 3>::position_type> &x_interp,
               CPU::Levelset<GALS::CPU::Grid<T, 3>, T> &levelset);

  /*! Overload operator to compute hermite interpolation of a scalar field.
   *
   * \param x_interp interpolation points.
   * \param levelset levelset field to interpolate.
   */
  void operator()(
      const GALS::CPU::Array<GALS::CPU::Grid<T, 3>, typename GALS::CPU::Grid<T, 3>::position_type> &x_interp,
      CPU::Levelset<GALS::CPU::Grid<T, 3>, T> &levelset)
  {
    compute(x_interp, levelset);
  }

  /*! Overload operator to compute hermite interpolation of a scalar field.
   *
   * This is a placeholder function to avoid compilation error.
   *
   * \param x_interp interpolation points.
   * \param alpha variable to interpolate.
   * \param alpha_interpolated interpolated values are written to this variable.
   */
  void operator()(
      const GALS::CPU::Array<GALS::CPU::Grid<T, 3>, typename GALS::CPU::Grid<T, 3>::position_type> &x_interp,
      const GALS::CPU::Array<GALS::CPU::Grid<T, 3>, T> &alpha,
      GALS::CPU::Array<GALS::CPU::Grid<T, 3>, T> &alpha_interpolated)
  {
    GALS_FUNCTION_NOT_IMPLEMENTED();
  }
};

}  // namespace INTERPOLATION
}  // namespace GALS
