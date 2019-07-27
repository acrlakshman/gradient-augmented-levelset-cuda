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

#include "gals/cpu/interpolation/hermite.h"

#include <math.h>

template <typename T>
GALS::CPU::InterpolatedFields<GALS::CPU::Vec3<T>> GALS::INTERPOLATION::Hermite<T, GALS::CPU::Grid<T, 3>>::interpolate(
    const GALS::CPU::Grid<typename GALS::CPU::Grid<T, 3>::value_type, GALS::CPU::Grid<T, 3>::dim>& grid,
    const typename GALS::CPU::Grid<T, 3>::position_type& x_interp,
    const GALS::CPU::Levelset<GALS::CPU::Grid<T, 3>, T>& levelset, const bool use_gradient_limiting)
{
  GALS::CPU::InterpolatedFields<GALS::CPU::Vec3<T>> hermite_fields;

  using T_GRID = GALS::CPU::Grid<T, 3>;
  const int dim = T_GRID::dim;
  const int axis_x = 0, axis_y = 1, axis_z = 2;
  const auto& axis_vectors = GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim>::axis_vectors;
  const GALS::CPU::Vec3<int> base_i_j_k = grid.baseNodeId(x_interp);
  const GALS::CPU::Vec3<int> base_ip1_j_k =
      GALS::CPU::Vec3<int>(base_i_j_k[0] + axis_vectors(0, 0), base_i_j_k[1], base_i_j_k[2]);
  const GALS::CPU::Vec3<int> base_i_jp1_k =
      GALS::CPU::Vec3<int>(base_i_j_k[0], base_i_j_k[1] + axis_vectors(1, 1), base_i_j_k[2]);
  const GALS::CPU::Vec3<int> base_i_j_kp1 =
      GALS::CPU::Vec3<int>(base_i_j_k[0], base_i_j_k[1], base_i_j_k[2] + axis_vectors(2, 2));
  const GALS::CPU::Vec3<int> base_ip1_jp1_k =
      GALS::CPU::Vec3<int>(base_i_j_k[0] + axis_vectors(0, 0), base_i_j_k[1] + axis_vectors(1, 1), base_i_j_k[2]);
  const GALS::CPU::Vec3<int> base_i_jp1_kp1 =
      GALS::CPU::Vec3<int>(base_i_j_k[0], base_i_j_k[1] + axis_vectors(1, 1), base_i_j_k[2] + axis_vectors(2, 2));
  const GALS::CPU::Vec3<int> base_ip1_j_kp1 =
      GALS::CPU::Vec3<int>(base_i_j_k[0] + axis_vectors(0, 0), base_i_j_k[1], base_i_j_k[2] + axis_vectors(2, 2));
  const GALS::CPU::Vec3<int> base_ip1_jp1_kp1 = GALS::CPU::Vec3<int>(
      base_i_j_k[0] + axis_vectors(0, 0), base_i_j_k[1] + axis_vectors(1, 1), base_i_j_k[2] + axis_vectors(2, 2));

  const typename T_GRID::position_type x_base = grid(base_i_j_k);
  const auto& dx = grid.dX();
  const auto& one_over_dx = grid.oneOverDX();

  GALS::CPU::Vec3<T> eta = (x_interp - x_base) * one_over_dx;

  const ControlPoints<T>& control_points_x = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_i_j_k), levelset.psiPrev()(base_i_j_k)[axis_x], levelset.phiPrev()(base_ip1_j_k),
      levelset.psiPrev()(base_ip1_j_k)[axis_x], dx[axis_x], use_gradient_limiting);
  const ControlPoints<T>& control_points_y = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_i_j_k), levelset.psiPrev()(base_i_j_k)[axis_y], levelset.phiPrev()(base_i_jp1_k),
      levelset.psiPrev()(base_i_jp1_k)[axis_y], dx[axis_y], use_gradient_limiting);
  const ControlPoints<T>& control_points_z = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_i_j_k), levelset.psiPrev()(base_i_j_k)[axis_z], levelset.phiPrev()(base_i_j_kp1),
      levelset.psiPrev()(base_i_j_kp1)[axis_z], dx[axis_z], use_gradient_limiting);
  const ControlPoints<T>& control_points_xy = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_ip1_j_k), levelset.psiPrev()(base_ip1_j_k)[axis_y], levelset.phiPrev()(base_ip1_jp1_k),
      levelset.psiPrev()(base_ip1_jp1_k)[axis_y], dx[axis_y], use_gradient_limiting);
  const ControlPoints<T>& control_points_yx = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_i_jp1_k), levelset.psiPrev()(base_i_jp1_k)[axis_x], levelset.phiPrev()(base_ip1_jp1_k),
      levelset.psiPrev()(base_ip1_jp1_k)[axis_x], dx[axis_x], use_gradient_limiting);
  const ControlPoints<T>& control_points_xz = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_ip1_j_k), levelset.psiPrev()(base_ip1_j_k)[axis_z], levelset.phiPrev()(base_ip1_j_kp1),
      levelset.psiPrev()(base_ip1_j_kp1)[axis_z], dx[axis_z], use_gradient_limiting);
  const ControlPoints<T>& control_points_zxy = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_ip1_j_kp1), levelset.psiPrev()(base_ip1_j_kp1)[axis_y],
      levelset.phiPrev()(base_ip1_jp1_kp1), levelset.psiPrev()(base_ip1_jp1_kp1)[axis_y], dx[axis_y],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_xyz = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_ip1_jp1_k), levelset.psiPrev()(base_ip1_jp1_k)[axis_z],
      levelset.phiPrev()(base_ip1_jp1_kp1), levelset.psiPrev()(base_ip1_jp1_kp1)[axis_z], dx[axis_z],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_zx = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_i_j_kp1), levelset.psiPrev()(base_i_j_kp1)[axis_x], levelset.phiPrev()(base_ip1_j_kp1),
      levelset.psiPrev()(base_ip1_j_kp1)[axis_x], dx[axis_x], use_gradient_limiting);
  const ControlPoints<T>& control_points_zy = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_i_j_kp1), levelset.psiPrev()(base_i_j_kp1)[axis_y], levelset.phiPrev()(base_i_jp1_kp1),
      levelset.psiPrev()(base_i_jp1_kp1)[axis_y], dx[axis_y], use_gradient_limiting);
  const ControlPoints<T>& control_points_yzx = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_i_jp1_kp1), levelset.psiPrev()(base_i_jp1_kp1)[axis_x],
      levelset.phiPrev()(base_ip1_jp1_kp1), levelset.psiPrev()(base_ip1_jp1_kp1)[axis_x], dx[axis_x],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_yz = GALS::INTERPOLATION::get_control_points(
      levelset.phiPrev()(base_i_jp1_k), levelset.psiPrev()(base_i_jp1_k)[axis_z], levelset.phiPrev()(base_i_jp1_kp1),
      levelset.psiPrev()(base_i_jp1_kp1)[axis_z], dx[axis_z], use_gradient_limiting);
  const ControlPoints<T>& control_points_y_z = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_j_k)[axis_z], levelset.phiMixedDerivativesPrev()(base_i_j_k)[1],
      levelset.psiPrev()(base_i_jp1_k)[axis_z], levelset.phiMixedDerivativesPrev()(base_i_jp1_k)[1], dx[axis_y],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_zy_z = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_j_kp1)[axis_z], levelset.phiMixedDerivativesPrev()(base_i_j_kp1)[1],
      levelset.psiPrev()(base_i_jp1_kp1)[axis_z], levelset.phiMixedDerivativesPrev()(base_i_jp1_kp1)[1], dx[axis_y],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_z_x = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_j_k)[axis_x], levelset.phiMixedDerivativesPrev()(base_i_j_k)[2],
      levelset.psiPrev()(base_i_j_kp1)[axis_x], levelset.phiMixedDerivativesPrev()(base_i_j_kp1)[2], dx[axis_z],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_x_y = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_j_k)[axis_y], levelset.phiMixedDerivativesPrev()(base_i_j_k)[0],
      levelset.psiPrev()(base_ip1_j_k)[axis_y], levelset.phiMixedDerivativesPrev()(base_ip1_j_k)[0], dx[axis_x],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_z_xy = GALS::INTERPOLATION::get_control_points(
      levelset.phiMixedDerivativesPrev()(base_i_j_k)[0], levelset.phiMixedDerivativesPrev()(base_i_j_k)[3],
      levelset.phiMixedDerivativesPrev()(base_i_j_kp1)[0], levelset.phiMixedDerivativesPrev()(base_i_j_kp1)[3],
      dx[axis_z], use_gradient_limiting);
  const ControlPoints<T>& control_points_z_y = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_j_k)[axis_y], levelset.phiMixedDerivativesPrev()(base_i_j_k)[1],
      levelset.psiPrev()(base_i_j_kp1)[axis_y], levelset.phiMixedDerivativesPrev()(base_i_j_kp1)[1], dx[axis_z],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_x_z = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_j_k)[axis_z], levelset.phiMixedDerivativesPrev()(base_i_j_k)[2],
      levelset.psiPrev()(base_ip1_j_k)[axis_z], levelset.phiMixedDerivativesPrev()(base_ip1_j_k)[2], dx[axis_x],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_zx_z = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_j_kp1)[axis_z], levelset.phiMixedDerivativesPrev()(base_i_j_kp1)[2],
      levelset.psiPrev()(base_ip1_j_kp1)[axis_z], levelset.phiMixedDerivativesPrev()(base_ip1_j_kp1)[2], dx[axis_x],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_zx_y = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_j_kp1)[axis_y], levelset.phiMixedDerivativesPrev()(base_i_j_kp1)[0],
      levelset.psiPrev()(base_ip1_j_kp1)[axis_y], levelset.phiMixedDerivativesPrev()(base_ip1_j_kp1)[0], dx[axis_x],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_yx_y = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_jp1_k)[axis_y], levelset.phiMixedDerivativesPrev()(base_i_jp1_k)[0],
      levelset.psiPrev()(base_ip1_jp1_k)[axis_y], levelset.phiMixedDerivativesPrev()(base_ip1_jp1_k)[0], dx[axis_x],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_yz_xy = GALS::INTERPOLATION::get_control_points(
      levelset.phiMixedDerivativesPrev()(base_i_jp1_k)[0], levelset.phiMixedDerivativesPrev()(base_i_jp1_k)[3],
      levelset.phiMixedDerivativesPrev()(base_i_jp1_kp1)[0], levelset.phiMixedDerivativesPrev()(base_i_jp1_kp1)[3],
      dx[axis_z], use_gradient_limiting);
  const ControlPoints<T>& control_points_yz_y = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_jp1_k)[axis_y], levelset.phiMixedDerivativesPrev()(base_i_jp1_k)[1],
      levelset.psiPrev()(base_i_jp1_kp1)[axis_y], levelset.phiMixedDerivativesPrev()(base_i_jp1_kp1)[1], dx[axis_z],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_yx_z = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_jp1_k)[axis_z], levelset.phiMixedDerivativesPrev()(base_i_jp1_k)[2],
      levelset.psiPrev()(base_ip1_jp1_k)[axis_z], levelset.phiMixedDerivativesPrev()(base_ip1_jp1_k)[2], dx[axis_x],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_yzx_y = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_jp1_kp1)[axis_y], levelset.phiMixedDerivativesPrev()(base_i_jp1_kp1)[0],
      levelset.psiPrev()(base_ip1_jp1_kp1)[axis_y], levelset.phiMixedDerivativesPrev()(base_ip1_jp1_kp1)[0], dx[axis_x],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_yzx_z = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_i_jp1_kp1)[axis_z], levelset.phiMixedDerivativesPrev()(base_i_jp1_kp1)[2],
      levelset.psiPrev()(base_ip1_jp1_kp1)[axis_z], levelset.phiMixedDerivativesPrev()(base_ip1_jp1_kp1)[2], dx[axis_x],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_xz_xy = GALS::INTERPOLATION::get_control_points(
      levelset.phiMixedDerivativesPrev()(base_ip1_j_k)[0], levelset.phiMixedDerivativesPrev()(base_ip1_j_k)[3],
      levelset.phiMixedDerivativesPrev()(base_ip1_j_kp1)[0], levelset.phiMixedDerivativesPrev()(base_ip1_j_kp1)[3],
      dx[axis_z], use_gradient_limiting);
  const ControlPoints<T>& control_points_xz_y = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_ip1_j_k)[axis_y], levelset.phiMixedDerivativesPrev()(base_ip1_j_k)[1],
      levelset.psiPrev()(base_ip1_j_kp1)[axis_y], levelset.phiMixedDerivativesPrev()(base_ip1_j_kp1)[1], dx[axis_z],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_xyz_xy = GALS::INTERPOLATION::get_control_points(
      levelset.phiMixedDerivativesPrev()(base_ip1_jp1_k)[0], levelset.phiMixedDerivativesPrev()(base_ip1_jp1_k)[3],
      levelset.phiMixedDerivativesPrev()(base_ip1_jp1_kp1)[0], levelset.phiMixedDerivativesPrev()(base_ip1_jp1_kp1)[3],
      dx[axis_z], use_gradient_limiting);
  const ControlPoints<T>& control_points_xyz_y = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_ip1_jp1_k)[axis_y], levelset.phiMixedDerivativesPrev()(base_ip1_jp1_k)[1],
      levelset.psiPrev()(base_ip1_jp1_kp1)[axis_y], levelset.phiMixedDerivativesPrev()(base_ip1_jp1_kp1)[1], dx[axis_z],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_xy_z = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_ip1_j_k)[axis_z], levelset.phiMixedDerivativesPrev()(base_ip1_j_k)[1],
      levelset.psiPrev()(base_ip1_jp1_k)[axis_z], levelset.phiMixedDerivativesPrev()(base_ip1_jp1_k)[1], dx[axis_y],
      use_gradient_limiting);
  const ControlPoints<T>& control_points_zxy_z = GALS::INTERPOLATION::get_control_points(
      levelset.psiPrev()(base_ip1_j_kp1)[axis_z], levelset.phiMixedDerivativesPrev()(base_ip1_j_kp1)[1],
      levelset.psiPrev()(base_ip1_jp1_kp1)[axis_z], levelset.phiMixedDerivativesPrev()(base_ip1_jp1_kp1)[1], dx[axis_y],
      use_gradient_limiting);

  const T& c_21_x = control_points_x.c_21;
  const T& c_12_x = control_points_x.c_12;
  const T& c_21_y = control_points_y.c_21;
  const T& c_12_y = control_points_y.c_12;
  const T& c_30_z = control_points_z.c_30;
  const T& c_21_z = control_points_z.c_21;
  const T& c_12_z = control_points_z.c_12;
  const T& c_03_z = control_points_z.c_03;
  const T& c_21_xy = control_points_xy.c_21;
  const T& c_12_xy = control_points_xy.c_12;
  const T& c_21_yx = control_points_yx.c_21;
  const T& c_12_yx = control_points_yx.c_12;
  const T& c_30_xz = control_points_xz.c_30;
  const T& c_21_xz = control_points_xz.c_21;
  const T& c_12_xz = control_points_xz.c_12;
  const T& c_03_xz = control_points_xz.c_03;
  const T& c_21_zxy = control_points_zxy.c_21;
  const T& c_12_zxy = control_points_zxy.c_12;
  const T& c_30_xyz = control_points_xyz.c_30;
  const T& c_21_xyz = control_points_xyz.c_21;
  const T& c_12_xyz = control_points_xyz.c_12;
  const T& c_03_xyz = control_points_xyz.c_03;
  const T& c_21_zx = control_points_zx.c_21;
  const T& c_12_zx = control_points_zx.c_12;
  const T& c_21_zy = control_points_zy.c_21;
  const T& c_12_zy = control_points_zy.c_12;
  const T& c_21_yzx = control_points_yzx.c_21;
  const T& c_12_yzx = control_points_yzx.c_12;
  const T& c_30_yz = control_points_yz.c_30;
  const T& c_21_yz = control_points_yz.c_21;
  const T& c_12_yz = control_points_yz.c_12;
  const T& c_03_yz = control_points_yz.c_03;
  const T& c_21_y_z = control_points_y_z.c_21;
  const T& c_12_y_z = control_points_y_z.c_12;
  const T& c_21_zy_z = control_points_zy_z.c_21;
  const T& c_12_zy_z = control_points_zy_z.c_12;
  const T& c_21_z_x = control_points_z_x.c_21;
  const T& c_12_z_x = control_points_z_x.c_12;
  const T& c_21_x_y = control_points_x_y.c_21;
  const T& c_12_x_y = control_points_x_y.c_12;
  const T& c_21_z_xy = control_points_z_xy.c_21;
  const T& c_12_z_xy = control_points_z_xy.c_12;
  const T& c_21_z_y = control_points_z_y.c_21;
  const T& c_12_z_y = control_points_z_y.c_12;
  const T& c_21_x_z = control_points_x_z.c_21;
  const T& c_12_x_z = control_points_x_z.c_12;
  const T& c_21_zx_z = control_points_zx_z.c_21;
  const T& c_12_zx_z = control_points_zx_z.c_12;
  const T& c_21_zx_y = control_points_zx_y.c_21;
  const T& c_12_zx_y = control_points_zx_y.c_12;
  const T& c_21_yx_y = control_points_yx_y.c_21;
  const T& c_12_yx_y = control_points_yx_y.c_12;
  const T& c_21_yz_xy = control_points_yz_xy.c_21;
  const T& c_12_yz_xy = control_points_yz_xy.c_12;
  const T& c_21_yz_y = control_points_yz_y.c_21;
  const T& c_12_yz_y = control_points_yz_y.c_12;
  const T& c_21_yx_z = control_points_yx_z.c_21;
  const T& c_12_yx_z = control_points_yx_z.c_12;
  const T& c_21_yzx_y = control_points_yzx_y.c_21;
  const T& c_12_yzx_y = control_points_yzx_y.c_12;
  const T& c_21_yzx_z = control_points_yzx_z.c_21;
  const T& c_12_yzx_z = control_points_yzx_z.c_12;
  const T& c_21_xz_xy = control_points_xz_xy.c_21;
  const T& c_12_xz_xy = control_points_xz_xy.c_12;
  const T& c_21_xz_y = control_points_xz_y.c_21;
  const T& c_12_xz_y = control_points_xz_y.c_12;
  const T& c_21_xyz_xy = control_points_xyz_xy.c_21;
  const T& c_12_xyz_xy = control_points_xyz_xy.c_12;
  const T& c_21_xyz_y = control_points_xyz_y.c_21;
  const T& c_12_xyz_y = control_points_xyz_y.c_12;
  const T& c_21_xy_z = control_points_xy_z.c_21;
  const T& c_12_xy_z = control_points_xy_z.c_12;
  const T& c_21_zxy_z = control_points_zxy_z.c_21;
  const T& c_12_zxy_z = control_points_zxy_z.c_12;

  // HERE
  const T& hx_by_three = dx[axis_x] * (T)one_third;
  const T& hy_by_three = dx[axis_y] * (T)one_third;
  const T& hz_by_three = dx[axis_z] * (T)one_third;

  T bx[] = {B0(eta[0]), B1(eta[0]), B2(eta[0]), B3(eta[0])};
  T by[] = {B0(eta[1]), B1(eta[1]), B2(eta[1]), B3(eta[1])};
  T bz[] = {B0(eta[2]), B1(eta[2]), B2(eta[2]), B3(eta[2])};
  T bx_prime[] = {B0_Prime(eta[0]), B1_Prime(eta[0]), B2_Prime(eta[0]), B3_Prime(eta[0])};
  T by_prime[] = {B0_Prime(eta[1]), B1_Prime(eta[1]), B2_Prime(eta[1]), B3_Prime(eta[1])};
  T bz_prime[] = {B0_Prime(eta[2]), B1_Prime(eta[2]), B2_Prime(eta[2]), B3_Prime(eta[2])};

  T control_points_all[] = {
      c_30_z,
      c_21_z,
      c_12_z,
      c_03_z,
      c_21_y,
      c_21_y + (dx[axis_z] * (T)one_third * c_21_y_z),
      c_21_zy - (dx[axis_z] * (T)one_third * c_21_zy_z),
      c_21_zy,
      c_12_y,
      c_12_y + (dx[axis_z] * (T)one_third * c_12_y_z),
      c_12_zy - (dx[axis_z] * (T)one_third * c_12_zy_z),
      c_12_zy,
      c_30_yz,
      c_21_yz,
      c_12_yz,
      c_03_yz,
      c_21_x,
      c_21_z + (dx[axis_x] * (T)one_third * c_21_z_x),
      c_12_z + (dx[axis_x] * (T)one_third * c_12_z_x),
      c_21_zx,
      c_21_x + (dx[axis_y] * (T)one_third * c_21_x_y),
      c_21_x + (hx_by_three * hy_by_three * c_21_z_xy) + (hy_by_three * c_21_z_y) + (hz_by_three * c_21_x_z),
      c_21_zx + (hx_by_three * hy_by_three * c_12_z_xy) + (hy_by_three * c_12_z_y) - (hz_by_three * c_21_zx_z),
      c_21_zx + (hy_by_three * c_21_zx_y),
      c_21_yx - (hy_by_three * c_21_yx_y),
      c_21_yx - (hx_by_three * hy_by_three * c_21_yz_xy) - (hy_by_three * c_21_yz_y) + (hz_by_three * c_21_yx_z),
      c_21_yzx - (hx_by_three * hy_by_three * c_12_yz_xy) - (hy_by_three * c_12_yz_y) - (hz_by_three * c_21_yzx_z),
      c_21_yzx - (hy_by_three * c_21_yzx_y),
      c_21_yx,
      c_21_yx + (hz_by_three * c_21_yx_z),
      c_21_yzx - (hz_by_three * c_21_yzx_z),
      c_21_yzx,
      c_12_x,
      c_12_x + (hz_by_three * c_12_x_z),
      c_12_zx - (hz_by_three * c_12_zx_z),
      c_12_zx,
      c_12_x + (hy_by_three * c_12_x_y),
      c_12_x - (hx_by_three * hy_by_three * c_21_xz_xy) + (hy_by_three * c_21_xz_y) + (hz_by_three * c_12_x_z),
      c_12_zx - (hx_by_three * hy_by_three * c_12_xz_xy) + (hy_by_three * c_12_xz_y) - (hz_by_three * c_12_zx_z),
      c_12_zx + (hy_by_three * c_12_zx_y),
      c_12_yx - (hy_by_three * c_12_yx_y),
      c_12_yx + (hx_by_three * hy_by_three * c_21_xyz_xy) - (hy_by_three * c_21_xyz_y) + (hz_by_three * c_12_yx_z),
      c_12_yzx + (hx_by_three * hy_by_three * c_12_xyz_xy) - (hy_by_three * c_12_xyz_y) - (hz_by_three * c_12_yzx_z),
      c_12_yzx - (hy_by_three * c_12_yzx_y),
      c_12_yx,
      c_12_yx + (hz_by_three * c_12_yx_z),
      c_12_yzx - (hz_by_three * c_12_yzx_z),
      c_12_yzx,
      c_30_xz,
      c_21_xz,
      c_12_xz,
      c_03_xz,
      c_21_xy,
      c_21_xy + (hz_by_three * c_21_xy_z),
      c_21_zxy - (hz_by_three * c_21_zxy_z),
      c_21_zxy,
      c_12_xy,
      c_12_xy + (hz_by_three * c_12_xy_z),
      c_12_zxy - (hz_by_three * c_12_zxy_z),
      c_12_zxy,
      c_30_xyz,
      c_21_xyz,
      c_12_xyz,
      c_03_xyz};

  for (int i = 0, l = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k) {
        hermite_fields.phi_interpolated += bx[i] * by[j] * bz[k] * control_points_all[l];
        ++l;
      }

  for (int i = 0, l = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k) {
        hermite_fields.psi_interpolated[axis_x] += bx_prime[i] * by[j] * bz[k] * control_points_all[l];
        ++l;
      }
  hermite_fields.psi_interpolated[axis_x] *= one_over_dx[axis_x];

  for (int i = 0, l = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k) {
        hermite_fields.psi_interpolated[axis_y] += bx[i] * by_prime[j] * bz[k] * control_points_all[l];
        ++l;
      }
  hermite_fields.psi_interpolated[axis_y] *= one_over_dx[axis_y];

  for (int i = 0, l = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k) {
        hermite_fields.psi_interpolated[axis_z] += bx[i] * by[j] * bz_prime[k] * control_points_all[l];
        ++l;
      }
  hermite_fields.psi_interpolated[axis_z] *= one_over_dx[axis_z];

  return hermite_fields;
}

template <typename T>
void GALS::INTERPOLATION::Hermite<T, GALS::CPU::Grid<T, 3>>::compute(
    const GALS::CPU::Array<GALS::CPU::Grid<T, 3>, typename GALS::CPU::Grid<T, 3>::position_type>& x_interp,
    GALS::CPU::Levelset<GALS::CPU::Grid<T, 3>, T>& levelset)
{
  typedef GALS::CPU::Grid<T, 3> T_GRID;

  const GALS::CPU::Vec3<int> num_cells_interp = x_interp.numCells();
  const T_GRID& grid = levelset.grid();
  const GALS::CPU::Vec3<typename T_GRID::value_type> dx = grid.dX();
  const auto& axis_vectors = GALS::CPU::Grid<typename T_GRID::value_type, T_GRID::dim>::axis_vectors;

  for (int i = 0; i < num_cells_interp[0]; ++i)
    for (int j = 0; j < num_cells_interp[1]; ++j)
      for (int k = 0; k < num_cells_interp[2]; ++k) {
        const auto& hermite_fields = this->interpolate(grid, x_interp(i, j, k), levelset);

        levelset.phiInterpPrev()(i, j, k) = hermite_fields.phi_interpolated;
        levelset.psiInterpPrev()(i, j, k) = hermite_fields.psi_interpolated;
      }
}

template class GALS::INTERPOLATION::Hermite<double, GALS::CPU::Grid<double, 3>>;
