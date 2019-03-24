/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#pragma once

#include "vec_n.h"

#include <vector>

namespace GALS
{
namespace CPU
{
template <typename T, int DIM = 3>
class Grid
{
 public:
  typedef T value_type;
  static const int dim = DIM;

  Grid(int nx, int ny, int nz);

  Grid(int nx, int ny);

  Grid(int nx);

  ~Grid();

  VecN<T, 3>& x(const int i, const int j = 1, const int k = 1);

  const int dimension() const;

  const int size() const;

  const int* getMask() const;

  const std::vector<int> getNumCells() const;

  const int getPadding() const;

  const std::size_t getIndex(const int i, const int j, const int k);

  const VecN<T, 3> dX() const;

  VecN<T, 3>& operator()(const int i, const int j, const int k);

  void setPadding(const int pad);

  void generate(T x_min, T x_max, T y_min, T y_max, T z_min, T z_max);

  void print(bool show_padding = false);

 private:
  int m_dimension, m_nx, m_ny, m_nz, m_pad;
  int m_mask[3];
  std::vector<VecN<T, 3>> m_grid;
  VecN<T, 3> m_dx;
};

}  // namespace CPU
}  // namespace GALS
