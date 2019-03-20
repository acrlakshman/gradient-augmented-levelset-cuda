/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#pragma once

#include <vector>

namespace GALS {
namespace CPU {

template <typename T, int DIM = 3>
class Grid {
 public:
  typedef T value_type;
  static const int dim = DIM;

  Grid(int nx, int ny, int nz);

  Grid(int nx, int ny);

  Grid(int nx);

  ~Grid();

  std::vector<T> x(int i, int j = 1, int k = 1);

  const int dimension() const;

  const int size() const;

  const int* getMask() const;

  const std::vector<int> getNumCells() const;

  const int getPadding() const;

  const std::size_t getIndex(int i, int j, int k);

  void setPadding(const int pad);

  void generate(T x_min, T x_max, T y_min, T y_max, T z_min, T z_max);

 private:
  int m_dimension, m_nx, m_ny, m_nz, m_pad;
  int m_mask[3];
  std::vector<T> m_grid;
};

}  // namespace CPU
}  // namespace GALS
