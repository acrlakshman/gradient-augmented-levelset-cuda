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

  Grid(int nx, int ny, int nz) : m_dimension(DIM), m_nx(nx), m_ny(ny), m_nz(nz), m_pad(1) {
    m_mask[0] = 1, m_mask[1] = 0, m_mask[2] = 0;

    if (dim == 2)
      m_mask[1] = 1;
    else if (dim == 3)
      m_mask[1] = 1, m_mask[2] = 1;
  }

  Grid(int nx, int ny) : Grid(nx, ny, 1) {}

  Grid(int nx) : Grid(nx, 1, 1) {}

  ~Grid() {
    m_grid.clear();
    m_grid.shrink_to_fit();
  }

  std::vector<T> x(int i, int j = 1, int k = 1) {
    std::vector<T> pos(m_dimension);

    for (int c = 0; c < m_dimension; ++c) pos[i] = m_grid[((k * m_nx * m_ny) + (j * m_nx) + i) * m_dimension + c];

    return pos;
  }

  const int* getMask() const { return m_mask; }

  const int dimension() const { return m_dimension; }

  const int size() const { return m_grid.size(); }

  const std::vector<int> getNumCells() const { return std::vector<int>{m_nx, m_ny, m_nz}; }

  const int getPadding() const { return m_pad; }

  const std::size_t getIndex(int i, int j, int k) {
    std::size_t idx = ((k + m_pad) * (m_nx + 2 * m_pad) * (m_ny + 2 * m_pad) * m_mask[2]) +
                      ((j + m_pad) * (m_nx + 2 * m_pad) * m_mask[1]) + (i + m_pad);

    return idx;
  }

  void setPadding(const int pad) { m_pad = pad; }

  void generate(T x_min, T x_max, T y_min, T y_max, T z_min, T z_max) {
    // TODO
  }

 private:
  int m_dimension, m_nx, m_ny, m_nz, m_pad;
  int m_mask[3];
  std::vector<T> m_grid;
  // Added Comment
};

}  // namespace CPU
}  // namespace GALS
