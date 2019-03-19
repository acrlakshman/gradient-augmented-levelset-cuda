/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include <vector>

namespace GALS {
namespace CPU {

template <typename T, int DIM = 3>
class Grid {
 public:
  typedef T value_type;

  Grid(int nx, int ny, int nz) : m_dimension(DIM), m_nx(nx), m_ny(ny), m_nz(nz) {
    m_grid.resize((m_nx + 2) * (m_ny + 2) * (m_nz + 2));
  }

  Grid(int nx, int ny) : Grid(nx, ny, 1) {}

  Grid(int nx) : Grid(nx, 1, 1) {}

  ~Grid() {}

  std::vector<T> x(int i, int j = 1, int k = 1) {
    std::vector<T> pos(m_dimension);

    for (int c = 0; c < m_dimension; ++c) pos[i] = m_grid[((k * m_nx * m_ny) + (j * m_nx) + i) * m_dimension + c];

    return pos;
  }

  const int dimension() const { return m_dimension; }

  void generate(T x_min, T x_max, T y_min, T y_max, T z_min, T z_max) {
    // TODO
  }

 private:
  int m_dimension, m_nx, m_ny, m_nz;
  std::vector<T> m_grid;
};

}  // namespace CPU
}  // namespace GALS
