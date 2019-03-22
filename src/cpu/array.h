/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#pragma once

#include "grid.h"

#include <vector>

namespace GALS {
namespace CPU {

template <typename T_GRID, typename T_ARRAY>
class Array {
 public:
  Array(Grid<typename T_GRID::value_type, T_GRID::dim> &grid)
      : m_grid(grid),
        m_nx(grid.getNumCells()[0]),
        m_ny(grid.getNumCells()[1]),
        m_nz(grid.getNumCells()[2]),
        m_pad(grid.getPadding()) {
    m_data.resize(m_grid.size());

    for (int i = 0; i < m_data.size(); ++i) {
      m_data[i] = T_ARRAY();
    }
  }

  ~Array() {
    m_data.clear();
    m_data.shrink_to_fit();
  }

  const std::size_t size() const { return m_data.size(); }

  const T_ARRAY &operator[](const std::size_t idx) const { return m_data[idx]; }

  T_ARRAY &operator[](const std::size_t idx) { return m_data[idx]; }

  const T_ARRAY operator()(int i, int j, int k) {
    const std::size_t idx = m_grid.getIndex(i, j, k);
    return m_data[idx];
  }

 private:
  const T_GRID &m_grid;
  const int m_nx, m_ny, m_nz, m_pad;
  std::vector<T_ARRAY> m_data;
};

}  // namespace CPU
}  // namespace GALS
