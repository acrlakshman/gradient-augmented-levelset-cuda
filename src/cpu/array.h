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
  Array(Grid<typename T_GRID::value_type, T_GRID::dim> &grid) : m_grid(grid), m_pad(grid.getPadding()) {
    m_num_cells = m_grid.getNumCells();
    m_data.resize(m_grid.size());

    for (int i = 0; i < m_data.size(); ++i) {
      m_data[i] = T_ARRAY();
    }
  }

  ~Array() {
    m_data.clear();
    m_data.shrink_to_fit();
  }

  const int size() const { return m_data.size(); }

  const T_ARRAY &operator[](int idx) const { return m_data[idx]; }

  T_ARRAY &operator[](int idx) { return m_data[idx]; }

  const T_ARRAY operator()(int i, int j, int k) {
    // TODO
    // return m_data[]
  }

 private:
  const T_GRID &m_grid;
  const int m_pad;
  std::vector<int> m_num_cells;
  std::vector<T_ARRAY> m_data;
};

}  // namespace CPU
}  // namespace GALS
