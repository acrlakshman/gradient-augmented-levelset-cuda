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

template <typename T_GRID>
class Array {
 public:
  typedef typename T_GRID::value_type value_type;
  typedef value_type T;

  Array(Grid<T, T_GRID::dim> &grid) : m_grid(grid) {
  }

  ~Array() {}

 private:
  T_GRID &m_grid;
  T *m_data;
};

}
}
