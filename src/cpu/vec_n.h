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

template <typename T, int SIZE = 3>
class VecN {
 public:
  typedef T value_type;

  VecN() {
    m_data.resize(SIZE);
    for (int i = 0; i < SIZE; ++i) m_data[i] = static_cast<T>(0);
  }

  ~VecN() {
    m_data.clear();
    m_data.shrink_to_fit();
  }

  const int size() const { return m_data.size(); }

  const T operator[](const int idx) const { return m_data[idx]; }

  T &operator[](const int idx) { return m_data[idx]; }

 private:
  std::vector<T> m_data;
};

}  // namespace CPU
}  // namespace GALS
