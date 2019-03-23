/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include "vec_n.h"

template <typename T, int SIZE>
GALS::CPU::VecN<T, SIZE>::VecN() {
  m_data.resize(SIZE);

  for (int i = 0; i < SIZE; ++i) m_data[i] = static_cast<T>(0);
}

template <typename T, int SIZE>
GALS::CPU::VecN<T, SIZE>::~VecN() {
  m_data.clear();
  m_data.shrink_to_fit();
}

template <typename T, int SIZE>
const int GALS::CPU::VecN<T, SIZE>::size() const {
  return m_data.size();
}

template <typename T, int SIZE>
const T GALS::CPU::VecN<T, SIZE>::operator[](const int idx) const {
  return m_data[idx];
}

template <typename T, int SIZE>
T& GALS::CPU::VecN<T, SIZE>::operator[](const int idx) {
  return m_data[idx];
}

template class GALS::CPU::VecN<int, 1>;
template class GALS::CPU::VecN<int, 2>;
template class GALS::CPU::VecN<int, 3>;
template class GALS::CPU::VecN<double, 1>;
template class GALS::CPU::VecN<double, 2>;
template class GALS::CPU::VecN<double, 3>;
