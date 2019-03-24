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

namespace GALS
{
namespace CPU
{
template <typename T, int SIZE = 3>
class VecN
{
 public:
  typedef T value_type;

  VecN();

  ~VecN();

  const int size() const;

  const T operator[](const int idx) const;

  T &operator[](const int idx);

 private:
  std::vector<T> m_data;
};

}  // namespace CPU
}  // namespace GALS
