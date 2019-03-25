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
/*! \class VecN
 *
 * Class to create varying size elements at a computational cell. For e.x. velocity, gradients, etc.
 */
template <typename T, int SIZE = 3>
class VecN
{
 public:
  typedef T value_type;

  /*! Default constructor
   */
  VecN();

  /*! Destructor
   */
  ~VecN();

  /*! Returns number of elements.
   *
   * \return number of elements.
   */
  const int size() const;

  /*! Overloaded subscript operator.
   *
   * \param idx zero based index of element.
   *
   * \return element at index (idx).
   */
  const T operator[](const int idx) const;

  /*! Overloaded subscript operator that returns a reference.
   *
   * \param idx zero based index of element.
   *
   * \return element at index (idx).
   */
  T &operator[](const int idx);

 private:
  std::vector<T> m_data;
};

}  // namespace CPU
}  // namespace GALS
