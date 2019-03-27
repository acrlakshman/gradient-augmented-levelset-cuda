/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#pragma once

#include <math.h>

namespace GALS
{
static double VSMALL = 1e-10;

/*! Check for equality within a tolerance.
 *
 * \return number of elements.
 */
template <typename T>
bool is_equal(T a, T b)
{
  return fabs(a - b) < static_cast<T>(VSMALL);
}

}  // namespace GALS
