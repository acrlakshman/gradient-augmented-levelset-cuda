/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include "grid.h"
#include "array.h"

#include <iostream>

int main()
{
  std::cout << "Running CPU version of GALS implementation" << std::endl;

  GALS::CPU::Grid<double, 3> grid(4, 4, 4);

  std::cout << "grid.x = " << grid.x(1)[0] << std::endl;
  std::cout << "dimension = " << grid.dimension() << std::endl;

  // Array
  GALS::CPU::Array<GALS::CPU::Grid<double, 3>> levelset(grid);

  return 0;
}
