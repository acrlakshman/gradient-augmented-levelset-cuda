/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include "array.h"
#include "grid.h"

#include <iostream>

int main()
{
  std::cout << "Running CPU version of GALS implementation" << std::endl;

  int n_x = 10, n_y = 10;

  GALS::CPU::Grid<double, 2> grid(n_x, n_y);

  grid.generate(-1, 1, -1, 1, -1, 1);

  grid.writeToFile();
  return 0;
}
