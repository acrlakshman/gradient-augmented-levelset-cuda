Overview {#mainpage}
========

[TOC]

## About

Implementation of Gradient Augmented Levelset method in CPU and GPU.

## Source code

[Gradient Augmented Levelset Implementation in CPU & GPU][src]

## Installation

#### Dependencies

* [CMake]
* Compiler that supports C++17

```sh
   git clone https://github.com/acrlakshman/gradient-augmented-levelset-cuda --recursive
   cd gradient-augmented-levelset-cuda
   mkdir -p build && cd build
   cmake .. && make -j 4
```

### Additional build options

#### Dependencies

* [Doxygen]
* [lcov]

```sh
cmake .. -DBUILD_DOCUMENTATION=ON -DBUILD_COVERAGE=ON
make -j 4
./tests/gals_unit_tests
lcov --directory . --base-directory ../src --capture --no-external --output-file coverage.info
genhtml coverage.info --output-directory ./docs/html/coverage
```

* Documentation can be found at `./docs/html/index.html`.

## Coverage

[Coverage]

### Primary developers

* [Lakshman Anumolu][Lakshman] (acrlakshman@yahoo.co.in)
* [Raunak Bardia][Raunak] (raunakbardia@gmail.com)

[src]:https://github.com/acrlakshman/gradient_augmented_levelset_cuda
[CMake]:https://github.com/Kitware/CMake
[Doxygen]:https://github.com/doxygen/doxygen
[lcov]:https://github.com/linux-test-project/lcov
[Coverage]:https://acrlakshman.github.io/gradient-augmented-levelset-cuda/coverage
[Lakshman]:https://lakshmananumolu.com
[Raunak]:https://raunakbardia.wordpress.com
