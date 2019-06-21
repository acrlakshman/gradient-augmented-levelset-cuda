## Gradient Augmented Level Set Method - CPU & CUDA
[![Build Status](https://travis-ci.org/acrlakshman/gradient-augmented-levelset-cuda.svg?branch=master)](https://travis-ci.org/acrlakshman/gradient-augmented-levelset-cuda)
[![Coverage Status](https://coveralls.io/repos/github/acrlakshman/gradient-augmented-levelset-cuda/badge.svg)](https://coveralls.io/github/acrlakshman/gradient-augmented-levelset-cuda)

* Implementation of Gradient Augmented Level Set method in both CPU and GPU (using CUDA).

### Build instructions

#### Dependencies

* [CMake]
* Compiler that supports C++11

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

### Documentation

* [Documentation]
* [Coverage]

License
-------

BSD 3-Clause License. Please check the accompanying [LICENSE] file

Acknowledgements
----------------

* [Pradeep Garigipati] for helping with cmake.
* [Ryan Krattiger] for helping with cmake.
* [John Van Gilder] for helping with coverage tools.

[CMake]:https://github.com/Kitware/CMake
[Doxygen]:https://github.com/doxygen/doxygen
[lcov]:https://github.com/linux-test-project/lcov
[Documentation]:https://acrlakshman.github.io/gradient-augmented-levelset-cuda
[Coverage]:https://acrlakshman.github.io/gradient-augmented-levelset-cuda/coverage
[LICENSE]:https://github.com/acrlakshman/gradient_augmented_levelset_cuda/blob/master/LICENSE
[Pradeep Garigipati]:https://github.com/9prady9
[Ryan Krattiger]:https://github.com/rjk9w5
[John Van Gilder]: https://github.com/JohnVanGilder
