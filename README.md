## Gradient Augmented Level Set Method - CPU & CUDA
[![Build Status](https://travis-ci.org/acrlakshman/gradient-augmented-levelset-cuda.svg?branch=master)](https://travis-ci.org/acrlakshman/gradient-augmented-levelset-cuda)

* Implementation of Gradient Augmented Level Set method in both CPU and GPU (using CUDA).

### Build instructions

```sh
git clone https://github.com/acrlakshman/gradient-augmented-levelset-cuda --recursive
cd gradient-augmented-levelset-cuda
mkdir -p build && cd build
cmake .. -GNinja && ninja -j 4
```


License
-------

BSD 2-Clause License. Please check the accompanying [LICENSE].txt file

Acknowledgements
----------------

[Pradeep Garigipati] for helping with cmake.

[Ryan Krattiger] for helping with cmake.

[LICENSE]:https://github.com/acrlakshman/gradient_augmented_levelset_cuda/blob/master/LICENSE
[Pradeep Garigipati]:https://github.com/9prady9
[Ryan Krattiger]:https://github.com/rjk9w5
