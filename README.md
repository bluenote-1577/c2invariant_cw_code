# c2invariant_cw_code

Code for computing c_2 invariants of completed primitive graphs in phi4 theory.

## Requirements

* [Giac : C++ library for symbolic computation](https://www-fourier.ujf-grenoble.fr/~parisse/giac.html) - The code was tested with version 1.4.9-59. The linear algebra library was used to compute symbolic determinants for Dodgson polynomials. The instructions for installing the library are [here](https://www-fourier.ujf-grenoble.fr/~parisse/giac_compile.html).
* C++ compiler - We compile with -std=c++0x. gcc version 4.4.7 was used.
* OpenMP - The OpenMP library that was packaged with gcc version 4.4.7 was used.
* Maple - Maple 18 was used but presumably any version should work.

## Usage
We require the graphs to be in the same format as in the Periods file of https://arxiv.org/abs/1603.04289.

1. Change the number of threads used by changing `const int NUM_THREADS = x` in server_c2.cpp.
2. `cd` into the directory containing this repository and `make gen`
3. `./gen (graphsfile) (number from 0-31)` to generate the edge sequence.
4. `maple make_recipe.m` to generate the recipe file.
5. `make c2` 
6. `./c2cw.out (recipe file)`
7. c2 will be calculated up to p=31. Results are written in "c2.txt".
