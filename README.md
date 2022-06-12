# BSplinebasis
This repository contains a small, template-based BSpline library in C++ (>= C++17), geared towards the use as basis functions in analytical problems. The main library can be found in the `include/` directory. To use the library,
all you have to do is to add this folder to the include path.

## Why BSplines?
### Advantages of using BSplines
BSplines have a number of advantagous properties for the use as basis functions. Some of these advantages are:

* They are very versatile and can be adapted to many different problems.
* The parametrization is very intuitive.
* They have a finite support, ensuring that each spline overlaps only with a few neighbour splines. They can therefore be used in sparse-matrix scenarios.
* Many integrals can easily be evaluated analytically which makes the calculations quite fast and reasonably accurate.

For a case study on the numerical properties of a BSpline basis, see [here](readme/convergence.md).

### Disadvantages of using BSplines

* They are only finitely many times continously differentiable. Make sure to also read the `Cautions` section.
* They are not orthogonal. Using the BSplines as basis functions in eigenvalue problems will usually result in generalized, algebraic eigenvalue problems of the kind A x = lambda B x (still self-adjoint, though).

## Usage of the library
### Generating BSplines
The first step to generate a basis set of BSplines is to define your knots. From the knots, a set of BSplines can be generated using the `bspline::BSplineGenerator`.
```C++
  #include <bspline/Core.h>
  
  static constexpr size_t SPLINE_ORDER = 3;
  const std::vector<double> knots{0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  bspline::BSplineGenerator gen(knots);
  const std::vector<bspline::Spline<double, SPLINE_ORDER>> = gen.template generateBSplines<SPLINE_ORDER + 1>();
```

The corresponding splines are shown in the following graphic.
![Third order BSplines.](readme/splines_normal.png?raw=true "Third order BSplines.")


## Dependencies
The **core library** does not have any additional dependencies beyond a C++ compiler supporting C++17. Everything that is directly or indirectly included by `include/bspline/Core.h` is considered part of the core library. Parts that are not part of the core library are:

* The contents of `include/bspline/interpolation/`.
  * The interpolation needs a linear algebra framework to solve the linear equation system arising during the interpolation. You can use your linear algebra framework of choice with relative ease, or use interpolation based on `armadillo` (only double precision) and `Eigen` provided by the library. For details, see the source or [docs](https://okruz.github.io/BSplinebasis/namespacebspline_1_1interpolation.html).
* The contents of `include/bspline/integration/numerical.h`.
  * The numerical integration routine implemented there is based on the Gauss-Legendre quadrature scheme provided by `Boost::math`.

Furthermore, the tests and examples require `Eigen` to compile and the tests are based on `Boost::test`.

The tests are currently only run regularly on an x64 linux plattform using gcc and clang. The main library should be usable with every standard-conformant C++ compiler supporting C++17. If you are using the library on a different plattform, I would be happy to receive your feedback.


## Cautions
The operators implemented in `include/bspline/operators/` assume that the result of their application to a spline can again be represented as a spline. As each spline is only finitely many times continously differentiable, this is not true for every application of a derivative operator to a spline. If the spline is not sufficiently many times continously differentiable, Dirac delta distributions may arise during the application, which are not handled at all by the library. It is left up to the user of the library to make sure that the forementioned assumption holds true.


## Docs
The Doxygen docs can be found [here](https://okruz.github.io/BSplinebasis/).

If you are not familiar with BSplines, a good place to start is the [Wikipedia article](https://en.wikipedia.org/wiki/B-spline).
