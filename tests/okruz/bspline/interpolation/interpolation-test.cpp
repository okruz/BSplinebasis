#define BOOST_TEST_MODULE SplineTest
#include <okruz/bspline/interpolation/interpolation.h>
#include <okruz/bspline/support/Support.h>

#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#ifdef MYSPLINE_INTERPOLATION_USE_EIGEN
template <typename T, size_t order>
void testInterpolationEigen(T tol) {
  using Spline = okruz::bspline::Spline<T, order>;
  using Support = okruz::bspline::support::Support<T>;
  using Grid = okruz::bspline::support::Grid<T>;
  using Construction = okruz::bspline::support::Construction;

  const Grid grid(std::vector<T>{-3.0l, -2.5l, -1.5l, -1.0l, 0.0l, 0.5l, 1.5l,
                                 2.5l, 3.5l, 4.0l, 5.0l});
  const Support x(grid, Construction::WHOLE_GRID);
  const std::vector<T> y{-3.0l, -2.5l, -1.5l, -1.0l, 0.0l, -0.5l,
                         -1.5l, -2.5l, -3.5l, -4.0l, 3.0l};
  Spline s =
      okruz::bspline::interpolation::interpolate_using_eigen<T, order>(x, y);
  for (size_t i = 0; i < x.size(); i++) {
    BOOST_CHECK_SMALL(s(x[i]) - y[i], tol);
  }
}

BOOST_AUTO_TEST_CASE(TestInterpolationEigen) {
  testInterpolationEigen<double, 1>(2.0e-14);
  testInterpolationEigen<double, 2>(2.0e-14);
  testInterpolationEigen<double, 3>(2.0e-14);
  testInterpolationEigen<double, 4>(2.0e-14);

  if constexpr (sizeof(long double) != sizeof(double)) {
    testInterpolationEigen<long double, 1>(1.0e-17l);
    testInterpolationEigen<long double, 2>(1.0e-17l);
    testInterpolationEigen<long double, 3>(1.0e-17l);
    testInterpolationEigen<long double, 4>(1.0e-17l);
  }
}
#endif

#ifdef MYSPLINE_INTERPOLATION_USE_ARMADILLO
template <size_t order>
void testInterpolationArmadillo(double tol) {
  using Spline = okruz::bspline::Spline<double, order>;
  using Support = okruz::bspline::support::Support<double>;
  using Grid = okruz::bspline::support::Grid<double>;
  using Construction = okruz::bspline::support::Construction;

  const Grid grid(std::vector<double>{-3.0l, -2.5l, -1.5l, -1.0l, 0.0l, 0.5l,
                                      1.5l, 2.5l, 3.5l, 4.0l, 5.0l});
  const Support x(grid, Construction::WHOLE_GRID);
  const std::vector<double> y{-3.0l, -2.5l, -1.5l, -1.0l, 0.0l, -0.5l,
                              -1.5l, -2.5l, -3.5l, -4.0l, 3.0l};
  Spline s =
      okruz::bspline::interpolation::interpolate_using_armadillo<order>(x, y);
  for (size_t i = 0; i < x.size(); i++) {
    BOOST_CHECK_SMALL(s(x[i]) - y[i], tol);
  }
}

BOOST_AUTO_TEST_CASE(TestInterpolationArmadillo) {
  testInterpolationArmadillo<1>(2.0e-14);
  testInterpolationArmadillo<2>(2.0e-14);
  testInterpolationArmadillo<3>(2.0e-14);
  testInterpolationArmadillo<4>(2.0e-14);
}
#endif
