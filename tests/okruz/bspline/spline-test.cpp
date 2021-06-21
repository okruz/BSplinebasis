#define BOOST_TEST_MODULE SplineTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <okruz/bspline/BSplineGenerator.h>
#include <okruz/bspline/integration/analytical.h>
#include <okruz/bspline/integration/numerical.h>

using okruz::bspline::BSplineGenerator;
using okruz::bspline::Spline;

using namespace okruz::bspline::support;

template <typename T> Spline<T, 0> getOne(const Grid<T> &grid) {
  T onet = static_cast<T>(1);
  Support support(grid, Construction::WHOLE_GRID);
  std::vector<std::array<T, 1>> coeffs(support.size() - 1, {onet});
  return Spline<T, 0>(std::move(support), std::move(coeffs));
}

template <typename T, size_t order> void testIntegration(T tol) {
  using namespace okruz::bspline::integration;

  using Spline = okruz::bspline::Spline<T, order>;
  using Spline0 = okruz::bspline::Spline<T, 0>;
  const BSplineGenerator generator(std::vector<T>{
      -7.0l,  -6.85l, -6.55l, -6.3l, -6.0l, -5.75l, -5.53l, -5.2l,
      -4.75l, -4.5l,  -3.0l,  -2.5l, -1.5l, -1.0l,  0.0l,   0.5l,
      1.5l,   2.5l,   3.5l,   4.0l,  4.35l, 4.55l,  4.95l,  5.4l,
      5.7l,   6.1l,   6.35l,  6.5l,  6.85l, 7.0l});

  const std::vector<Spline> splines =
      generator.template generateBSplines<order + 1>();
  const Spline0 one = getOne(generator.getGrid());

  const auto f1 = [](const T &x) { return static_cast<T>(1); };
  const auto fx = [](const T &x) { return x; };

  for (const auto &s1 : splines) {
    for (const auto &s2 : splines) {
      auto s2dx2 = s2.dx2();
      BOOST_CHECK_SMALL(overlap<T>(s1, s2) - integrate<T>(s1 * s2), tol);
      BOOST_CHECK_SMALL(overlap<T>(s1, s2) - overlap<T>(one, s1 * s2), tol);
      BOOST_CHECK_SMALL(overlap<T>(s1, s2.timesx()) - integrate_x<T>(s1, s2),
                        tol);
      BOOST_CHECK_SMALL(overlap<T>(s1.timesx(), s2) - integrate_x<T>(s1, s2),
                        tol);
      BOOST_CHECK_SMALL(overlap<T>(s1, s2.timesx().timesx()) -
                            integrate_x2<T>(s1, s2),
                        static_cast<T>(5) * tol);
      BOOST_CHECK_SMALL(overlap<T>(s1.timesx(), s2.timesx()) -
                            integrate_x2<T>(s1, s2),
                        static_cast<T>(5) * tol);
      BOOST_CHECK_SMALL(overlap<T>(s1.timesx().timesx(), s2) -
                            integrate_x2<T>(s1, s2),
                        static_cast<T>(5) * tol);
      BOOST_CHECK_SMALL(overlap<T>(s1, s2.dx()) - integrate_dx<T>(s1, s2), tol);
      BOOST_CHECK_SMALL(
          overlap<T>(s1.timesx(), s2.dx()) - integrate_x_dx<T>(s1, s2), tol);
      BOOST_CHECK_SMALL(overlap<T>(s1, s2.dx().dx()) - integrate_dx2<T>(s1, s2),
                        tol);
      BOOST_CHECK_SMALL(overlap<T>(s1, s2dx2) - integrate_dx2<T>(s1, s2), tol);
      BOOST_CHECK_SMALL(overlap<T>(s1.timesx(), s2.dx().dx()) -
                            integrate_x_dx2<T>(s1, s2),
                        static_cast<T>(11) * tol);
      BOOST_CHECK_SMALL(overlap<T>(s1, s2.dx().dx().timesx()) -
                            integrate_x_dx2<T>(s1, s2),
                        static_cast<T>(8) * tol);
      BOOST_CHECK_SMALL(overlap<T>(s1.timesx().timesx(), s2.dx().dx()) -
                            integrate_x2_dx2<T>(s1, s2),
                        static_cast<T>(60) * tol);
      BOOST_CHECK_SMALL(overlap<T>(s1, s2.dx().dx().timesx().timesx()) -
                            integrate_x2_dx2<T>(s1, s2),
                        static_cast<T>(60) * tol);
      BOOST_CHECK_SMALL(overlap<T>(s1, s2dx2.timesx().timesx()) -
                            integrate_x2_dx2<T>(s1, s2),
                        static_cast<T>(60) * tol);
      BOOST_CHECK_SMALL(overlap<T>(s1, s2) - integrate<2 * order>(f1, s1, s2),
                        static_cast<T>(10) * tol);
      BOOST_CHECK_SMALL(integrate_x<T>(s1, s2) -
                            integrate<2 * order>(fx, s1, s2),
                        static_cast<T>(10) * tol);
    }

    auto s1_d_order = s1.template dx<order>();
    BOOST_TEST(!s1_d_order.isZero());

    auto s1_d_orderp1 = s1.template dx<order + 1>();
    BOOST_CHECK_SMALL(integrate<T>(s1_d_orderp1), tol);
    BOOST_TEST(s1_d_orderp1.isZero());

    auto s1_d_orderp2 = s1.template dx<order + 2>();
    BOOST_CHECK_SMALL(integrate<T>(s1_d_orderp2), tol);
    BOOST_TEST(s1_d_orderp2.isZero());
  }
}

BOOST_AUTO_TEST_CASE(TestIntegration) {
  constexpr double TOL = 1.0e-15;
  testIntegration<double, 2>(TOL);
  testIntegration<double, 3>(TOL);
  testIntegration<double, 4>(TOL);
  testIntegration<double, 5>(TOL);
  testIntegration<double, 6>(TOL);
  testIntegration<double, 7>(TOL);
  testIntegration<double, 8>(TOL);
  testIntegration<double, 9>(TOL);
  testIntegration<double, 10>(TOL);

  if constexpr (sizeof(long double) != sizeof(double)) {
    constexpr long double TOLL = 1.0e-18l;
    testIntegration<long double, 2>(TOLL);
    testIntegration<long double, 3>(TOLL);
    testIntegration<long double, 4>(TOLL);
    testIntegration<long double, 5>(TOLL);
    testIntegration<long double, 6>(TOLL);
    testIntegration<long double, 7>(TOLL);
    testIntegration<long double, 8>(TOLL);
    testIntegration<long double, 9>(TOLL);
    testIntegration<long double, 10>(TOLL);
  }
}

template <typename T, size_t order> void testArithmetic(T tol) {
  static_assert(order >= 2, "For this test, order must be at least 2");
  using Spline = okruz::bspline::Spline<T, order>;
  using Spline6 = okruz::bspline::Spline<T, 2 * order>;
  using Spline0 = okruz::bspline::Spline<T, 0>;
  using Spline1 = okruz::bspline::Spline<T, order - 2>;

  BSplineGenerator<T> generator(std::vector<T>{
      -7.0l,  -6.85l, -6.55l, -6.3l, -6.0l, -5.75l, -5.53l, -5.2l,
      -4.75l, -4.5l,  -3.0l,  -2.5l, -1.5l, -1.0l,  0.0l,   0.5l,
      1.5l,   2.5l,   3.5l,   4.0l,  4.35l, 4.55l,  4.95l,  5.4l,
      5.7l,   6.1l,   6.35l,  6.5l,  6.85l, 7.0l});

  const std::vector<Spline> splines =
      generator.template generateBSplines<order + 1>();
  const Spline0 one = getOne(generator.getGrid());

  for (const auto &s : splines) {
    Spline s2 = s * static_cast<T>(2);
    Spline sm = -s;
    Spline6 sprod = s * s;
    Spline s22 = s;
    s22 *= static_cast<T>(2);
    Spline shalf = s / static_cast<T>(2);
    Spline shalf2 = s;
    shalf2 /= static_cast<T>(2);
    Spline s5half = s2 + shalf;
    Spline s5half2 = s2;
    s5half2 += shalf;
    Spline s3half = s2 - shalf;
    Spline s3half2 = s2;
    s3half2 -= shalf;
    Spline splusone = s + one;
    Spline1 sdx2 = s.dx().dx();
    Spline1 sdx22 = s.dx2();
    Spline sdx0 = s.template dx<0>();

    for (T x = s.start(); x <= s.end(); x += 0.01L) {
      BOOST_CHECK_SMALL(sm(x) + s(x), tol);
      BOOST_CHECK_SMALL(s2(x) - static_cast<T>(2) * s(x), tol);
      BOOST_CHECK_SMALL(s22(x) - static_cast<T>(2) * s(x),
                        tol); // Tests *= operator
      BOOST_CHECK_SMALL(shalf(x) - s(x) / static_cast<T>(2),
                        tol); // Tests / operator
      BOOST_CHECK_SMALL(shalf2(x) - s(x) / static_cast<T>(2),
                        tol); // Tests /= operator
      BOOST_CHECK_SMALL(s5half(x) -
                            static_cast<T>(5) * s(x) / static_cast<T>(2),
                        tol); // Tests + operator
      BOOST_CHECK_SMALL(s5half2(x) -
                            static_cast<T>(5) * s(x) / static_cast<T>(2),
                        tol); // Tests += operator
      BOOST_CHECK_SMALL(s3half(x) -
                            static_cast<T>(3) * s(x) / static_cast<T>(2),
                        tol); // Tests - operator
      BOOST_CHECK_SMALL(s3half2(x) -
                            static_cast<T>(3) * s(x) / static_cast<T>(2),
                        tol); // Tests -= operator
      BOOST_CHECK_SMALL(splusone(x) - s(x) - static_cast<T>(1),
                        tol);                     // Tests + operator
      BOOST_CHECK_SMALL(sdx2(x) - sdx22(x), tol); // Tests dx method
      BOOST_CHECK_SMALL(s(x) - sdx0(x), tol);     // Tests dx method
    }
  }
}

BOOST_AUTO_TEST_CASE(TestArithmetic) {
  constexpr double TOL = 1.0e-15;
  testArithmetic<double, 2>(TOL);
  testArithmetic<double, 3>(TOL);
  testArithmetic<double, 4>(TOL);
  testArithmetic<double, 5>(TOL);
  testArithmetic<double, 6>(TOL);
  testArithmetic<double, 7>(TOL);
  testArithmetic<double, 8>(TOL);
  testArithmetic<double, 9>(TOL);
  testArithmetic<double, 10>(TOL);

  if constexpr (sizeof(long double) != sizeof(double)) {
    constexpr long double TOLL = 1.0e-18l;
    testArithmetic<long double, 2>(TOLL);
    testArithmetic<long double, 3>(TOLL);
    testArithmetic<long double, 4>(TOLL);
    testArithmetic<long double, 5>(TOLL);
    testArithmetic<long double, 6>(TOLL);
    testArithmetic<long double, 7>(TOLL);
    testArithmetic<long double, 8>(TOLL);
    testArithmetic<long double, 9>(TOLL);
    testArithmetic<long double, 10>(TOLL);
  }
}
