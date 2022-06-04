/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <bspline/support/Support.h>

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

template <typename T>
void testSupport() {
  using Support = bspline::support::Support<T>;
  using Construction = bspline::support::Construction;
  using BSplineException = bspline::exceptions::BSplineException;
  using Grid = bspline::support::Grid<T>;
  const T tol = static_cast<T>(1.0e-15l);

  Grid grid1(std::vector<T>{-7.0l,  -6.85l, -6.55l, -6.3l, -6.0l, -5.75l,
                            -5.53l, -5.2l,  -4.75l, -4.5l, -3.0l, -2.5l,
                            -1.5l,  -1.0l,  0.0l,   0.5l,  1.5l,  2.5l,
                            3.5l,   4.0l,   4.35l,  4.55l, 4.95l, 5.4l,
                            5.7l,   6.1l,   6.35l,  6.5l,  6.85l, 7.0l});

  // Represents the same grid as grid1
  Grid grid11(std::vector<T>{-7.0l,  -6.85l, -6.55l, -6.3l, -6.0l, -5.75l,
                             -5.53l, -5.2l,  -4.75l, -4.5l, -3.0l, -2.5l,
                             -1.5l,  -1.0l,  0.0l,   0.5l,  1.5l,  2.5l,
                             3.5l,   4.0l,   4.35l,  4.55l, 4.95l, 5.4l,
                             5.7l,   6.1l,   6.35l,  6.5l,  6.85l, 7.0l});

  // Represents a different grid than grid1
  Grid grid2(std::vector<T>{-7.0l,  -6.85l, -6.55l, -6.3l, -6.0l, -5.75l,
                            -5.53l, -5.2l,  -4.75l, -4.5l, -3.0l, -2.5l,
                            -1.5l,  -1.0l,  0.0l,   0.5l,  1.53l, 2.5l,
                            3.5l,   4.0l,   4.35l,  4.55l, 4.95l, 5.4l,
                            5.7l,   6.1l,   6.35l,  6.5l,  6.85l, 7.0l});

  // Represents a different grid than grid1
  Grid grid3(std::vector<T>{-7.0l, -6.85l, -6.55l, -6.3l, -6.0l, -5.75l, -5.53l,
                            -5.2l, -4.75l, -4.5l,  -3.0l, -2.5l, -1.5l,  -1.0l,
                            0.0l,  0.5l,   1.5l,   2.5l,  3.5l,  4.0l,   4.35l,
                            4.55l, 4.95l,  5.4l,   5.7l,  6.1l,  6.35l,  6.5l,
                            6.85l, 7.0l,   8.0l});

  for (const Grid &grid : {grid1, grid11, grid2, grid3}) {
    BOOST_REQUIRE_NO_THROW(grid.at(grid.size() - 1));
    BOOST_REQUIRE_THROW(grid.at(grid.size()), BSplineException);
  }

  Support s1(grid1, 0, grid1.size());
  Support s11(grid11, 0, grid11.size());
  Support s12(grid1, Construction::WHOLE_GRID);
  Support s2(grid1);
  Support s3(grid1, 3, 5);
  Support s32(grid1, 0, 2);
  Support spl(grid1, 1, 2);  // point-like support
  Support s3i = s3.calcIntersection(s32);
  Support s3u = s3.calcUnion(s32);

  Support s4(grid2, 0, grid2.size());
  Support s5(grid3);

  for (const Support &support : {s1, s11, s12, s2, s3, s32, spl, s3i, s3u}) {
    if (!support.empty()) {
      BOOST_REQUIRE_NO_THROW(support.at(support.size() - 1));
      BOOST_REQUIRE_THROW(support.at(support.size()), BSplineException);
    }
  }

  BOOST_CHECK_SMALL(s1.front() - grid1.front(), tol);
  BOOST_CHECK_SMALL(s1.back() - grid1.back(), tol);
  BOOST_CHECK_SMALL(spl.front() - grid1.at(1), tol);
  BOOST_CHECK_SMALL(spl.back() - grid1.at(1), tol);
  BOOST_TEST(s1.hasSameGrid(s11));
  BOOST_TEST((s1 == s1));
  BOOST_TEST((s1 == s11));
  BOOST_TEST((s1 == s12));
  BOOST_TEST((s1.calcUnion(s11) == s1));
  BOOST_TEST((s1.calcUnion(s12) == s1));
  BOOST_TEST((s2.calcUnion(s12) == s1));
  BOOST_TEST((s1.calcIntersection(s11) == s1));
  BOOST_TEST((s1.calcIntersection(s2) == s2));
  BOOST_TEST((s1.calcUnion(s2) == s1));
  BOOST_TEST((s1.calcIntersection(s3) == s3));
  BOOST_TEST((s1.calcUnion(s3) == s1));
  BOOST_TEST((s1.calcIntersection(spl) == spl));
  BOOST_TEST((s1.calcUnion(spl) == s1));
  BOOST_TEST((s3i.empty() && s3i.hasSameGrid(s3)));
  BOOST_TEST((s3u.getStartIndex() == 0 && s3u.getEndIndex() == 5 &&
              s3u.hasSameGrid(s3)));
  BOOST_TEST(s1.hasSameGrid(s11));
  BOOST_TEST(s1.hasSameGrid(s2));
  BOOST_TEST(s1.hasSameGrid(s3));
  BOOST_TEST(s11.hasSameGrid(s3i));
  BOOST_TEST(s11.hasSameGrid(s3u));
  BOOST_TEST(!s1.hasSameGrid(s4));
  BOOST_TEST(!s1.hasSameGrid(s5));
  BOOST_TEST(s3.relativeFromAbsolute(3).value() == 0);
  BOOST_TEST(s3.relativeFromAbsolute(4).value() == 1);
  BOOST_TEST(!s3.relativeFromAbsolute(2));
  BOOST_TEST(!s3.relativeFromAbsolute(8));
  BOOST_TEST(s3.absoluteFromRelative(0) == 3);
  BOOST_TEST(s3.absoluteFromRelative(1) == 4);
}


BOOST_AUTO_TEST_SUITE(SupportTestSuite)
BOOST_AUTO_TEST_CASE(TestSupport) {
  testSupport<float>();
  testSupport<double>();
  testSupport<long double>();
}
BOOST_AUTO_TEST_SUITE_END()
