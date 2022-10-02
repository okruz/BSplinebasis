/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <bspline/support/Support.h>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <type_traits>

static_assert(!std::is_move_constructible_v<bspline::support::Grid<double>> &&
                  !std::is_move_assignable_v<bspline::support::Grid<double>>,
              "Grid is movable.");
static_assert(
    std::is_nothrow_move_constructible_v<bspline::support::Support<double>> &&
        std::is_nothrow_move_assignable_v<bspline::support::Support<double>>,
    "Support is not nothrow move constructible.");

void testCopyGrid() {
  // A Grid cannot be moved, only copied.
  using Grid = bspline::support::Grid<double>;

  Grid grid1(std::vector<double>{-7.0l,  -6.85l, -6.55l, -6.3l, -6.0l, -5.75l,
                                 -5.53l, -5.2l,  -4.75l, -4.5l, -3.0l, -2.5l,
                                 -1.5l,  -1.0l,  0.0l,   0.5l,  1.5l,  2.5l,
                                 3.5l,   4.0l,   4.35l,  4.55l, 4.95l, 5.4l,
                                 5.7l,   6.1l,   6.35l,  6.5l,  6.85l, 7.0l});
  auto grid2 = grid1;
  BOOST_TEST((grid1.getData() != nullptr && grid2.getData() != nullptr));
  BOOST_TEST((grid1 == grid2));
}

void testMoveSupport() {
  using Grid = bspline::support::Grid<double>;
  using Support = bspline::support::Support<double>;

  Grid grid1(std::vector<double>{-7.0l,  -6.85l, -6.55l, -6.3l, -6.0l, -5.75l,
                                 -5.53l, -5.2l,  -4.75l, -4.5l, -3.0l, -2.5l,
                                 -1.5l,  -1.0l,  0.0l,   0.5l,  1.5l,  2.5l,
                                 3.5l,   4.0l,   4.35l,  4.55l, 4.95l, 5.4l,
                                 5.7l,   6.1l,   6.35l,  6.5l,  6.85l, 7.0l});
  auto support1 = Support::createWholeGrid(grid1);
  auto support2 = std::move(support1);

  BOOST_TEST(support1.hasSameGrid(support2));
  BOOST_TEST(support1.empty());
}

void testIteration() {
  using Grid = bspline::support::Grid<double>;
  using Support = bspline::support::Support<double>;

  const std::vector<double> gridData{
      -7.0l, -6.85l, -6.55l, -6.3l, -6.0l, -5.75l, -5.53l, -5.2l, -4.75l, -4.5l,
      -3.0l, -2.5l,  -1.5l,  -1.0l, 0.0l,  0.5l,   1.5l,   2.5l,  3.5l,   4.0l,
      4.35l, 4.55l,  4.95l,  5.4l,  5.7l,  6.1l,   6.35l,  6.5l,  6.85l,  7.0l};

  Grid grid1(gridData);

  std::vector<double> gridData2;
  for (const auto &val : grid1) {
    gridData2.push_back(val);
  }
  BOOST_TEST((gridData == gridData2));
  std::vector<double> gridData3{grid1.begin(), grid1.end()};
  BOOST_TEST((gridData == gridData3));

  const Support support1{grid1, 5, 9};
  const std::vector<double> supportData{gridData.begin() + 5,
                                        gridData.begin() + 9};
  std::vector<double> supportData2;
  for (const auto &val : support1) {
    supportData2.push_back(val);
  }
  BOOST_TEST((supportData == supportData2));
  std::vector<double> supportData3{support1.begin(), support1.end()};
  BOOST_TEST((supportData == supportData3));
}

template <typename T>
void testSupport() {
  using Support = bspline::support::Support<T>;
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

  const Support s1(grid1, 0, grid1.size());
  const Support s11(grid11, 0, grid11.size());
  const auto s12 = Support::createWholeGrid(grid1);
  const auto s2 = Support::createEmpty(grid1);
  const Support s3(grid1, 3, 5);
  const Support s32(grid1, 0, 2);
  const Support spl(grid1, 1, 2);  // point-like support
  const Support s3i = s3.calcIntersection(s32);
  const Support s3u = s3.calcUnion(s32);

  const Support s4(grid2, 0, grid2.size());
  const auto s5 = Support::createEmpty(grid3);

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

BOOST_AUTO_TEST_CASE(TestCopyGrid) { testCopyGrid(); }

BOOST_AUTO_TEST_CASE(TestMoveSupport) { testMoveSupport(); }

BOOST_AUTO_TEST_CASE(TestIteration) { testIteration(); }

BOOST_AUTO_TEST_SUITE_END()
