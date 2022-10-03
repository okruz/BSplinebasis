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

static const std::vector<double> DEFAULT_GRID_DATA{
    -7.0l, -6.85l, -6.55l, -6.3l, -6.0l, -5.75l, -5.53l, -5.2l, -4.75l, -4.5l,
    -3.0l, -2.5l,  -1.5l,  -1.0l, 0.0l,  0.5l,   1.5l,   2.5l,  3.5l,   4.0l,
    4.35l, 4.55l,  4.95l,  5.4l,  5.7l,  6.1l,   6.35l,  6.5l,  6.85l,  7.0l};

static void testCopyGrid() {
  // A Grid cannot be moved, only copied.
  using Grid = bspline::support::Grid<double>;

  Grid grid1(DEFAULT_GRID_DATA);
  auto grid2 = grid1;
  BOOST_TEST((grid1.getData() != nullptr && grid2.getData() != nullptr));
  BOOST_TEST((grid1 == grid2));
}

static void testIteration() {
  using Grid = bspline::support::Grid<double>;

  Grid grid1(DEFAULT_GRID_DATA);

  std::vector<double> gridData2;
  for (const auto &val : grid1) {
    gridData2.push_back(val);
  }
  BOOST_TEST((DEFAULT_GRID_DATA == gridData2));
  std::vector<double> gridData3{grid1.begin(), grid1.end()};
  BOOST_TEST((DEFAULT_GRID_DATA == gridData3));
}

BOOST_AUTO_TEST_SUITE(GridTestSuite)

BOOST_AUTO_TEST_CASE(TestCopyGrid) { testCopyGrid(); }

BOOST_AUTO_TEST_CASE(TestIteration) { testIteration(); }

BOOST_AUTO_TEST_SUITE_END()
