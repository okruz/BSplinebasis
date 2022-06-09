/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <harmonic-oscillator.h>
#include <hydrogen.h>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

using bspline::examples::data_t;

BOOST_AUTO_TEST_SUITE(ExampleTestSuite)
BOOST_AUTO_TEST_CASE(HarmonicOscillatorTest) {
  const auto eigenSpaces =
      bspline::examples::harmonic_oscillator::solveHarmonicOscillator();
  static constexpr data_t TOL = 1.0e-12;

  for (size_t i = 0; i < eigenSpaces.size(); i++) {
    const data_t eigenvalue = eigenSpaces.at(i).energy;
    const data_t expected = static_cast<data_t>(2 * i + 1) / 2;
    BOOST_CHECK_SMALL((eigenvalue - expected) / abs(expected), TOL);
  }
}

BOOST_AUTO_TEST_CASE(HydrogenTest) {
  const auto eigenSpaces = bspline::examples::hydrogen::solveRadialHydrogen();
  static constexpr data_t TOL = 5.0e-12;

  for (size_t i = 0; i < eigenSpaces.size(); i++) {
    const data_t eigenvalue = eigenSpaces.at(i).energy;

    const size_t n = i + 2;
    const data_t expected = static_cast<data_t>(-1) / (n * n);

    BOOST_CHECK_SMALL((eigenvalue - expected) / abs(expected), TOL);
  }
}

BOOST_AUTO_TEST_SUITE_END()
