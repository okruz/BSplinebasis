/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include "harmonic-oscillator.h"

#include <bspline/BSplineGenerator.h>

#include <algorithm>
#include <eigen3/Eigen/Eigenvalues>

namespace bspline::examples::harmonic_oscillator {

using namespace bspline;

// Hamiltonian operator -1/2 d^2/dx^2 + 1/2 x^2
static const auto hamiltonOperator =
    (static_cast<data_t>(1) / 2) * (-operators::Dx<2>{} + operators::X<2>{});

/**
 * @brief setUpKnotsVector Sets up the knots vector on which the basis splines
 * are defined.
 * @return A vector of data_t representing the knots.
 */
static std::vector<data_t> setUpKnotsVector() {
  std::vector<data_t> ret{0};

  for (int i = 1; i < 100; i++) {
    const data_t val = static_cast<data_t>(i * i) / 1000;
    ret.push_back(val);
    ret.push_back(-val);
  }
  std::sort(ret.begin(), ret.end());
  return ret;
}

/**
 * @brief setUpBasis Sets up the BSpline basis.
 * @return A vector of BSplines representing the basis.
 */
static std::vector<Spline> setUpBasis() {
  BSplineGenerator gen(setUpKnotsVector());
  return gen.template generateBSplines<SPLINE_ORDER + 1>();
}

std::vector<Eigenspace> solveHarmonicOscillator() {
  // Get the basis.
  std::vector<Spline> basis = setUpBasis();

  DeMat hamiltonian =
      setUpSymmetricMatrix(integration::BilinearForm{hamiltonOperator}, basis);

  // overlap matrix
  DeMat overlapMatrix =
      setUpSymmetricMatrix(integration::ScalarProduct{}, basis);

  // Solve the generalized eigenvalue problem A.x = lambda B.x
  Eigen::GeneralizedSelfAdjointEigenSolver<DeMat> ges{hamiltonian,
                                                      overlapMatrix};

  // Retrieve the (real) eigenvalues and eigenvectors.
  const auto &eigenvalues = ges.eigenvalues();
  const auto &eigenvectors = ges.eigenvectors();

  std::vector<Eigenspace> ret;
  ret.reserve(10);
  // Return the eienvalues and eigenfunctions correspondig to the ten lowest
  // eigenvalues.
  for (size_t i = 0; i < 10; i++) {
    const auto eigenvalue = eigenvalues(i);
    const auto eigenvector = toStdVector(eigenvectors.col(i));
    auto wavefunction = bspline::linearCombination(eigenvector, basis);

    ret.push_back({eigenvalue, std::move(wavefunction)});
  }
  return ret;
}

}  // namespace bspline::examples::harmonic_oscillator
