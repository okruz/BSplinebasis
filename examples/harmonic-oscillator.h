#ifndef OKRUZ_BSPLINE_EXAMPLES_HARMONICOSCILLATOR_H
#define OKRUZ_BSPLINE_EXAMPLES_HARMONICOSCILLATOR_H
/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <okruz/bspline/Spline.h>

#include <vector>

#include "misc.h"

namespace okruz::bspline::examples::harmonic_oscillator {
/**
 * Solves the quantum mechanical harmonic oscillator problem with a BSpline
 * basis.
 * @return Returns the (real) eigenenergies and the corresponding wavefunctions.
 */
std::vector<Eigenspace> solveHarmonicOscillator();

}  // namespace okruz::bspline::examples::harmonic_oscillator
#endif  // OKRUZ_BSPLINE_EXAMPLES_HARMONICOSCILLATOR_H
