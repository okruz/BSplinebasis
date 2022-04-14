#ifndef BSPLINE_EXAMPLES_HARMONICOSCILLATOR_SPLINE_H
#define BSPLINE_EXAMPLES_HARMONICOSCILLATOR_SPLINE_H
/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <bspline/Spline.h>

#include <vector>

#include "misc.h"

namespace bspline::examples::harmonic_oscillator_spline {
/**
 * Solves the quantum mechanical harmonic oscillator problem with a BSpline
 * basis and with the potential expressed as a Bspline.
 * @return Returns the (real) eigenenergies and the corresponding wavefunctions.
 */
std::vector<Eigenspace> solveHarmonicOscillator();

}  // namespace bspline::examples::harmonic_oscillator_spline
#endif  // BSPLINE_EXAMPLES_HARMONICOSCILLATOR_SPLINE_H
