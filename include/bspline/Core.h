/*
 * ########################################################################
 * The contents of this file is free and unencumbered software released into the
 * public domain. For more information, please refer to <http://unlicense.org/>
 * ########################################################################
 */

#include <bspline/BSplineGenerator.h>
#include <bspline/Spline.h>
#include <bspline/integration/BilinearForm.h>
#include <bspline/integration/LinearForm.h>
#include <bspline/operators/CompoundOperators.h>
#include <bspline/operators/Derivative.h>
#include <bspline/operators/Position.h>
#include <bspline/operators/ScalarOperators.h>
#include <bspline/operators/SplineOperator.h>

/*!
 * @mainpage BSplineBasis: A BSpline library for use in basis expansions.
 *
 * @section overview Overview
 * This library contains a small template-based C++ BSpline-library geared
 * towards the use as basis functions for a basis-expansion ansatz to certain
 * analytical problems. It provides facilities for the generation, evaluation
 * and integration of the BSplines.
 *
 * The Spline class as well as the generation facilities can be found in the
 * top-level @link bspline @endlink namespace. The namespace
 * bspline::integration contains the code for the evaluation of integrals over
 * the Splines, most notably the BilinearForm and LinearForm. The operators
 * defined in the bspline::operators namespace can be applied to the Splines
 * either directly or during the evaluation of the integrals.
 */
