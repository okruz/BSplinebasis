// compile via: g++ -std=c++17 -O3 -I ../include/ -I /usr/include/eigen3/ -Wall -fsanitize=address -fsanitize=undefined -DMYSPLINE_INTERPOLATION_USE_ARMADILLO -DMYSPLINE_INTERPOLATION_USE_EIGEN -DARMA_DONT_USE_WRAPPER -o st spline-interpolation-test.cpp -lopenblas
#define BOOST_TEST_MODULE SplineTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <spline-interpolation.h>


#ifdef MYSPLINE_INTERPOLATION_USE_EIGEN
template<typename T, size_t order>
void testInterpolationEigen(T tol) {
   using spline = myspline::myspline<T, order>;
   const std::vector<T> x{-3.0l, -2.5l, -1.5l, -1.0l, 0.0l, 0.5l, 1.5l, 2.5l, 3.5l, 4.0l, 5.0l};
   const std::vector<T> y{-3.0l, -2.5l, -1.5l, -1.0l, 0.0l, -0.5l, -1.5l, -2.5l, -3.5l, -4.0l, 3.0l};
   spline s = myspline::interpolate_using_eigen<T, order>(x,y); 
   for (size_t i = 0; i < x.size(); i++) {
       BOOST_CHECK_SMALL(s(x[i]) - y[i], tol);
   }
}


BOOST_AUTO_TEST_CASE (TestInterpolationEigen)
{
    testInterpolationEigen<double, 1>(2.0e-14);
    testInterpolationEigen<long double, 1>(1.0e-17l);
    testInterpolationEigen<double, 2>(2.0e-14);
    testInterpolationEigen<long double, 2>(1.0e-17l);
    testInterpolationEigen<double, 3>(2.0e-14);
    testInterpolationEigen<long double, 3>(1.0e-17l);
    testInterpolationEigen<double, 4>(2.0e-14);
    testInterpolationEigen<long double, 4>(1.0e-17l);
}
#endif

#ifdef MYSPLINE_INTERPOLATION_USE_ARMADILLO
template<size_t order>
void testInterpolationArmadillo(double tol) {
   using spline = myspline::myspline<double, order>;
   const std::vector<double> x{-3.0l, -2.5l, -1.5l, -1.0l, 0.0l, 0.5l, 1.5l, 2.5l, 3.5l, 4.0l, 5.0l};
   const std::vector<double> y{-3.0l, -2.5l, -1.5l, -1.0l, 0.0l, -0.5l, -1.5l, -2.5l, -3.5l, -4.0l, 3.0l};
   spline s = myspline::interpolate_using_armadillo<order>(x,y); 
   for (size_t i = 0; i < x.size(); i++) {
       BOOST_CHECK_SMALL(s(x[i]) - y[i], tol);
   }
}


BOOST_AUTO_TEST_CASE (TestInterpolationArmadillo)
{
    testInterpolationArmadillo<1>(2.0e-14);
    testInterpolationArmadillo<2>(2.0e-14);
    testInterpolationArmadillo<3>(2.0e-14);
    testInterpolationArmadillo<4>(2.0e-14);
}
#endif
