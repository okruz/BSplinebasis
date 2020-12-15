#define BOOST_TEST_MODULE SplineTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <spline-template.h>


template<typename T>
std::vector<myspline::myspline<T, 3>> getSplines(const std::vector<T> &grid) {
   std::vector<myspline::myspline<T, 3>> splines;
   for(size_t i = 0; i + 4 < grid.size(); i++) splines.push_back(myspline::generateBspline<T, 4>(grid, i));
   return splines;
}

template<typename T>
myspline::myspline<T,0 > getOne(const std::vector<T> &grid) {
   T onet = static_cast<T>(1);
   std::vector<std::array<T, 1>> coeffs(grid.size()-1, {onet});
   return myspline::myspline<T, 0>(grid, std::move(coeffs));
}



template<typename T>
void testIntegration(T tol) {
   using spline = myspline::myspline<T, 3>;
   using spline0 = myspline::myspline<T, 0>;
   const std::vector<T> grid{-3.0l, -2.5l, -1.5l, -1.0l, 0.0l, 0.5l, 1.5l, 2.5l, 3.5l, 4.0l};
   const std::vector<spline> splines = getSplines(grid);
   const spline0 one = getOne(grid);
   for (const auto &s1:splines) {
       for (const auto &s2: splines) {
           auto s2dx2 = s2.dx2();
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1, s2)- myspline::integrate<T>(s1*s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1, s2)- myspline::overlap<T>(one, s1*s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1, s2.timesx())- myspline::integrate_x<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1.timesx(), s2)- myspline::integrate_x<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1, s2.timesx().timesx())- myspline::integrate_x2<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1.timesx(), s2.timesx())- myspline::integrate_x2<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1.timesx().timesx(), s2)- myspline::integrate_x2<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1, s2.dx())- myspline::integrate_dx<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1.timesx(), s2.dx())- myspline::integrate_x_dx<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1, s2.dx().dx())- myspline::integrate_dx2<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1, s2dx2)- myspline::integrate_dx2<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1.timesx(), s2.dx().dx())- myspline::integrate_x_dx2<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1, s2.dx().dx().timesx())- myspline::integrate_x_dx2<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1.timesx().timesx(), s2.dx().dx())- myspline::integrate_x2_dx2<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1, s2.dx().dx().timesx().timesx())- myspline::integrate_x2_dx2<T>(s1, s2), tol);
           BOOST_CHECK_SMALL(myspline::overlap<T>(s1, s2dx2.timesx().timesx())- myspline::integrate_x2_dx2<T>(s1, s2), tol);
       }
       spline s1_inv_inv = s1.invert().invert();
       BOOST_CHECK_SMALL(myspline::overlap<T>(s1 - s1_inv_inv, s1 - s1_inv_inv), tol);
   }

}


BOOST_AUTO_TEST_CASE (TestIntegration)
{
    testIntegration<double>(1.0e-15);
    testIntegration<long double>(1.0e-18l);
}

template<typename T>
void testArithmetic(T tol) {
   using spline = myspline::myspline<T, 3>;
   using spline6 = myspline::myspline<T, 6>;
   using spline0 = myspline::myspline<T, 0>;
   using spline1 = myspline::myspline<T, 1>;
   const std::vector<T> grid{-3.0l, -2.5l, -1.5l, -1.0l, 0.0l, 0.5l, 1.5l, 2.5l, 3.5l, 4.0l};
   const std::vector<spline> splines = getSplines(grid);
   const spline0 one = getOne(grid);
   for (const auto &s: splines) {
       spline s2 = s * static_cast<T>(2);
       spline6 sprod = s * s;
       spline s22 = s;
       s22 *= static_cast<T>(2);
       spline shalf = s / static_cast<T>(2);
       spline shalf2 = s;
       shalf2 /= static_cast<T>(2);
       spline s5half = s2 + shalf;
       spline s5half2 = s2;
       s5half2 += shalf;
       spline s3half = s2 - shalf;
       spline s3half2 = s2;
       s3half2 -= shalf;
       spline splusone = s + one;
       spline sinverted = s.invert();
       spline1 sdx2 = s.dx().dx();
       spline1 sdx22 = s.dx2();
       spline sdx0 = s.template dx<0>();

       for(T x = -3.0L; x <= 4.0L; x+= 0.01L) {
           BOOST_CHECK_SMALL(sprod(x) - s(x) * s(x), tol); // Tests * operator
           BOOST_CHECK_SMALL(s2(x) - static_cast<T>(2) * s(x), tol); // Tests * operator
           BOOST_CHECK_SMALL(s22(x) - static_cast<T>(2) * s(x), tol); // Tests *= operator
           BOOST_CHECK_SMALL(shalf(x) - s(x)/static_cast<T>(2), tol); // Tests / operator
           BOOST_CHECK_SMALL(shalf2(x) - s(x)/static_cast<T>(2), tol); // Tests /= operator
           BOOST_CHECK_SMALL(s5half(x) - static_cast<T>(5)*s(x)/static_cast<T>(2), tol); // Tests + operator
           BOOST_CHECK_SMALL(s5half2(x) - static_cast<T>(5)*s(x)/static_cast<T>(2), tol); // Tests += operator
           BOOST_CHECK_SMALL(s3half(x) - static_cast<T>(3)*s(x)/static_cast<T>(2), tol); // Tests - operator
           BOOST_CHECK_SMALL(s3half2(x) - static_cast<T>(3)*s(x)/static_cast<T>(2), tol); // Tests -= operator
           BOOST_CHECK_SMALL(splusone(x) - s(x) - static_cast<T>(1), tol); // Tests + operator
           BOOST_CHECK_SMALL(sinverted(x) - s(-x), tol); // Tests invert method
           BOOST_CHECK_SMALL(sdx2(x) - sdx22(x), tol); // Tests dx method
           BOOST_CHECK_SMALL(s(x) - sdx0(x), tol); // Tests dx method
         
       }
   }
}


BOOST_AUTO_TEST_CASE (TestArithmetic)
{
    testArithmetic<double>(1.0e-15);
    testArithmetic<long double>(1.0e-18l);
}
