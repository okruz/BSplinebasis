#define BOOST_TEST_MODULE SupportTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <spline-internal.h>


template<typename T>
void testSupport() {
   using support = myspline::internal::support<T>;
   using grid = std::shared_ptr<const std::vector<T>>;

   grid grid1 = std::make_shared<const std::vector<T>>(std::vector<T>{-7.0l,-6.85l,  -6.55l, -6.3l, -6.0l, -5.75l, -5.53l, -5.2l, -4.75l, -4.5l, -3.0l, -2.5l, -1.5l, -1.0l, 0.0l, 0.5l, 1.5l, 2.5l, 3.5l, 4.0l, 4.35l, 4.55l, 4.95l, 5.4l, 5.7l, 6.1l, 6.35l, 6.5l, 6.85l, 7.0l});

   // Represents the same grid as grid1
   grid grid11 = std::make_shared<const std::vector<T>>(std::vector<T>{-7.0l,-6.85l,  -6.55l, -6.3l, -6.0l, -5.75l, -5.53l, -5.2l, -4.75l, -4.5l, -3.0l, -2.5l, -1.5l, -1.0l, 0.0l, 0.5l, 1.5l, 2.5l, 3.5l, 4.0l, 4.35l, 4.55l, 4.95l, 5.4l, 5.7l, 6.1l, 6.35l, 6.5l, 6.85l, 7.0l});

   // Represents a different grid than grid1
   grid grid2 = std::make_shared<const std::vector<T>> (std::vector<T>{-7.0l,-6.85l,  -6.55l, -6.3l, -6.0l, -5.75l, -5.53l, -5.2l, -4.75l, -4.5l, -3.0l, -2.5l, -1.5l, -1.0l, 0.0l, 0.5l, 1.53l, 2.5l, 3.5l, 4.0l, 4.35l, 4.55l, 4.95l, 5.4l, 5.7l, 6.1l, 6.35l, 6.5l, 6.85l, 7.0l});

   // Represents a different grid than grid1
   grid grid3 = std::make_shared<const std::vector<T>>(std::vector<T>{-7.0l,-6.85l,  -6.55l, -6.3l, -6.0l, -5.75l, -5.53l, -5.2l, -4.75l, -4.5l, -3.0l, -2.5l, -1.5l, -1.0l, 0.0l, 0.5l, 1.5l, 2.5l, 3.5l, 4.0l, 4.35l, 4.55l, 4.95l, 5.4l, 5.7l, 6.1l, 6.35l, 6.5l, 6.85l, 7.0l, 8.0l});
    
   
   support s1(grid1, 0, grid1->size());
   support s11(grid11, 0, grid11->size());
   support s2(grid1);
   support s3(grid1, 3, 5);
   support s32(grid1, 0, 2);
   support s3i = s3.calcIntersection(s32);
   support s3u = s3.calcUnion(s32);
   
   support s4(grid2, 0, grid2->size());
   support s5(grid3);
  
   BOOST_TEST(s1.hasSameGrid(s11));
   BOOST_TEST((s1 == s1));
   BOOST_TEST((s1 == s11));
   BOOST_TEST((s1.calcUnion(s11) == s1));
   BOOST_TEST((s1.calcIntersection(s11) == s1));
   BOOST_TEST((s1.calcIntersection(s2) == s2));
   BOOST_TEST((s1.calcUnion(s2) == s1));
   BOOST_TEST((s1.calcIntersection(s3) == s3));
   BOOST_TEST((s1.calcUnion(s3) == s1));
   BOOST_TEST((s3i.empty() && s3i.hasSameGrid(s3)));
   BOOST_TEST((s3u.getStartIndex() == 0 && s3u.getEndIndex() == 5 && s3u.hasSameGrid(s3)));
   BOOST_TEST(s1.hasSameGrid(s11));
   BOOST_TEST(s1.hasSameGrid(s2));
   BOOST_TEST(s1.hasSameGrid(s3));
   BOOST_TEST(s11.hasSameGrid(s3i));
   BOOST_TEST(s11.hasSameGrid(s3u));
   BOOST_TEST(!s1.hasSameGrid(s4));
   BOOST_TEST(!s1.hasSameGrid(s5));
}


BOOST_AUTO_TEST_CASE (TestSupport)
{

    testSupport<float>();
    testSupport<double>();
    testSupport<long double>();
}
