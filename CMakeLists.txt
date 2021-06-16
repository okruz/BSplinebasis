cmake_minimum_required(VERSION 3.15)
project(BSplinebasis)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

#Add support test
add_executable(supportTest tests/okruz/bspline/support/support-test.cpp)
set_property(TARGET supportTest PROPERTY CXX_STANDARD 17)
target_include_directories(supportTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(supportTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(supportTest supportTest)

#Add spline test
add_executable(splineTest tests/okruz/bspline/spline-test.cpp)
set_property(TARGET splineTest PROPERTY CXX_STANDARD 17)
target_include_directories(splineTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(splineTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(splineTest splineTest)

#Add spline-interpolation test
#add_executable(splineInterpolationTest tests/okruz/bspline/interpolation/interpolation-test.cpp)
#set_property(TARGET splineInterpolationTest PROPERTY CXX_STANDARD 17)
#target_include_directories(splineInterpolationTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_link_libraries(splineInterpolationTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
#add_test(splineInterpolationTest splineInterpolationTest)

# TODO: Check whether armadillo and eigen are installed and set MYSPLINE_INTERPOLATION_USE_* accordingly
#target_compile_definitions(splineInterpolationTest PUBLIC 
#                                          MYSPLINE_INTERPOLATION_USE_ARMADILLO
#                                          MYSPLINE_INTERPOLATION_USE_EIGEN)
#target_link_libraries(splineInterpolationTest armadillo)
#target_include_directories(splineInterpolationTest PUBLIC "/usr/include/eigen3/")


set(CMAKE_CXX_FLAGS "-O3 -Wall")

enable_testing()