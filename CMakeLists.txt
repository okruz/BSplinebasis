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
add_executable(splineInterpolationTest tests/okruz/bspline/interpolation/interpolation-test.cpp)
set_property(TARGET splineInterpolationTest PROPERTY CXX_STANDARD 17)
target_include_directories(splineInterpolationTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(splineInterpolationTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(splineInterpolationTest splineInterpolationTest)

# TODO: Check whether armadillo and eigen are installed and set MYSPLINE_INTERPOLATION_USE_* accordingly
target_compile_definitions(splineInterpolationTest PUBLIC 
                                          OKRUZ_BSPLINE_INTERPOLATION_USE_ARMADILLO
                                          OKRUZ_BSPLINE_INTERPOLATION_USE_EIGEN)
target_link_libraries(splineInterpolationTest armadillo)



add_executable(examples examples/main.cpp 
                        examples/harmonic-oscillator.cpp
                        examples/hydrogen.cpp)

set_property(TARGET examples PROPERTY CXX_STANDARD 17)
target_include_directories(examples PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(examples PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/examples)

# Might be needed, if you want to use boost multiprecision types as data_t in the examples.
# target_include_directories(examples PUBLIC /usr/include/eigen3)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Werror -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

enable_testing()
