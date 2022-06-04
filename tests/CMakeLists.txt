find_package(Boost COMPONENTS unit_test_framework)
find_package(Armadillo)
find_package(Eigen3)

if (Boost_unit_test_framework_FOUND)

  #Add support test
  add_executable(supportTest bspline/support/support-test.cpp)
  set_property(TARGET supportTest PROPERTY CXX_STANDARD 17)
  target_include_directories(supportTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
  target_include_directories(supportTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(supportTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  add_test(supportTest supportTest)

  #Add spline test
  add_executable(splineTest bspline/spline-test.cpp)
  set_property(TARGET splineTest PROPERTY CXX_STANDARD 17)
  target_include_directories(splineTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
  target_include_directories(splineTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(splineTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  add_test(splineTest splineTest)

  if (ARMADILLO_FOUND OR Eigen3_FOUND)
    #Add spline-interpolation test
    add_executable(splineInterpolationTest bspline/interpolation/interpolation-test.cpp)
    set_property(TARGET splineInterpolationTest PROPERTY CXX_STANDARD 17)
    target_include_directories(splineInterpolationTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
    target_include_directories(splineInterpolationTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(splineInterpolationTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    add_test(splineInterpolationTest splineInterpolationTest)

    if (ARMADILLO_FOUND)                                     
      target_link_libraries(splineInterpolationTest armadillo)
      target_include_directories(splineInterpolationTest PUBLIC ${ARMADILLO_INCLUDE_DIRS})
      target_compile_definitions(splineInterpolationTest PUBLIC 
             BSPLINE_INTERPOLATION_USE_ARMADILLO
      )
    endif (ARMADILLO_FOUND)
  
    if (Eigen3_FOUND)
      target_link_libraries(splineInterpolationTest Eigen3::Eigen)
      target_compile_definitions(splineInterpolationTest PUBLIC 
             BSPLINE_INTERPOLATION_USE_EIGEN
      )
    endif (Eigen3_FOUND)
  endif (ARMADILLO_FOUND OR Eigen3_FOUND)
  
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Werror -Wno-deprecated-declarations")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
  enable_testing()
endif (Boost_unit_test_framework_FOUND)