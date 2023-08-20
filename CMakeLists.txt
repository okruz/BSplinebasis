cmake_minimum_required(VERSION 3.15)
project(BSplinebasis)

set(CMAKE_CXX_STANDARD 17)  
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

if (CMAKE_BUILD_TYPE STREQUAL "Debug" AND BSPLINE_COVERAGE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  add_custom_target(coverage-report
      COMMAND gcovr --html-details report.html -r ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/tests
      DEPENDS run-test
  )
endif(CMAKE_BUILD_TYPE STREQUAL "Debug" AND BSPLINE_COVERAGE)


add_library(main_library INTERFACE)
target_include_directories(main_library INTERFACE
   include
)

add_subdirectory(examples)
add_subdirectory(tests)

if(BSPLINE_GENERATE_README)
  add_subdirectory(readme)
endif(BSPLINE_GENERATE_README)
