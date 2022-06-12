cmake_minimum_required(VERSION 3.15)
project(BSplinebasis)

add_library(main_library INTERFACE)
target_include_directories(main_library INTERFACE
   include
)

add_subdirectory(examples)
add_subdirectory(tests)

if(BSPLINE_GENERATE_README)
  add_subdirectory(readme)
endif(BSPLINE_GENERATE_README)
