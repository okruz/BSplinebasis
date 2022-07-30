cmake_minimum_required(VERSION 3.15)
project(BSplinebasis)

set(CMAKE_CXX_STANDARD 17)  
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")


add_library(main_library INTERFACE)
target_include_directories(main_library INTERFACE
   include
)

add_subdirectory(examples)
add_subdirectory(tests)

if(BSPLINE_GENERATE_README)
  add_subdirectory(readme)
endif(BSPLINE_GENERATE_README)
