find_package(Eigen3)

if (Eigen3_FOUND)
  add_executable(accuracy
    accuracy/harmonic-oscillator.cpp
  )

  set_property(TARGET accuracy PROPERTY CXX_STANDARD 17)

  target_link_libraries(accuracy PUBLIC main_library)
  target_include_directories(accuracy PUBLIC /usr/include/eigen3)

  target_compile_options(accuracy PRIVATE -pthread)

  add_custom_target(generate-accuracy-plots
    COMMAND
      ${CMAKE_CURRENT_BINARY_DIR}/accuracy ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND
      python3 ${CMAKE_CURRENT_SOURCE_DIR}/accuracy/plot.py ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/accuracy.png
    DEPENDS accuracy
  )

else()
  message(WARNING "Eigen3 framework not found. Not building accuracy test.")
endif (Eigen3_FOUND)


add_executable(generation
  generation/generation.cpp
)

set_property(TARGET generation PROPERTY CXX_STANDARD 17)

target_link_libraries(generation PUBLIC main_library)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Werror -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

add_custom_target(generation-plots
  COMMAND
    ${CMAKE_CURRENT_BINARY_DIR}/generation ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    python3 ${CMAKE_CURRENT_SOURCE_DIR}/generation/plot.py ${CMAKE_CURRENT_BINARY_DIR}/splines_normal.txt ${CMAKE_CURRENT_SOURCE_DIR}/splines_normal.png
  COMMAND
    python3 ${CMAKE_CURRENT_SOURCE_DIR}/generation/plot.py ${CMAKE_CURRENT_BINARY_DIR}/splines_non_continuous.txt ${CMAKE_CURRENT_SOURCE_DIR}/splines_non_contiuous.png 3
  DEPENDS generation
)

