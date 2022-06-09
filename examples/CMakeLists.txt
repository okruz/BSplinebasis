find_package(Eigen3)

if (Eigen3_FOUND)
  add_library(examples_objects OBJECT
    harmonic-oscillator.cpp
    spline-potential.cpp
    hydrogen.cpp
  )

  target_compile_definitions(examples_objects PUBLIC 
      BSPLINE_INTERPOLATION_USE_EIGEN
  )

  target_link_libraries(examples_objects PUBLIC Eigen3::Eigen)
  set_property(TARGET examples_objects PROPERTY CXX_STANDARD 17)


  target_include_directories(examples_objects PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}
  )

  add_executable(examples 
    main.cpp 
  )

  set_property(TARGET examples PROPERTY CXX_STANDARD 17)

  target_link_libraries(examples PUBLIC examples_objects)

  # Might be needed, if you want to use boost multiprecision types as data_t in the examples.
  # target_include_directories(examples_objects PUBLIC /usr/include/eigen3)
  
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Werror")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
else()
  message(WARNING "Eigen3 framework not found. Not building examples.")
endif (Eigen3_FOUND)
