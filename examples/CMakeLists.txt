find_package(Eigen3)
find_package(BLAS)

if (Eigen3_FOUND)
  add_executable(examples 
    main.cpp 
    harmonic-oscillator.cpp
    spline-potential.cpp
    hydrogen.cpp
  )

  set_property(TARGET examples PROPERTY CXX_STANDARD 17)
  target_include_directories(examples PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_link_libraries(examples Eigen3::Eigen)

  if (BLAS_FOUND)
    target_link_libraries(examples ${BLAS_LIBRARIES})
    target_compile_definitions(examples PUBLIC 
      EIGEN_USE_BLAS
    )
  endif(BLAS_FOUND)

  # Might be needed, if you want to use boost multiprecision types as data_t in the examples.
  # target_include_directories(examples PUBLIC /usr/include/eigen3)
  
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Werror")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -flto")
else()
  message(WARNING "Eigen3 framework not found. Not building examples.")
endif (Eigen3_FOUND)
