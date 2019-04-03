add_library(lambrex STATIC AmrSim.cpp)

target_compile_features(lambrex PUBLIC cxx_std_11)

# clang/GCC warnings
target_compile_options(lambrex
  PRIVATE $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall>
  )
target_include_directories(lambrex
  PUBLIC ${CMAKE_SOURCE_DIR}/include
  )
target_link_libraries(lambrex
  PUBLIC AMReX::amrex gfortran
  )
