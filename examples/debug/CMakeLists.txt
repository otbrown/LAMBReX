add_executable(amr_debug amr_debug.cpp)

# clang/GCC warnings
target_compile_options(amr_debug
  PRIVATE $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall>
  )

target_include_directories(amr_debug
  PUBLIC ${CMAKE_SOURCE_DIR}/tests
  )

target_link_libraries(amr_debug PRIVATE lambrex)
