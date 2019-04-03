add_executable(amr_pulse amr_pulse.cpp)

# clang/GCC warnings
target_compile_options(amr_pulse
  PRIVATE $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall>
  )
target_link_libraries(amr_pulse PRIVATE lambrex)
