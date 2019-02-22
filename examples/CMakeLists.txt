
add_executable(pulse pulse.cpp)

# clang/GCC warnings
target_compile_options(pulse
  PRIVATE $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall>
  )
target_link_libraries(pulse PRIVATE lambrex)
