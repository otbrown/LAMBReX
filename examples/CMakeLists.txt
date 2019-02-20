add_executable(
  c2InitTests
  catch2Main.cpp catch2InitTests.cpp ${CMAKE_SOURCE_DIR}/include/lambrex.h
)
target_compile_features(c2InitTests PUBLIC cxx_std_11)
target_compile_definitions(c2InitTests PRIVATE AMREX_SPACEDIM=${DIM})
target_include_directories(
  c2InitTests
  PRIVATE
    ${CMAKE_SOURCE_DIR}/libamrex/include
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)
target_link_directories(
  c2InitTests
  PRIVATE ${CMAKE_SOURCE_DIR}/liblambrex
)
target_link_libraries(c2InitTests PRIVATE lambrex AMReX::amrex gfortran Catch2::Catch2)
set_target_properties(
  c2InitTests
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bins
)
catch_discover_tests(c2InitTests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bins)

add_executable(pulse pulse.cpp ${CMAKE_SOURCE_DIR}/include/lambrex.h)
target_compile_features(pulse PUBLIC cxx_std_11)
target_compile_options(
  pulse
  PRIVATE
    # clang/GCC warnings
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall>
    # MSVC warnings
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)
target_compile_definitions(pulse PRIVATE AMREX_SPACEDIM=${DIM})
target_include_directories(
  pulse
  PRIVATE
    ${CMAKE_SOURCE_DIR}/libamrex/include
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)
target_link_directories(
  pulse
  PRIVATE ${CMAKE_SOURCE_DIR}/liblambrex
)
target_link_libraries(pulse PRIVATE lambrex AMReX::amrex gfortran)
set_target_properties(
  pulse
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bins
)

add_custom_command(
  TARGET pulse
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E echo "Running $<TARGET_FILE_NAME:pulse> ..."
  COMMAND pulse > pulse.out
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bins
)
