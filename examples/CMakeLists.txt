add_executable(pulse pulse.cpp ${CMAKE_SOURCE_DIR}/include/lambrex.h)
target_compile_features(pulse PUBLIC cxx_std_11)
target_compile_options(pulse
                        PRIVATE
                          # clang/GCC warnings
                            $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall>
                          # MSVC warnings
                            $<$<CXX_COMPILER_ID:MSVC>:/W4>
                      )
target_compile_definitions(pulse PRIVATE AMREX_SPACEDIM=${DIM})
target_include_directories(pulse
                            PRIVATE
                              ${CMAKE_SOURCE_DIR}/libamrex/include
                              ${CMAKE_SOURCE_DIR}/include
                              ${CMAKE_SOURCE_DIR}/src
                          )
target_link_directories(pulse
                         PRIVATE
                          ${CMAKE_SOURCE_DIR}/libamrex/lib
                          ${CMAKE_SOURCE_DIR}/liblambrex
                       )
target_link_libraries(pulse PRIVATE lambrex amrex gfortran)
set_target_properties(pulse PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bins
                      )
