add_library(lambrex STATIC Simulation.cpp ${CMAKE_SOURCE_DIR}/include/Simulation.h)
target_compile_features(lambrex PUBLIC cxx_std_11)
target_compile_options(lambrex
                        PRIVATE
                          # clang/GCC warnings
                            $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall>
                          # MSVC warnings
                            $<$<CXX_COMPILER_ID:MSVC>:/W4>
                       )
target_compile_definitions(lambrex PRIVATE AMREX_SPACEDIM=${DIM})
target_include_directories(lambrex
                            PUBLIC
                              ${CMAKE_SOURCE_DIR}/libamrex/include
                              ${CMAKE_SOURCE_DIR}/include
                           )
set_target_properties(lambrex PROPERTIES
                      ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/liblambrex
                      )
