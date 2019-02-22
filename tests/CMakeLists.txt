add_executable(
  c2InitTests
  catch2Main.cpp catch2InitTests.cpp
  )

target_link_libraries(c2InitTests Catch2::Catch2 lambrex)
catch_discover_tests(c2InitTests)
