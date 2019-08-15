add_executable(
  c2InitTests
  catch2Main.cpp catch2InitTests.cpp
  )

target_link_libraries(c2InitTests Catch2::Catch2 lambrex)
ParseAndAddCatchTests(c2InitTests)

add_executable(
  c2RegressionTests
  catch2Main.cpp catch2RegressionTests.cpp
)

target_link_libraries(c2RegressionTests Catch2::Catch2 lambrex)
ParseAndAddCatchTests(c2RegressionTests)

add_executable(
  c2AMRTests
  catch2Main.cpp catch2AMRTests.cpp AmrTest.cpp
)

target_link_libraries(c2AMRTests Catch2::Catch2 lambrex)
ParseAndAddCatchTests(c2AMRTests)
