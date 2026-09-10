list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT TARGET "test-drive::test-drive")
  find_package("test-drive" REQUIRED)
endif()


# Unit testing
set(
  tests
  "random"
  "math"
  "common"
  "MCMC"
  "APMCMC"
  "DEMCz"
  "model"
  #"wrappers_for_R"
)
set(
  test-srcs
  "test_functions.f90"
  "main.f90"
)
foreach(t IN LISTS tests)
  string(MAKE_C_IDENTIFIER ${t} t) 
  list(APPEND test-srcs "test_${t}.f90")
endforeach()

add_executable(
  "${PROJECT_NAME}-tester"
  "${test-srcs}"
)

target_link_libraries(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "cardamom-lib"
  "test-drive::test-drive"
)

foreach(t IN LISTS tests)
  add_test("${PROJECT_NAME}/${t}" "${PROJECT_NAME}-tester" "${t}")
endforeach()
