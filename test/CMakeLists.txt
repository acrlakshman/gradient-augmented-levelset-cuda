### Below code was used to verify GNU compiler
### version on travis-ci servers
IF (CMAKE_COMPILER_IS_GNUCC)
    EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -dumpversion
        OUTPUT_VARIABLE GCC_VERSION)
    STRING(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS
        ${GCC_VERSION})
    LIST(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
    LIST(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)

    MESSAGE(STATUS "Using GNUCC Compiler Verson ${GCC_MAJOR}.${GCC_MINOR}")
ENDIF()

IF (NOT TARGET gtest)
  # gtest targets cmake version 2.6 which throws warnings for policy CMP0042 on
  # newer cmakes. This sets the default global setting for that policy.
  SET(CMAKE_POLICY_DEFAULT_CMP0042 NEW)
  IF(WIN32)
    SET(gtest_force_shared_crt ON
      CACHE INTERNAL "Required so that the libs Runtime is not set to MT DLL")
  ENDIF()

  SET(GOOGLETEST_VERSION 1.8.1) #DUMMY to suppress cmake bug
  ADD_SUBDIRECTORY(googletest/googletest EXCLUDE_FROM_ALL)
  SET_TARGET_PROPERTIES(gtest gtest_main
    PROPERTIES
      FOLDER "ExternalProjectTargets/gtest")

  # Hide gtest project variables
  MARK_AS_ADVANCED(
    BUILD_SHARED_LIBS
    gtest_build_samples
    gtest_build_tests
    gtest_disable_pthreads
    gtest_force_shared_crt
    gtest_hide_internal_symbols)
ENDIF()

FILE(GLOB FILES "*.cpp" "*.cc")

SET(TEST_NAME gals_unit_test)
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR} TEST_FILE)
ADD_EXECUTABLE(${TEST_NAME} ${TEST_FILE})

SET_PROPERTY(TARGET ${TEST_NAME} PROPERTY CXX_STANDARD 14)

TARGET_INCLUDE_DIRECTORIES(${TEST_NAME}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/src
)

TARGET_LINK_LIBRARIES(${TEST_NAME}
  PUBLIC
    cpu_test
  PRIVATE
    gtest
    gtest_main
)

IF (BUILD_COVERAGE)
    TARGET_COMPILE_OPTIONS(${TEST_NAME}
      PRIVATE
        -g -O0 -fprofile-arcs -ftest-coverage
      )
    TARGET_LINK_LIBRARIES(${TEST_NAME}
      PRIVATE
        -g -O0 -fprofile-arcs -ftest-coverage
      )
ENDIF (BUILD_COVERAGE)

SET_TARGET_PROPERTIES(${TEST_NAME} PROPERTIES FOLDER "Tests")

ADD_TEST(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
