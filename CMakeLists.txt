CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

PROJECT(gals_cuda)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

OPTION(BUILD_DOCUMENTATION "Build Documentation" OFF)
OPTION(BUILD_EXAMPLES "Build examples" OFF)
OPTION(BUILD_TESTS "Build Tests" ON)
OPTION(BUILD_COVERAGE "Create test coverage report" OFF)
OPTION(USE_SYSTEM_GTEST "Use GTEST from system libraries" OFF)
OPTION(USE_RELATIVE_TEST_DIR "Use relative paths for the test data directory(For continious integration(CI) purposes only)" OFF)

# TODO (lakshman): Move googletest to ext/
IF ((NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/googletest"))
  MESSAGE(WARNING "third party libraries are missing.")
  MESSAGE("Did you miss the --recursive option when cloning?")
  MESSAGE("Run the following commands to correct this:")
  MESSAGE("git submodule init")
  MESSAGE("git submodule update")
  MESSAGE("git submodule foreach git pull origin master")
ENDIF ()

# Set a default build type if none was specified
IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
      "MinSizeRel" "RelWithDebInfo")
ENDIF()

MARK_AS_ADVANCED(BUILD_COVERAGE)

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

IF (BUILD_TESTS)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILD_TESTS")
ENDIF()

IF (BUILD_COVERAGE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage -DBUILD_COVERAGE")
ENDIF()

# Generate documentation
IF (BUILD_DOCUMENTATION)
    ADD_SUBDIRECTORY(docs)
ENDIF(BUILD_DOCUMENTATION)

IF (BUILD_EXAMPLES)
    ADD_SUBDIRECTORY(examples)
ENDIF()

IF (BUILD_TESTS)
    ENABLE_TESTING()
    ADD_SUBDIRECTORY(tests)
ENDIF()

# Build applications
ADD_SUBDIRECTORY(applications)

FILE(GLOB_RECURSE GALS_CPU_SRC
  src/cpu/*.cc
  src/cpu/*.h
  )

# Add external dependencies
ADD_SUBDIRECTORY(ext)

INCLUDE_DIRECTORIES(${YAML_CPP_INCLUDE_DIR} ${SPDLOG_INCLUDE_DIR} ${FILESYSTEM_INCLUDE_DIR})

SET(GALS_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/include")
INCLUDE_DIRECTORIES(${GALS_INCLUDE})

ADD_LIBRARY(gals_cpu SHARED ${GALS_CPU_SRC})

IF (APPLE)
  TARGET_LINK_LIBRARIES(gals_cpu PRIVATE yaml-cpp spdlog filesystem)
ENDIF()

ADD_EXECUTABLE(gals_cuda src/main.cc)
TARGET_INCLUDE_DIRECTORIES(gals_cuda PUBLIC src/cpu)
TARGET_LINK_LIBRARIES(gals_cuda PUBLIC gals_cpu PRIVATE yaml-cpp spdlog filesystem)

INSTALL(TARGETS gals_cuda DESTINATION .)
