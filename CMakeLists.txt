CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(gals_cuda)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

OPTION(BUILD_DOCUMENTATION "Build Documentation" OFF)
OPTION(BUILD_EXAMPLES "Build examples" OFF)
OPTION(BUILD_TESTS "Build Tests" ON)
OPTION(BUILD_COVERAGE "Create test coverage report" OFF)
OPTION(USE_SYSTEM_GTEST "Use GTEST from system libraries" OFF)
OPTION(USE_RELATIVE_TEST_DIR "Use relative paths for the test data directory(For continious integration(CI) purposes only)" OFF)

IF ((NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/googletest") OR
  (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/yaml-cpp"))

  MESSAGE(WARNING "third party libraries are missing.")
  MESSAGE("Did you miss the --recursive option when cloning?")
  MESSAGE("Run the following commands to correct this:")
  MESSAGE("git submodule init")
  MESSAGE("git submodule update")
  MESSAGE("git submodule foreach git pull origin master")

ENDIF ()

IF( NOT TARGET yaml-cpp )
  # Turn off YAML-CPP targets
  set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Turn off YAML-CPP tools target" FORCE)
  set(YAML_CPP_INSTALL     OFF CACHE BOOL "Turn off YAML-CPP install target" FORCE)
  set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Turn off YAML-CPP tests target" FORCE)

  ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/third_party/yaml-cpp)
  SET_TARGET_PROPERTIES(yaml-cpp
    PROPERTIES
      FOLDER "ExternalProjectTargets/yaml-cpp")
ENDIF()

# Set a default build type if none was specified
IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
      "MinSizeRel" "RelWithDebInfo")
ENDIF()

MARK_AS_ADVANCED(BUILD_COVERAGE)

IF(UNIX)
    IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    ELSEIF("${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS 4.7)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    ELSE()
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    ENDIF()
ENDIF()

IF (BUILD_COVERAGE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
ENDIF()

# Generate documentation
IF(BUILD_DOCUMENTATION)
    ADD_SUBDIRECTORY(doc)
ENDIF(BUILD_DOCUMENTATION)

IF(BUILD_EXAMPLES)
    ADD_SUBDIRECTORY(examples)
ENDIF()

IF(BUILD_TESTS)
    ENABLE_TESTING()
    ADD_SUBDIRECTORY(test)
ENDIF()

FILE(GLOB_RECURSE GALS_CPU_SRC src/cpu/*.cc src/cpu/*.h)

# Add third party dependencies
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/third_party/yaml-cpp/include/)

ADD_LIBRARY(gals_cpu SHARED ${GALS_CPU_SRC})
ADD_EXECUTABLE(gals_cuda src/main.cc)
TARGET_INCLUDE_DIRECTORIES(gals_cuda PUBLIC src/cpu)
TARGET_LINK_LIBRARIES(gals_cuda PUBLIC gals_cpu PRIVATE yaml-cpp)

INSTALL(TARGETS gals_cuda DESTINATION .)
