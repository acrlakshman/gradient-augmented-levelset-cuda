INCLUDE(CheckCXXCompilerFlag)

IF(NOT IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/yaml-cpp")
  MESSAGE(FATAL_ERROR "Missing dependencies. Please clone the project with --recursive."
    "It is possible to recover by calling \"git submodule update --init --recursive\"")
ENDIF()

IF(NOT IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/spdlog")
  MESSAGE(FATAL_ERROR "Missing dependencies. Please clone the project with --recursive."
    "It is possible to recover by calling \"git submodule update --init --recursive\"")
ENDIF()

# yaml-cpp
SET(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Turn off YAML-CPP tools target" FORCE)
SET(YAML_CPP_INSTALL     OFF CACHE BOOL "Turn off YAML-CPP install target" FORCE)
SET(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Turn off YAML-CPP tests target" FORCE)

ADD_SUBDIRECTORY(yaml-cpp)

SET_PROPERTY(TARGET yaml-cpp PROPERTY FOLDER "dependencies")

SET(YAML_CPP_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/yaml-cpp/include
  )

# spdlog
ADD_SUBDIRECTORY(spdlog)

SET_PROPERTY(TARGET spdlog PROPERTY FOLDER "dependencies")

SET(SPDLOG_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/spdlog/include)

SET(CompilerFlags
  CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
  CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
  YAML_CPP_INCLUDE_DIR
  SPDLOG_INCLUDE_DIR
  )

FOREACH(CompilerFlag ${CompilerFlags})
  SET(${CompilerFlag} "${${CompilerFlag}}" PARENT_SCOPE)
ENDFOREACH()
