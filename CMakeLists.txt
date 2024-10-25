cmake_minimum_required(VERSION 3.9...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.12)
endif()

project(elm-pb-mri LANGUAGES CXX C)

# Extra CMake scripts in cmake/ subdirectory
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

option(ELMPB_UPDATE_GIT_SUBMODULE "Check submodules are up-to-date during build" ON)
# Adapted from https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
# Update submodules as needed
function(elmpb_update_submodules)
  if(NOT ELMPB_UPDATE_GIT_SUBMODULE)
    return()
  endif()
  find_package(Git QUIET)
  if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} -c submodule.recurse=false submodule update --init --recursive
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
endfunction()
elmpb_update_submodules()

# Get the Git revision
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC ELMPB_REVISION)
if(ELMPB_REVISION STREQUAL "GITDIR-NOTFOUND")
  set(ELMPB_REVISION "Unknown")
endif()
message(STATUS "Git revision: ${ELMPB_REVISION}")

# BOUT++ is a dependency
option(ELMPB_BUILD_BOUT "Build BOUT++ in external/BOUT-dev" ON)
if(ELMPB_BUILD_BOUT)
  set(BOUT_BUILD_EXAMPLES, OFF) # Don't create example makefiles
  add_subdirectory(external/BOUT-dev)
else()
  find_package(bout++ REQUIRED)
endif()

set(ELMPB_SOURCES
    src/elm_pb_mri.cxx
    )

# The main executable target
add_executable(elm-pb-mri
               ${ELMPB_SOURCES}
               ${CMAKE_CURRENT_BINARY_DIR}/include/revision.hxx)

target_link_libraries(elm-pb-mri PRIVATE bout++::bout++)

target_include_directories(elm-pb-mri PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )

# Build the file containing just the commit hash
# This will be rebuilt on every commit!
configure_file(
  "${PROJECT_SOURCE_DIR}/include/revision.hxx.in"
  "${PROJECT_BINARY_DIR}/include/revision.hxx")

# Once built, copy the data and test directories
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/examples $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/tests $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests)

# Tests
option(ELMPB_TESTS "Build the tests" ON)
if(ELMPB_TESTS)
  enable_testing()

  # Integrated tests
  function(elmpb_add_integrated_test TESTNAME)
    add_test(NAME ${TESTNAME}
           WORKING_DIRECTORY tests/${TESTNAME}
           COMMAND runtest)
  endfunction()

  # Example of how to add a test
  # elmpb_add_integrated_test(sometest)
endif()

# Generate the build config header

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/elmpb_build_config.hxx")
  # If we do in source builds, this is fine
  if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
    message(FATAL_ERROR "Generated elmpb_build_config.hxx header already exists; please remove '${CMAKE_CURRENT_SOURCE_DIR}/include/elmpb_build_config.hxx' before continuing")
  endif()
endif()

configure_file(include/elmpb_build_config.hxx.in include/elmpb_build_config.hxx)

# Installation
install(TARGETS elm-pb-mri)
