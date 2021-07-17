cmake_minimum_required(VERSION 3.10.2)

project(IBSLibDocs)
find_package(Doxygen REQUIRED)
find_program(SPHINX_EXECUTABLE
             NAMES sphinx-build
             DOC "Path to sphinx-build executable")
include(FindPackageHandleStandardArgs)
get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

#Handle standard arguments to find_package like REQUIRED and QUIET
find_package_handle_standard_args(Sphinx
                                "Failed to find sphinx-build executable"
                                SPHINX_EXECUTABLE)

# Find all the public headers
file(GLOB_RECURSE IBS_PUBLIC_HEADERS ${PARENT_DIR}/cpp/include/ibs_bits/*.hpp)
message(STATUS "IBS_PUBLIC_HEADERS = ${IBS_PUBLIC_HEADERS}")

#This will be the main output of our command
set(DOXYGEN_INDEX_FILE ${CMAKE_CURRENT_SOURCE_DIR}/html/index.html)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/xml)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD  ${PARENT_DIR}/docs/build)
message(STATUS "CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "SPHINX_BUILD = ${SPHINX_BUILD}")

add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${IBS_PUBLIC_HEADERS}
                   COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                   MAIN_DEPENDENCY Doxyfile
                   COMMENT "Generating docs")

add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

add_custom_target(Sphinx ALL
                  COMMAND
                python -msphinx -b html
                  -Dbreathe_projects.ibs=${DOXYGEN_OUTPUT_DIR}
                  ${SPHINX_SOURCE} ${SPHINX_BUILD}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  COMMENT "Generating documentation with Sphinx")