cmake_minimum_required(VERSION 3.10.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(IBSLib VERSION 0.2.0)

# include dir
include_directories(/usr/local/include)

# Src
AUX_SOURCE_DIRECTORY(src SRC_FILES)

# Headers
set(PROJECT_SOURCE_DIR "src")
set(PROJECT_INCLUDE_DIR "include/ibs_bits")

# Source files
set(SOURCE_FILES
    ${PROJECT_INCLUDE_DIR}/twiss.hpp
    ${PROJECT_SOURCE_DIR}/twiss.cpp
)


# Set up such that XCode organizes the files correctly
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES})

# Add Library
add_library(${PROJECT_NAME}  SHARED ${SOURCE_FILES})

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE include/)

# Install
install(TARGETS ${PROJECT_NAME} DESTINATION lib)

# Install the headers
install(FILES include/ibs DESTINATION include)

# create base directory
install(DIRECTORY include/ DESTINATION include)