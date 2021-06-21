cmake_minimum_required(VERSION 3.10.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# setting build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# setting compiler flags
set(GCC_OPENMP_COMPILE_FLAGS "-fopenmp -D use_openmp -Ofast -march=native -ffast-math")
set(GCC_OPENMP_LINK_FLAGS "-fopenmp ")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_OPENMP_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} ${GCC_OPENMP_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_OPENMP_LINK_FLAGS}")


project(IBSLib VERSION 0.2.0)

# find python
find_program(
    PYTHON_EXECUTABLE
    NAMES python
)
message("Python Executable ${PYTHON_EXECUTABLE}")
# if python found 
# get the PYTHON_LIBRARY_DIR to allow to install files
if(PYTHON_EXECUTABLE)
    message("Python used: ${PYTHON_EXECUTABLE}")
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" -c
        "from distutils import sysconfig;print(sysconfig.get_config_var('LIBDIR'))"
            OUTPUT_VARIABLE PYTHON_LIBRARY_DIR
            RESULT_VARIABLE PYTHON_LIBRARY_DIR_RESULT
            ERROR_QUIET)    
    string(STRIP ${PYTHON_LIBRARY_DIR} PYTHON_LIBRARY_DIR)
    if(NOT PYTHON_LIBRARY_DIR_RESULT MATCHES 0)
        message(SEND_ERROR "Failed to determine PYTHON_LIBRARY")
    endif()
else()
    message(FATAL_ERROR "python executable not found.")
endif()

# adding pybind11
# git submodule add https://github.com/pybind/pybind11.git third_party/pybind11-2.5.0
SET(MY_PYBIND ${CMAKE_CURRENT_SOURCE_DIR}/third_party/pybind11-2.5.0)
add_subdirectory(${MY_PYBIND})


include_directories("${CMAKE_SOURCE_DIR}/cpp/include/ibs_bits")
include_directories("${CMAKE_SOURCE_DIR}/python")

file (GLOB SOURCE_FILES "cpp/src/*.cpp")
file (GLOB HEADER_FILES "cpp/include/ibs_bits/*.hpp")
file (GLOB PYTHON_FILES "python/*.cpp" "python/*.hpp")
message("Source files ${SOURCE_FILES}")
message("Header files ${HEADER_FILES}")
message("Python files ${PYTHON_FILES}")
# Set up such that XCode organizes the files
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES} ${PYTHON_FILES} )

# include(pybind11.cmake)
pybind11_add_module(${PROJECT_NAME} 
    ${SOURCE_FILES}
    ${HEADER_FILES}
    ${PYTHON_FILES}
)


target_link_libraries(${PROJECT_NAME}  PUBLIC)

# show the PYTHON_LIBRARY_DIR
message("PYTHON_LIBRARY_DIR : ${PYTHON_LIBRARY_DIR}")
install(TARGETS ${PROJECT_NAME} 
    COMPONENT python
    LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
)