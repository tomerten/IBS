cmake_minimum_required(VERSION 3.10.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


project(test_cpp)

include_directories(/usr/local/include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/../bin)

find_library(IBSLIB_LIB IBSLib HINTS /usr/local/lib/)

add_executable(test_cpp src/test_cpp.cpp)
add_executable(test_radiation_damping_cpp src/DemoRadDamping.cpp)
add_executable(test_coulomblog_functions_cpp src/DemoCoulombLog.cpp)
add_executable(test_integrators_cpp src/DemoIntegrators.cpp)
add_executable(test_ibs_models_cpp src/DemoIBS.cpp)
add_executable(test_ibs_ode_cpp src/DemoODE.cpp)


target_link_libraries(test_cpp PUBLIC ${IBSLIB_LIB})
target_link_libraries(test_radiation_damping_cpp PUBLIC ${IBSLIB_LIB})
target_link_libraries(test_coulomblog_functions_cpp PUBLIC ${IBSLIB_LIB})
target_link_libraries(test_integrators_cpp PUBLIC ${IBSLIB_LIB})
target_link_libraries(test_ibs_models_cpp PUBLIC ${IBSLIB_LIB})
target_link_libraries(test_ibs_ode_cpp PUBLIC ${IBSLIB_LIB})