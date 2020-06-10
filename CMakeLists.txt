# Minimum version of CMake required to build this project
cmake_minimum_required(VERSION 3.0)

# Name of the project
project(girih)

# Add compile options
add_compile_options(-g -Wall)

# Add all the source files needed to build the executable
add_executable(girih src/main.cpp 
                     src/utilities/affinity/affinity.cpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Link the executable and the library together
target_link_libraries(girih -lpthread)