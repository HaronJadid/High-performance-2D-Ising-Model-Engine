cmake_minimum_required(VERSION 3.20)
project(IsingModelSim)

set(CMAKE_CXX_STANDARD 17)

# 1. Find OpenMP Package
find_package(OpenMP REQUIRED)

# 2. Define Executable (ONLY ONCE)
add_executable(IsingModelSim src/main.cpp src/Lattice.cpp src/Lattice.h)

# 3. Compiler Options (Warnings + Optimization)
if(MSVC)
    target_compile_options(IsingModelSim PRIVATE /W4 /O2)
else()
    target_compile_options(IsingModelSim PRIVATE -Wall -Wextra -O3)
endif()

# 4. Link OpenMP if found
if(OpenMP_CXX_FOUND)
    target_link_libraries(IsingModelSim PUBLIC OpenMP::OpenMP_CXX)
endif()