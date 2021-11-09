# see https://internalpointers.com/post/modern-cmake-beginner-introduction
cmake_minimum_required(VERSION 3.15.4)
project(cfd LANGUAGES CXX)

find_package(LAPACK)

function(add_proj name)
    add_executable(${name} ${ARGN})
    target_compile_features(${name} PRIVATE cxx_std_11)
    target_compile_options(${name} PRIVATE -g -Wall -Wextra -Wpedantic -O3)
    target_link_libraries(${name} ${LAPACK_LIBRARIES})
endfunction()

# add targets
add_proj(heat_conduction  
         heat_conduction/main.cpp 
         heat_conduction/solution.cpp 
         heat_conduction/grid.cpp 
         heat_conduction/matrix.cpp
         heat_conduction/heat_equation.cpp)
