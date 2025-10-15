cmake_minimum_required(VERSION 4.0)
project(optymalizacja)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(optymalizacja
        main.cpp
        matrix.cpp
        matrix.h
        ode_solver.cpp
        ode_solver.h
        opt_alg.cpp
        opt_alg.h
        solution.cpp
        solution.h
        user_funs.cpp
        user_funs.h)
