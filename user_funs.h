#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1T(matrix, matrix = NAN, matrix = NAN);
long long fib_num(int k);
matrix ff1S(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1C(matrix, matrix = NAN, matrix = NAN);
matrix gram_schmidt(const matrix& Q_star, int n);
matrix ff2T(matrix x, matrix ud1, matrix ud2);
matrix ff3T(matrix x, matrix ud1, matrix ud2);
matrix ff4R(matrix x, matrix ud1, matrix ud2);
matrix df4(double t, matrix Y, matrix ud1, matrix ud2);






