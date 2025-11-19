#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix gram_schmidt(const matrix& Q_star, int n);
matrix ff2T(matrix x, matrix ud1, matrix ud2);
matrix ff2R(matrix x, matrix ud1, matrix ud2);
matrix df2R(double t, matrix Y, matrix ud1, matrix ud2);