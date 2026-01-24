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
matrix ff5T(matrix x, matrix ud1, matrix ud2);
matrix gf5T(matrix x, matrix ud1, matrix ud2);
matrix Hf5T(matrix x, matrix ud1, matrix ud2);
matrix ff5T_1D(matrix h, matrix ud1, matrix ud2);
matrix ff_logistic(matrix theta, matrix X_data, matrix Y_data);
matrix gf_logistic(matrix theta, matrix X_data, matrix Y_data);
double classification_accuracy(const matrix& theta, const matrix& X_data, const matrix& Y_data);
void set_lab5_params(double a, double w);
matrix ff5R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff5R_1D(matrix h, matrix ud1, matrix ud2);
double calculate_f1(matrix x, double a);
double calculate_f2(matrix x, double a);
double get_f1(matrix x, double a);
double get_f2(matrix x, double a);
matrix ff5_1D(matrix h, matrix ud1, matrix ud2);
matrix ff5(matrix x, matrix ud1, matrix ud2);
matrix ff6T(matrix x, matrix ud1, matrix ud2);
matrix df6R(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff6R(matrix x, matrix ud1, matrix ud2);


