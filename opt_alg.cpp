#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wej�ciowe:
	// ff - wska�nik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i g�rne ograniczenie
	// epslion - zak��dana dok�adno�� rozwi�zania
	// Nmax - maksymalna liczba wywo�a� funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy warto�� funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi�zanie z zadan� dok�adno�ci�
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo�a� funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double y0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution X0(x0);
		solution X1(x0 + d);

		double* p = new double[2] { 0, 0 };
		
		int i = 0;

		X0.fit_fun(ff);
		X1.fit_fun(ff);

		double f0 = m2d(X0.y);
		double f1 = m2d(X1.y);

		if (f1 == f0)
		{
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}
		if (f1 > f0)
		{
			d = -d;
			X1.x = x0 + d;
			X1.fit_fun(ff);
			f1 = m2d(X1.y);

			if (f1 >= f0)
			{
				p[0] = m2d(X1.x);
				p[1] = x0 - d;
				return p;
			}
		}

		solution X_prev(X1.x);
		double f_prev = f1;

		while (true)
		{
			if (solution::f_calls > Nmax)
			{
				delete[] p;
				throw std::string("expansion: function calls exceeded Nmax");
			}
			i++;

			solution X_next(X0.x + pow(alpha, i) * d);
			X_next.fit_fun(ff);
			double f_next = m2d(X_next.y);

			if (f_prev <= f_next)
			{
				if (d > 0)
				{
					p[0] = m2d(X_prev.x);
					p[1] = m2d(X_next.x);
				}
				else
				{
					p[0] = m2d(X_next.x);
					p[1] = m2d(X_prev.x);
				}
				return p;
			}
			X_prev = X_next;
			f_prev = f_next;
		}

		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution::clear_calls();
		int n = get_len(x0);
		solution xB(x0);
		xB.fit_fun(ff, ud1, ud2);

		solution x = xB;
		int exit_flag = -1;

		do
		{
			xB = x;
			x = HJ_trial(ff, xB, s, ud1, ud2);

			if (x.y < xB.y)
			{
				do
				{
					solution xB_prev = xB;
					xB = x;
					matrix x_new = 2.0 * xB.x - xB_prev.x; // 2*xB - xB_prev (pattern move)
					x = solution(x_new);
					x = HJ_trial(ff, x, s, ud1, ud2); // PRÓBUJ on the new point
					
					if (solution::f_calls > Nmax)
					{
						exit_flag = 0;
						break;
					}
				} while (x.y < xB.y);

				x = xB;
			}
			else
			{
				s = alpha * s;
			}

			if (solution::f_calls > Nmax)
			{
				exit_flag = 0;
				break;
			}
		} while (s >= epsilon);

		xB.flag = (exit_flag == 0) ? 0 : 1;
		return xB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_len(XB.x);
		solution X = XB;
		matrix e(n, 1, 0.0);

		for (int j = 0; j < n; ++j)
		{
			e = matrix(n, 1, 0.0);
			e(j) = 1.0;

			solution Xp(X.x + s * e);
			Xp.fit_fun(ff, ud1, ud2);

			if (Xp.y < X.y)
			{
				X = Xp;
			}
			else
			{
				solution Xm(X.x - s * e);
				Xm.fit_fun(ff, ud1, ud2);

				if (Xm.y < X.y)
					X = Xm;
			}
		}
		return X;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

matrix gram_schmidt(const matrix& Q_star, int n)
{
	matrix D_next(n, n);

	matrix v1 = get_col(Q_star, 0);

	double norm_v1 = norm(v1);
	if (std::abs(norm_v1) < 1e-12) {
		v1 = matrix(n, 1);
		v1(0) = 1.0;
		norm_v1 = norm(v1);
	}
	matrix d1 = v1 * (1.0 / norm_v1);
	D_next.set_col(d1, 0);

	for (int j = 1; j < n; ++j)
	{
		matrix vj = get_col(Q_star, j);

		for (int k = 0; k < j; ++k)
		{
			matrix dk = get_col(D_next, k);
			matrix product = trans(vj) * dk;
			double scalar = m2d(product);

			vj = vj - dk * scalar;
		}

		double norm_vj = norm(vj);
		if (std::abs(norm_vj) < 1e-12) {
			matrix ej = matrix(n, 1);
			ej((j + 1) % n) = 1.0;
			D_next.set_col(ej, j);

		}
		else {
			matrix dj = vj * (1.0 / norm_vj);
			D_next.set_col(dj, j);
		}
	}

	return D_next;
}


solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int n = get_len(x0);
		solution::clear_calls();

		int i = 0;

		matrix D = ident_mat(n); // D(i) przechowuje kierunki d_j(i) jako kolumny
		matrix Lambda(n, 1, 0.0);
		matrix P(n, 1, 0.0);
		solution xB(x0); // xB to wektor n x 1, przechowuje najlepszy punkt w bie��cej iteracji

		xB.fit_fun(ff, ud1, ud2);

		matrix S = s0; // Wektor d�ugo�ci krok�w s(i)

		int exit_flag = -1;

		do
		{
			double max_s_abs = 0.0;
			for (int j = 0; j < n; ++j) {
				max_s_abs = std::max(max_s_abs, std::abs(S(j)));
			}
			if (max_s_abs < epsilon) {
				exit_flag = 2; // Zbie�no��
				break;
			}

			solution x_prev = xB;
			int successful_steps = 0; // Licznik udanych krok�w w p�tli wewn�trznej

			for (int j = 0; j < n; ++j)
			{
				// Aktualny kierunek d_j(i) i d�ugo�� kroku s_j(i)
				matrix dj = get_col(D, j);
				double sj = S(j);

				solution x_new(xB.x + dj * sj);
				x_new.fit_fun(ff, ud1, ud2);

				if (x_new.y < xB.y)
				{
					xB = x_new;
					Lambda(j) += sj;
					S(j) *= alpha;
					successful_steps++; //???
				}
				else
				{
					S(j) *= (-beta);
					P(j) += 1.0;
				}
			}

			i++;
			bool re_orthogonalize = false;
			for (int j = 0; j < n; ++j) {
				if (std::abs(Lambda(j)) >= 1e-12 && P(j) >= 1.0) { // Pj musi by� wi�ksze lub r�wne 1
					re_orthogonalize = true;
					break;
				}
			}

			if (re_orthogonalize)
			{
				matrix L(n, n, 0.0);
				for (int j = 0; j < n; ++j) {
					for (int k = 0; k <= j; ++k) {
						L(j, k) = Lambda(k);
					}
				}

				// Q(i) = D(i) * L. D(i) to aktualna macierz D
				matrix Q = D * L;

				matrix Q_star = Q;

				D = gram_schmidt(Q_star, n);

				Lambda = matrix(n, 1, 0.0);
				P = matrix(n, 1, 0.0);
				S = s0;
			}

			if (solution::f_calls > Nmax)
			{
				exit_flag = 0; // Przekroczenie Nmax
				break;
			}

		} while (true);

		Xopt = xB;
		Xopt.flag = exit_flag;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
