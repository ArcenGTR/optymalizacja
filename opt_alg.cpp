#include"opt_alg.h"
#include"user_funs.h"
#include<algorithm>
#include<random>

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejœciowe:
	// ff - wskaŸnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zak³¹dana dok³adnoœæ rozwi¹zania
	// Nmax - maksymalna liczba wywo³añ funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj¹c rozk³ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi¹zanie do przedzia³u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy wartoœæ funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi¹zanie z zadan¹ dok³adnoœci¹
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo³añ funkcji celu
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

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };

		solution x1(x0);
		solution x2(x0 + d);

		x1.fit_fun(ff, ud1, ud2);
		x2.fit_fun(ff, ud1, ud2);

		if (solution::f_calls >= Nmax)
		{
			p[0] = m2d(x1.x);
			p[1] = m2d(x2.x);
			if (p[0] > p[1]) std::swap(p[0], p[1]);
			return p;
		}

		if (m2d(x1.y) <= m2d(x2.y))
		{
			d = -d;
			// Zamieniamy miejscami x1 i x2, aby x2 zawsze by³ "lepszym" punktem
			std::swap(x1, x2);
		}

		while (true)
		{
			solution x3(m2d(x2.x) + d);
			x3.fit_fun(ff, ud1, ud2);

			//Sprawdzenie, czy "przeskoczyliœmy" minimum
			if (m2d(x3.y) > m2d(x2.y))
			{
				// Znaleziono przedzia³. Minimum jest miêdzy x1 a x3.
				p[0] = m2d(x1.x);
				p[1] = m2d(x3.x);
				// Porz¹dkowanie wyniku, aby p[0] < p[1]
				if (p[0] > p[1]) std::swap(p[0], p[1]);
				break;
			}

			// Sprawdzenie warunku stopu
			if (solution::f_calls >= Nmax)
			{
				// Zwracamy najlepszy znaleziony dot¹d przedzia³
				p[0] = m2d(x1.x);
				p[1] = m2d(x2.x);
				if (p[0] > p[1]) std::swap(p[0], p[1]);
				break;
			}

			// Przejœcie do nastêpnej iteracji
			x1 = x2;
			x2 = x3;
			d *= alpha; // Zwiêkszenie kroku
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

       // Krok 1: ZnajdŸ najmniejsz¹ liczbê k spe³niaj¹c¹ nierównoœæ Fk > (b - a) / ?
       int k = 0;
       double length_ratio = (b - a) / epsilon;
       while (fib_num(k) <= length_ratio)
       {
           k++;
           // Dodatkowe sprawdzenie, aby unikn¹æ przepe³nienia long long dla zbyt du¿ego k
           if (k > 90) throw string("fib(): Za duze k. Zmien dlugosc kroku lub epsilon.");
       }

       // Oblicz wymagane liczby Fibonacciego
       double Fk_minus_1 = (double)fib_num(k - 1);
       double Fk = (double)fib_num(k);

       // Krok 2: a(0) = a, b(0) = b
       double a_i = a;
       double b_i = b;

       // Krok 3 & 4: c(0), d(0)
       double c_i = b_i - (Fk_minus_1 / Fk) * (b_i - a_i);
       double d_i = a_i + b_i - c_i;

       // Krok i=0, i=1
       solution Xc(c_i), Xd(d_i);
       Xc.fit_fun(ff, ud1, ud2);
       Xd.fit_fun(ff, ud1, ud2);

       // Glowny cykl
       // Krok 5: for i = 0 to k – 3 (wykonujemy k-2 iteracji)
       for (int i = 0; i <= k - 3; ++i)
       {
           // Sprawdzenie kryterium stopu (opcjonalnie, ale dobre praktycznie)
           if (solution::f_calls >= 1000000) // Dodatkowy Nmax, jeœli jest taka potrzeba
           {
               Xopt.x = (a_i + b_i) / 2.0;
               Xopt.fit_fun(ff, ud1, ud2);
               Xopt.flag = 0; // Oznaczenie, ¿e osi¹gniêto Nmax
               return Xopt;
           }

           // Krok 6: if f(c(i)) < f(d(i)) then
           if (m2d(Xc.y) < m2d(Xd.y))
           {
               // Krok 7-8: a(i+1) = a(i), b(i+1) = d(i)
               b_i = m2d(Xd.x); // b(i+1) = d(i)
               // Wykorzystujemy Xc jako nowy Xd, a c(i+1) obliczamy od nowa
               Xd = Xc;

               // Krok 13 (obliczenie nowego c):
               // Fk-i-2 i Fk-i-1 sa potrzebne do nowego c
               double Fk_i_minus_2 = (double)fib_num(k - i - 2);
               double Fk_i_minus_1 = (double)fib_num(k - i - 1);
               c_i = b_i - (Fk_i_minus_2 / Fk_i_minus_1) * (b_i - a_i);

               Xc = solution(c_i);
               Xc.fit_fun(ff, ud1, ud2);
           }
           else
           {
               // Krok 10-11: b(i+1) = b(i), a(i+1) = c(i)
               a_i = m2d(Xc.x); // a(i+1) = c(i)
               // Wykorzystujemy Xd jako nowy Xc, a d(i+1) obliczamy od nowa
               Xc = Xd;

               // Krok 14 (obliczenie nowego d):
               // Fk-i-2 i Fk-i-1 sa potrzebne do nowego d
               double Fk_i_minus_2 = (double)fib_num(k - i - 2);
               double Fk_i_minus_1 = (double)fib_num(k - i - 1);
               // Nowe d to a(i+1) + b(i+1) - c(i+1), gdzie c(i+1) jest obliczone dla nowej d³. przedzia³u
               // Po uproszczeniu, nowy d to b(i+1) - (Fk-i-2 / Fk-i-1) * (b(i+1) - a(i+1))
               d_i = a_i + b_i - c_i;

               Xd = solution(d_i);
               Xd.fit_fun(ff, ud1, ud2);
           }
       		//std::cout << i + 1 << " " << (b_i - a_i) << std::endl;
           // Na koncu kazdej iteracji: dlugosc przedzialu [a_i, b_i] jest wlasciwa
       }

       // Krok 16: return x* = c(i+1);
       // Po pêtli ostatni punkt c lub d jest w zasadzie rozwi¹zaniem.
       // Minimum znajduje siê w ostatnim przedziale [a_i, b_i].
       // Wg pseudokodu c(k-2) jest rozwi¹zaniem
       Xopt = Xc; // Zwracamy c(k-2) - czyli Xc po ostatniej iteracji

       // Finalna korekta (wg literatury, minimum to œrodek ostatniego przedzia³u)
       Xopt.x = (a_i + b_i) / 2.0;
       Xopt.fit_fun(ff, ud1, ud2);
       Xopt.flag = 1;
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

       // Krok 2: a(0) = a, b(0) = b.
       // Krok 1: Deklaracja kluczowych zmiennych PRZED pêtl¹, aby by³y dostêpne w ca³ym zakresie.
       double a_i = a;
       double b_i = b;
       // Za³o¿enie dla punktu wewnêtrznego c(0) (wg pseudokodu musi byæ podany)
       // Jeœli nie jest podany, u¿ywamy œrodka:
       double c_i = (a + b) / 2.0;

       // Deklaracja zmiennych d_i i d_i_prev przed pêtl¹
       double d_i;
       double d_i_prev = c_i; // Wymagane do kryterium |d(i) - d(i-1)| < gamma

       // Deklaracja obiektów solution przed pêtl¹
       solution Xa(a_i), Xb(b_i), Xc(c_i), Xd;

       // Wymagane obliczenia pocz¹tkowe dla a, b, c
       Xa.fit_fun(ff, ud1, ud2);
       Xb.fit_fun(ff, ud1, ud2);
       Xc.fit_fun(ff, ud1, ud2);


       // Krok 3: repeat (i = 0 do nieskonczonoœci, z kryteriami stopu)
       for (int i = 0; solution::f_calls < Nmax; ++i)
       {
           // Sprawdzenie kryterium stopu na pocz¹tku pêtli
           // Krok 39: until b(i) – a(i) < ? or |d(i) – d(i-1)| < ?
           if (i > 0 && ((b_i - a_i < epsilon) || (fabs(d_i - d_i_prev) < gamma)))
           {
               Xopt.flag = 1; // Optymalizacja zakoñczona sukcesem
               break;
           }

           // Wymagane s¹ aktualne wartoœci funkcji w punktach a, b, c
           double fa = m2d(Xa.y), fb = m2d(Xb.y), fc = m2d(Xc.y);

           // Krok 4 & 5: Obliczenie l i m
           double l = fa * (pow(b_i, 2) - pow(c_i, 2)) + fb * (pow(c_i, 2) - pow(a_i, 2)) + fc * (pow(a_i, 2) - pow(b_i, 2));
           double m = fa * (b_i - c_i) + fb * (c_i - a_i) + fc * (a_i - b_i);

           // Krok 6 & 7: if m <= 0 then return error
           if (m <= 0.0)
           {
               Xopt.flag = -1;
               Xopt.x = c_i;
               Xopt.fit_fun(ff, ud1, ud2);
               return Xopt;
           }

           // Krok 9: d(i) = 0,5 * l / m
           d_i_prev = d_i; // Zapis d(i-1) przed obliczeniem d(i)
           d_i = 0.5 * l / m;
           Xd.x = d_i;

           // Sprawdzenie, czy d(i) jest w przedziale [a,b]
           if (a_i < d_i && d_i < b_i)
           {
               Xd.fit_fun(ff, ud1, ud2);

               // Krok 10: if a(i) < d(i) < c(i) then
               if (d_i < c_i)
               {
                   // Krok 11-14: if f(d(i)) < f(c(i))
                   if (m2d(Xd.y) < m2d(Xc.y))
                   {
                       // a(i+1)=a(i), c(i+1)=d(i), b(i+1)=c(i)
                       b_i = c_i;
                       c_i = d_i;
                       Xb = Xc;
                       Xc = Xd;
                   }
                   else // Krok 15-18
                   {
                       // a(i+1)=d(i), c(i+1)=c(i), b(i+1)=b(i)
                       a_i = d_i;
                       Xa = Xd;
                   }
               }
               // Krok 21: else if c(i) < d(i) < b(i) then
               else // (d_i > c_i)
               {
                   // Krok 22-25: if f(d(i)) < f(c(i))
                   if (m2d(Xd.y) < m2d(Xc.y))
                   {
                       // a(i+1)=c(i), c(i+1)=d(i), b(i+1)=b(i)
                       a_i = c_i;
                       c_i = d_i;
                       Xa = Xc;
                       Xc = Xd;
                   }
                   else // Krok 26-29
                   {
                       // a(i+1)=a(i), c(i+1)=c(i), b(i+1)=d(i)
                       b_i = d_i;
                       Xb = Xd;
                   }
               }
           }
           // Krok 31-33: else return error (d(i) jest poza [a, b])
           else
           {
               Xopt.flag = -1;
               Xopt.x = c_i;
               Xopt.fit_fun(ff, ud1, ud2);
               return Xopt;
           }

           // Krok 36: Sprawdzenie Nmax
           if (solution::f_calls >= Nmax)
           {
               Xopt.flag = 0; // Oznaczenie osi¹gniêcia Nmax
               break;
           }
       		//std::cout << i + 1 << " " << (b_i - a_i) << std::endl;
       } // end repeat

       // Krok 40: return x* = d(i)
       Xopt.x = d_i;
       Xopt.fit_fun(ff, ud1, ud2);
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
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
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
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
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
       solution xB(x0); // xB to wektor n x 1, przechowuje najlepszy punkt w bie¿¹cej iteracji

       xB.fit_fun(ff, ud1, ud2);

       matrix S = s0; // Wektor d³ugoœci kroków s(i)

       int exit_flag = -1;

       do
       {
           double max_s_abs = 0.0;
           for (int j = 0; j < n; ++j) {
               max_s_abs = std::max(max_s_abs, std::abs(S(j)));
           }
           if (max_s_abs < epsilon) {
               exit_flag = 2; // Zbie¿noœæ
               break;
           }

           solution x_prev = xB;
           int successful_steps = 0; // Licznik udanych kroków w pêtli wewnêtrznej

           for (int j = 0; j < n; ++j)
           {
           		// Aktualny kierunek d_j(i) i d³ugoœæ kroku s_j(i)
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
               if (std::abs(Lambda(j)) >= 1e-12 && P(j) >= 1.0) { // Pj musi byæ wiêksze lub równe 1
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

       }
       while (true);

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
		solution x_curr(x0);
		double current_c = c;
		matrix current_ud1 = ud1;

		while (true)
		{
			// Update penalty coefficient in user data
			// We assume ud1(0) is reserved for 'c'
			current_ud1(0) = current_c;

			// Run Simplex
			// Parameters for NM: s=0.5 (initial step), alpha=1, beta=0.5, gamma=2, delta=0.5
			double s = 0.5;
			double alpha = 1.0, beta = 0.5, gamma = 2.0, delta = 0.5;

			// We run NM with the current penalized function
			solution x_next = sym_NM(ff, x_curr.x, s, alpha, beta, gamma, delta, epsilon, Nmax, current_ud1, ud2);

			// Check convergence of the penalty loop
			// If change in x is very small, we stop.
			double diff = norm(x_next.x - x_curr.x);

			x_curr = x_next;

			if (diff < epsilon || solution::f_calls > Nmax)
			{
				Xopt = x_curr;
				break;
			}

			// Increase penalty coefficient
			current_c *= dc;
		}

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
		int n = get_dim(solution(x0));

		// 1. Create initial simplex
		// p[0] = x0
		std::vector<solution> p;
		p.push_back(solution(x0));
		p[0].fit_fun(ff, ud1, ud2);

		// p[i] = p[0] + s * ei
		matrix I = ident_mat(n);
		for (int i = 0; i < n; ++i)
		{
			solution new_pt(x0 + get_col(I, i) * s);
			new_pt.fit_fun(ff, ud1, ud2);
			p.push_back(new_pt);
		}

		while (true)
		{
			// 6. Sort vertices by function value
			// p[0] is min (best), p[n] is max (worst)
			std::sort(p.begin(), p.end(), [](const solution& a, const solution& b) {
				return m2d(a.y) < m2d(b.y);
			});

			solution p_min = p[0];
			solution p_max = p[n];
			solution p_next_max = p[n - 1]; // Second worst

			// Check stop criterion: max || p_min - p_i || < epsilon
			bool stop = true;
			for (int i = 1; i <= n; ++i)
			{
				if (norm(p[i].x - p_min.x) >= epsilon)
				{
					stop = false;
					break;
				}
			}
			if (stop || solution::f_calls >= Nmax)
			{
				Xopt = p_min;
				Xopt.flag = (solution::f_calls >= Nmax) ? 0 : 1;
				break;
			}

			// 8. Centroid of all except max
			matrix sum_x(n, 1);
			for (int i = 0; i < n; ++i) // sum 0 to n-1
			{
				sum_x = sum_x + p[i].x;
			}
			solution p_bar(sum_x * (1.0 / n));

			// 9. Reflection
			solution p_ref(p_bar.x + (p_bar.x - p_max.x) * alpha);
			p_ref.fit_fun(ff, ud1, ud2);

			if (m2d(p_ref.y) < m2d(p_min.y))
			{
				// 11. Expansion
				solution p_exp(p_bar.x + (p_ref.x - p_bar.x) * gamma);
				p_exp.fit_fun(ff, ud1, ud2);

				if (m2d(p_exp.y) < m2d(p_ref.y))
					p[n] = p_exp; // Accept expansion
				else
					p[n] = p_ref; // Accept reflection
			}
			else
			{
				if (m2d(p_ref.y) < m2d(p_next_max.y)) // Better than second worst?
				{
					p[n] = p_ref; // Accept reflection
				}
				else
				{
					// Contraction
					solution p_con;
					bool accepted_con = false;

					if (m2d(p_ref.y) < m2d(p_max.y))
					{
						// Outside contraction (towards reflection)
						p_con = solution(p_bar.x + (p_ref.x - p_bar.x) * beta);
					}
					else
					{
						// Inside contraction (towards max)
						p_con = solution(p_bar.x + (p_max.x - p_bar.x) * beta);
					}
					p_con.fit_fun(ff, ud1, ud2);

					if (m2d(p_con.y) < std::min(m2d(p_ref.y), m2d(p_max.y)))
					{
						p[n] = p_con;
					}
					else
					{
						// Reduction
						for (int i = 1; i <= n; ++i)
						{
							p[i].x = p[0].x + (p[i].x - p[0].x) * delta;
							p[i].fit_fun(ff, ud1, ud2);
						}
					}
				}
			}
		}

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
		solution x(x0); // Current point
		solution x_prev(x0);

		// Initial evaluation
		x.fit_fun(ff, ud1, ud2);

		while(true)
		{
			// 3. Compute gradient
			matrix g = gf(x.x, ud1, ud2);
			solution::g_calls++;

			// Direction d = -gradient
			matrix d = -g;

			double h = h0;

			// 4. Determine step size h
			if (h0 == 0) // Variable step size using Golden Section
			{
				// Search range for h: [0, 1.0] (arbitrary, but sufficient for this problem)
				solution h_sol = golden(ff5T_1D, 0.0, 1.0, epsilon, Nmax, x.x, d);
				h = m2d(h_sol.x);
			}

			// 5. Update position
			x_prev = x;
			x.x = x.x + d * h;
			x.fit_fun(ff, ud1, ud2);

			// Check stops
			if (solution::f_calls > Nmax || solution::g_calls > Nmax)
			{
				Xopt = x;
				Xopt.flag = 0;
				break;
			}

			// Stop criterion: ||x_i - x_i-1|| < epsilon
			if (norm(x.x - x_prev.x) < epsilon)
			{
				Xopt = x;
				Xopt.flag = 1;
				break;
			}
		}
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
		solution x(x0);
		solution x_prev(x0);

		x.fit_fun(ff, ud1, ud2);

		matrix g = gf(x.x, ud1, ud2);
		solution::g_calls++;

		// Initial direction d(0) = -g(0)
		matrix d = -g;
		matrix g_prev = g;

		while(true)
		{
			// Line Search / Step size
			double h = h0;
			if (h0 == 0)
			{
				solution h_sol = golden(ff5T_1D, 0.0, 1.0, epsilon, Nmax, x.x, d);
				h = m2d(h_sol.x);
			}

			// Update x
			x_prev = x;
			x.x = x.x + d * h;
			x.fit_fun(ff, ud1, ud2);

			// Check Convergence
			if (solution::f_calls > Nmax || solution::g_calls > Nmax)
			{
				Xopt = x;
				Xopt.flag = 0;
				break;
			}
			if (norm(x.x - x_prev.x) < epsilon)
			{
				Xopt = x;
				Xopt.flag = 1;
				break;
			}

			// Compute new gradient
			g = gf(x.x, ud1, ud2);
			solution::g_calls++;

			// Calculate Beta (Fletcher-Reeves)
			double beta = pow(norm(g), 2) / pow(norm(g_prev), 2);

			// Update direction: d(i) = -g(i) + beta * d(i-1)
			d = -g + d * beta;

			g_prev = g;
		}
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
		solution x(x0);
		solution x_prev(x0);

		x.fit_fun(ff, ud1, ud2);

		while(true)
		{
			matrix g = gf(x.x, ud1, ud2);
			solution::g_calls++;

			matrix H = Hf(x.x, ud1, ud2);
			solution::H_calls++;

			// Direction: d = -H^-1 * g
			matrix d = -inv(H) * g;

			double h = h0;
			if (h0 == 0)
			{
				solution h_sol = golden(ff5T_1D, 0.0, 1.0, epsilon, Nmax, x.x, d);
				h = m2d(h_sol.x);
			}

			x_prev = x;
			x.x = x.x + d * h;
			x.fit_fun(ff, ud1, ud2);

			if (solution::f_calls > Nmax || solution::g_calls > Nmax) // or H_calls
			{
				Xopt = x;
				Xopt.flag = 0;
				break;
			}
			if (norm(x.x - x_prev.x) < epsilon)
			{
				Xopt = x;
				Xopt.flag = 1;
				break;
			}
		}
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
		double alpha = (sqrt(5.0) - 1.0) / 2.0;

		double a_i = a;
		double b_i = b;
		double c_i = b_i - alpha * (b_i - a_i);
		double d_i = a_i + alpha * (b_i - a_i);

		solution Xc{matrix(c_i)};
		Xc.fit_fun(ff, ud1, ud2);

		solution Xd{matrix(d_i)};
		Xd.fit_fun(ff, ud1, ud2);

		while(true)
		{
			if (solution::f_calls > Nmax)
			{
				Xopt = Xc;
				Xopt.flag = 0;
				break;
			}

			if (m2d(Xc.y) < m2d(Xd.y))
			{
				b_i = d_i;
				d_i = c_i;
				Xd = Xc;

				c_i = b_i - alpha * (b_i - a_i);
				Xc = solution(matrix(c_i)); // Tutaj przypisanie jest bezpieczne
				Xc.fit_fun(ff, ud1, ud2);
			}
			else
			{
				a_i = c_i;
				c_i = d_i;
				Xc = Xd;

				d_i = a_i + alpha * (b_i - a_i);
				Xd = solution(matrix(d_i)); // Tutaj przypisanie jest bezpieczne
				Xd.fit_fun(ff, ud1, ud2);
			}

			if ((b_i - a_i) < epsilon)
			{
				Xopt.x = (a_i + b_i) / 2.0;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 1;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}



solution SD_Logged(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2, ofstream& trajectory_file, const std::string& method_name, double h_val)
{
    try
    {
       solution Xopt;
       solution x(x0); // Current point
       solution x_prev(x0);
       int iter = 0; // Dodano licznik iteracji

       // Pocz¹tkowa ewaluacja (wartoœæ y dla x0 zosta³a ju¿ zapisana w lab6_trajectory)
       x.fit_fun(ff, ud1, ud2);
       // UWAGA: Nie logujemy tutaj, bo iteracja 0 jest ju¿ zalogowana w lab6_trajectory

       while(true)
       {
          iter++; // Licznik iteracji

          // 3. Obliczenie gradientu
          matrix g = gf(x.x, ud1, ud2);
          solution::g_calls++;

          // Kierunek d = -gradient
          matrix d = -g;

          double h = h0;

          // 4. Okreœlenie d³ugoœci kroku h
          if (h0 == 0) // Zmienny krok (Golden Section)
          {
             // Za³o¿enie: golden(ff5T_1D, ...) jest dostêpne
             solution h_sol = golden(ff5T_1D, 0.0, 1.0, epsilon, Nmax, x.x, d);
             h = m2d(h_sol.x);
          }

          // 5. Aktualizacja pozycji
          x_prev = x;
          x.x = x.x + d * h;
          x.fit_fun(ff, ud1, ud2); // W tym miejscu inkrementowane jest solution::f_calls

          // --- LOGOWANIE TEJ ITERACJI ---
          trajectory_file << method_name << ";" << h_val << ";" << iter << ";"
                          << x.x(0) << ";" << x.x(1) << ";"
                          << m2d(x.y) << ";"
                          << solution::f_calls << ";"
                          << solution::g_calls << ";"
                          << solution::H_calls << "\n";
          // -----------------------------

          // Kryteria stopu
          if (solution::f_calls > Nmax || solution::g_calls > Nmax)
          {
             Xopt = x;
             Xopt.flag = 0;
             break;
          }

          // Kryterium stopu: ||x_i - x_i-1|| < epsilon
          if (norm(x.x - x_prev.x) < epsilon)
          {
             Xopt = x;
             Xopt.flag = 1;
             break;
          }
       }
       return Xopt;
    }
    catch (string ex_info)
    {
       throw ("solution SD_Logged(...):\n" + ex_info);
    }
}

solution CG_Logged(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2, ofstream& trajectory_file, const std::string& method_name, double h_val)
{
    try
    {
       solution Xopt;
       solution x(x0);
       solution x_prev(x0);
       int iter = 0; // Dodano licznik iteracji

       x.fit_fun(ff, ud1, ud2);

       matrix g = gf(x.x, ud1, ud2);
       solution::g_calls++;

       // Pocz¹tkowy kierunek d(0) = -g(0)
       matrix d = -g;
       matrix g_prev = g;

       while(true)
       {
          iter++; // Licznik iteracji

          // Liniowe wyszukiwanie / D³ugoœæ kroku
          double h = h0;
          if (h0 == 0)
          {
             solution h_sol = golden(ff5T_1D, 0.0, 1.0, epsilon, Nmax, x.x, d);
             h = m2d(h_sol.x);
          }

          // Aktualizacja x
          x_prev = x;
          x.x = x.x + d * h;
          x.fit_fun(ff, ud1, ud2);

          // --- LOGOWANIE TEJ ITERACJI ---
          trajectory_file << method_name << ";" << h_val << ";" << iter << ";"
                          << x.x(0) << ";" << x.x(1) << ";"
                          << m2d(x.y) << ";"
                          << solution::f_calls << ";"
                          << solution::g_calls << ";"
                          << solution::H_calls << "\n";
          // -----------------------------

          // Kryteria konwergencji
          if (solution::f_calls > Nmax || solution::g_calls > Nmax)
          {
             Xopt = x;
             Xopt.flag = 0;
             break;
          }
          if (norm(x.x - x_prev.x) < epsilon)
          {
             Xopt = x;
             Xopt.flag = 1;
             break;
          }

          // Obliczenie nowego gradientu
          g = gf(x.x, ud1, ud2);
          solution::g_calls++;

          // Obliczenie Beta (Fletcher-Reeves)
          double beta = pow(norm(g), 2) / pow(norm(g_prev), 2);

          // Aktualizacja kierunku: d(i) = -g(i) + beta * d(i-1)
          d = -g + d * beta;

          g_prev = g;
       }
       return Xopt;
    }
    catch (string ex_info)
    {
       throw ("solution CG_Logged(...):\n" + ex_info);
    }
}

solution Newton_Logged(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
    matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2, ofstream& trajectory_file, const std::string& method_name, double h_val)
{
    try
    {
       solution Xopt;
       solution x(x0);
       solution x_prev(x0);
       int iter = 0; // Dodano licznik iteracji

       x.fit_fun(ff, ud1, ud2);

       while(true)
       {
       	iter++;

       	matrix g = gf(x.x, ud1, ud2);
       	solution::g_calls++;

       	matrix H = Hf(x.x, ud1, ud2);
       	solution::H_calls++;

       	// 1. Kierunek: d = -H^-1 * g
       	matrix d = -inv(H) * g;

       	// ----------------------------------------------------
       	// >>> OSTATECZNA POPRAWKA BLOKU REGULARYZACJI KIERUNKU <<<
       	// Poprawne u¿ycie funkcji trans() z pliku matrix.cpp
       	if (m2d(trans(g) * d) >= 0) { // U¿ycie trans(g)
       		// Kierunek d nie jest kierunkiem zstêpu (descent).
       		// Stosujemy kierunek Najszybszego Spadku.
       		d = -g;
       	}
       	// ----------------------------------------------------

       	double h = h0;
       	if (h0 == 0)
       	{
       		// 2. Liniowe wyszukiwanie (Z³oty Podzia³)
       		solution h_sol = golden(ff5T_1D, 0.0, 1.0, epsilon, Nmax, x.x, d);
       		h = m2d(h_sol.x);
       	}

       	x_prev = x;
       	x.x = x.x + d * h;
       	x.fit_fun(ff, ud1, ud2);


          // --- LOGOWANIE TEJ ITERACJI ---
          trajectory_file << method_name << ";" << h_val << ";" << iter << ";"
                          << x.x(0) << ";" << x.x(1) << ";"
                          << m2d(x.y) << ";"
                          << solution::f_calls << ";"
                          << solution::g_calls << ";"
                          << solution::H_calls << "\n";
          // -----------------------------

          if (solution::f_calls > Nmax || solution::g_calls > Nmax) // lub H_calls
          {
             Xopt = x;
             Xopt.flag = 0;
             break;
          }
          if (norm(x.x - x_prev.x) < epsilon)
          {
             Xopt = x;
             Xopt.flag = 1;
             break;
          }
       }
       return Xopt;
    }
    catch (string ex_info)
    {
       throw ("solution Newton_Logged(...):\n" + ex_info);
    }
}



solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
    {
        solution Xopt;
        solution x(x0);
        int n = get_dim(x);
        matrix D = ident_mat(n); // Step 2: Initialize directions as identity matrix
        solution::clear_calls(); // Reset counters

        // We use a specialized 1D function wrapper that expects:
        // ud1 = current point x, ud2 = direction d.
        // The actual function ff is effectively ff5R (via ff5R_1D).

        while (true) // Step 3: repeat
        {
            matrix p0 = x.x; // Step 4: p0 = x(i)

            // Step 5: for j = 1 to n
            for (int j = 0; j < n; ++j)
            {
                matrix dj = get_col(D, j);

                // Step 6: Determine h_j (Line Search)
                // We use expansion to find a bracket, then golden to find min
                double* range = expansion(ff5R_1D, 0.0, 1.0, 2.0, Nmax, x.x, dj);
                solution h_sol = golden(ff5R_1D, range[0], range[1], epsilon, Nmax, x.x, dj);
                double h = m2d(h_sol.x);
                delete[] range;

                // Step 7: p_j = p_{j-1} + h*d_j (Update x immediately)
                x.x = x.x + dj * h;
            }

            // Step 9: Check convergence ||pn - p0|| < epsilon
            if (norm(x.x - p0) < epsilon || solution::f_calls > Nmax)
            {
                Xopt = x;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = (solution::f_calls > Nmax) ? 0 : 1;
                break; // Step 10: return x*
            }

            // Step 12-14: Shift directions
            // d_j = d_{j+1} for j=1 to n-1
            for (int j = 0; j < n - 1; ++j)
            {
                matrix next_d = get_col(D, j + 1);
                D.set_col(next_d, j);
            }

            // Step 15: d_n = p_n - p_0
            matrix dn = x.x - p0;

            // Step 16: Determine h_{n+1} (Line Search along new direction)
            double* range = expansion(ff5R_1D, 0.0, 1.0, 2.0, Nmax, x.x, dn);
            solution h_sol = golden(ff5R_1D, range[0], range[1], epsilon, Nmax, x.x, dn);
            double h = m2d(h_sol.x);
            delete[] range;

            // Step 17: p_{n+1} = p_n + h * d_{n+1} -> x_new
            x.x = x.x + dn * h;

            // Update direction set: D col n-1 becomes dn
            D.set_col(dn, n - 1);

            // Step 20: until f_calls > Nmax
            if (solution::f_calls > Nmax)
            {
                Xopt = x;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 0;
                break;
            }
        }

        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution Powell(...):\n" + ex_info);
    }
}



#include <random>

// Globalny generator
static std::random_device rd;
static std::mt19937 gen(rd());

double randn() {
	std::normal_distribution<double> d(0, 1);
	return d(gen);
}

double randd() {
	std::uniform_real_distribution<double> d(0, 1);
	return d(gen);
}

// Zmienione nazwy, by unikn¹æ konfliktu "ambiguous call"
matrix random_matrix(int n) {
	matrix res(n, 1);
	std::uniform_real_distribution<double> d(0, 1);
	for (int i = 0; i < n; ++i) res(i, 0) = d(gen);
	return res;
}

matrix randn_matrix(int n) {
	matrix res(n, 1);
	std::normal_distribution<double> d(0, 1);
	for (int i = 0; i < n; ++i) res(i, 0) = d(gen);
	return res;
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try {
        solution::clear_calls();
        double alpha = pow(N, -0.5);
        double beta = pow(2 * N, -0.25);

        // 1. Inicjalizacja populacji
        solution* P = new solution[mi];
        for (int i = 0; i < mi; ++i) {
            P[i].x = matrix(N, 1);
            for (int k = 0; k < N; ++k)
                P[i].x(k) = lb(k) + (ub(k) - lb(k)) * randd();
            P[i].ud = sigma0;
            P[i].fit_fun(ff, ud1, ud2);
        }

        solution Xopt = P[0];

        while (solution::f_calls < Nmax) {
            // 2. Selekcja turniejowa (Ko³o ruletki)
            double* phi = new double[mi];
            double Phi_sum = 0;
            for (int j = 0; j < mi; ++j) {
                phi[j] = 1.0 / (m2d(P[j].y) + 1e-12);
                Phi_sum += phi[j];
            }
            double* q = new double[mi + 1];
            q[0] = 0;
            for (int j = 1; j <= mi; ++j) q[j] = q[j - 1] + phi[j - 1] / Phi_sum;

            double a = randn();
            solution* T = new solution[lambda];

            for (int j = 0; j < lambda; ++j) {
                // Wybór rodziców
                int idxA = 0, idxB = 0;
                double r1 = randd(), r2 = randd();
                for (int k = 1; k <= mi; ++k) {
                    if (r1 > q[k - 1] && r1 <= q[k]) idxA = k - 1;
                    if (r2 > q[k - 1] && r2 <= q[k]) idxB = k - 1;
                }

                // 3. Krzy¿owanie
                double r = randd();
                T[j].x = matrix(N, 1);
                for (int k = 0; k < N; ++k)
                    T[j].x(k) = P[idxA].x(k) * r + P[idxB].x(k) * (1.0 - r);

                double sigma_avg = m2d(P[idxA].ud) * r + m2d(P[idxB].ud) * (1.0 - r);

                // 4. Mutacja (z pilnowaniem ograniczeñ)
                double b_rand = randn();
                T[j].ud = sigma_avg * exp(alpha * a + beta * b_rand);
                double s = m2d(T[j].ud);

                for (int k = 0; k < N; ++k) {
                    T[j].x(k) = T[j].x(k) + s * randn();
                    // Bariery
                    if (T[j].x(k) < lb(k)) T[j].x(k) = lb(k);
                    if (T[j].x(k) > ub(k)) T[j].x(k) = ub(k);
                }
                T[j].fit_fun(ff, ud1, ud2);
            }

            // 5. Selekcja (mi + lambda) - Najlepsze przechodz¹ dalej
            solution* total_pop = new solution[mi + lambda];
            for (int k = 0; k < mi; ++k) total_pop[k] = P[k];
            for (int k = 0; k < lambda; ++k) total_pop[mi + k] = T[k];

            for (int k = 0; k < mi + lambda - 1; ++k) {
                for (int m = 0; m < mi + lambda - k - 1; ++m) {
                    if (m2d(total_pop[m].y) > m2d(total_pop[m + 1].y))
                        std::swap(total_pop[m], total_pop[m + 1]);
                }
            }

            for (int k = 0; k < mi; ++k) P[k] = total_pop[k];
            Xopt = P[0];

            delete[] phi; delete[] q; delete[] T; delete[] total_pop;

            if (m2d(Xopt.y) < epsilon) break;
        }
        delete[] P;
        return Xopt;
    }
    catch (string ex) { throw ("EA error: " + ex); }
}