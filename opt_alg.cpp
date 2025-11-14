#include"opt_alg.h"
#include"user_funs.h"

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
