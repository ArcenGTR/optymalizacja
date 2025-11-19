#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto�� funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsp�rz�dne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto�� funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	int n = get_len(Y[0]);									// d�ugo�� rozwi�zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad�a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto�� funkcji celu (ud1 to za�o�one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po�o�enia to pr�dko��
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z pr�dko�ci to przyspieszenie
	return dY;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	double x1 = x(0);
	double x2 = x(1);

	const double PI = 3.14159265358979323846;

	double y = x1 * x1 + x2 * x2
		- std::cos(2.5 * PI * x1)
		- std::cos(2.5 * PI * x2)
		+ 2.0;

	return matrix(y);
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	// x(0) = k1, x(1) = k2 - współczynniki wzmocnienia regulatora
	// ud1(0) = alpha_ref, ud1(1) = omega_ref
	// ud2 nie jest używane w tej funkcji (k1 i k2 są w x)
	
	// Parametry problemu
	double alpha_ref = ud1(0);  // π rad
	double omega_ref = ud1(1);  // 0 rad/s
	
	// Warunki początkowe: [alpha(0), omega(0)]
	matrix Y0(2, 1);
	Y0(0) = 0.0;  // alpha(0) = 0
	Y0(1) = 0.0;  // omega(0) = 0
	
	// Parametry symulacji
	double t0 = 0.0;
	double dt = 0.1;
	double tend = 100.0;
	
	// Rozwiązanie równania różniczkowego
	// ud1 przekazujemy dalej (alpha_ref, omega_ref)
	// x przekazujemy jako ud2 do df2R (k1, k2)
	matrix* Y = solve_ode(df2R, t0, dt, tend, Y0, ud1, x);
	
	// Obliczenie funkcjonału jakości metodą prostokątów
	int n = get_len(Y[0]);
	double Q = 0.0;
	
	for (int i = 0; i < n; ++i)
	{
		double t = Y[0](i);
		double alpha = Y[1](i, 0);
		double omega = Y[1](i, 1);
		
		// Moment siły M(t)
		double M = x(0) * (alpha_ref - alpha) + x(1) * (omega_ref - omega);
		
		// Integrand: 10*(alpha_ref - alpha)^2 + (omega_ref - omega)^2 + M^2
		double integrand = 10.0 * pow(alpha_ref - alpha, 2)
		                 + pow(omega_ref - omega, 2)
		                 + pow(M, 2);
		
		Q += integrand * dt;  // Metoda prostokątów
	}
	
	// Zwolnienie pamięci
	Y[0].~matrix();
	Y[1].~matrix();
	
	return matrix(Q);
}

matrix df2R(double t, matrix Y, matrix ud1, matrix ud2)
{
	// Y(0) = alpha (kąt)
	// Y(1) = omega (prędkość kątowa)
	// ud1(0) = alpha_ref, ud1(1) = omega_ref
	// ud2(0) = k1, ud2(1) = k2
	
	matrix dY(2, 1);
	
	// Parametry fizyczne
	double l = 2.0;      // długość ramienia [m]
	double mr = 1.0;     // masa ramienia [kg]
	double mc = 5.0;     // masa ciężarka [kg]
	double b = 0.25;     // współczynnik tarcia [Nms]
	
	// Moment bezwładności
	double I = (1.0/3.0) * mr * l * l + mc * l * l;
	
	// Wartości referencyjne
	double alpha_ref = ud1(0);
	double omega_ref = ud1(1);
	
	// Współczynniki wzmocnienia
	double k1 = ud2(0);
	double k2 = ud2(1);
	
	// Aktualny stan
	double alpha = Y(0);
	double omega = Y(1);
	
	// Moment siły
	double M = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);
	
	// Równania ruchu
	dY(0) = omega;  // d(alpha)/dt = omega
	dY(1) = (M - b * omega) / I;  // d(omega)/dt = (M - b*omega) / I
	
	return dY;
}