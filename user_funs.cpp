#include"user_funs.h"

#include "solution.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera wartoœæ funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wspó³rzêdne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera wartoœæ funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi¹zujemy równanie ró¿niczkowe
	int n = get_len(Y[0]);									// d³ugoœæ rozwi¹zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad³a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// wartoœæ funkcji celu (ud1 to za³o¿one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po³o¿enia to prêdkoœæ
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z prêdkoœci to przyspieszenie
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
	double val = m2d(x);
	double term1 = -cos(0.1 * val) * exp(-pow(0.1 * val - 2 * M_PI, 2));
	double term2 = 0.002 * pow(0.1 * val, 2);
	return matrix(term1 + term2);
}

long long fib_num(int k)
{
	if (k <= 0) return 0;
	if (k == 1) return 1;

	long long prev = 1; // F1
	long long curr = 1; // F2
	for (int i = 3; i <= k; ++i)
	{
		long long next = prev + curr;
		prev = curr;
		curr = next;
	}
	return curr;
}

matrix ff1S(double t, matrix Y, matrix ud1, matrix ud2) {
	const double g = 9.81;
	const double a = 0.98;
	const double b = 0.63;
	const double PA = 2.0;
	const double PB = 1.0;
	const double TB_IN = 20.0;
	const double FB_IN = 0.01;           // 10 l/s = 0.01 m3/s
	const double DB = 0.00365665;        // 36.5665 cm2 = 0.00365665 m2
	const double epsV = 1e-9;

	// Y = [VA, TA, VB, TB]^T
	double VA = Y(0);
	double TA = Y(1);
	double VB = Y(2);
	double TB = Y(3);

	double DA = m2d(ud1);

	// zabezpieczenia
	if (VA < 0) VA = 0;
	if (VB < 0) VB = 0;

	// PRZEP£YWY
	double hA = (PA > 0 ? VA / PA : 0.0);
	double hB = (PB > 0 ? VB / PB : 0.0);

	double F_out_A = 0.0;
	if (VA > epsV) {
		F_out_A = a * b * DA * sqrt(2.0 * g * hA);
	}

	double F_out_B = 0.0;
	if (VB > epsV) {
		F_out_B = a * b * DB * sqrt(2.0 * g * hB);
	}

	double F_in_A = F_out_A; // przep³yw z A do B

	// dV/dt dla zbiorników
	double dVA_dt = -F_out_A;
	double dVB_dt = F_in_A + FB_IN - F_out_B;

	// ZACHOWAWCZY BILANS ENERGII DLA ZBIORNIKA B
	// d(V*T)/dt = sum(Q_in * T_in) - sum(Q_out * T_out)
	// RHS = F_in_A*TA + FB_IN*TB_IN - F_out_B*TB
	double RHS_energy_B = F_in_A * TA + FB_IN * TB_IN - F_out_B * TB;

	double dVT_dt_B = RHS_energy_B; // to jest d(V*T)/dt
	double dTB_dt;
	if (VB > epsV) {
		// d(V T)/dt = V dT/dt + T dV/dt  => dT/dt = (d(VT)/dt - T dV/dt) / V
		dTB_dt = (dVT_dt_B - TB * dVB_dt) / VB;
	} else {
		dTB_dt = 0.0;
	}

	// Zbiornik A: temperatura trzymana sta³a (brak doplywu), wiêc dTA/dt = 0
	matrix dY(4, 1);
	dY(0) = dVA_dt;
	dY(1) = 0.0;
	dY(2) = dVB_dt;
	dY(3) = dTB_dt;

	return dY;
}

matrix ff1C(matrix x, matrix ud1, matrix ud2) {
	// x = [DA] w [cm^2]
	double DA_cm2 = m2d(x);
	// Konwersja na [m^2]
	double DA_m2 = DA_cm2 / 10000.0;

	// Warunki pocz¹tkowe Y0 = [VA0, TA0, VB0, TB0]^T
	matrix Y0(4, 1);
	Y0(0) = 5.0;   // VA0 = 5 m3
	Y0(1) = 95.0;  // TA0 = 95 C
	Y0(2) = 1.0;   // VB0 = 1 m3
	Y0(3) = 20.0;  // TB0 = 20 C

	// Czas symulacji
	double t0 = 0.0;
	double tend = 2000.0;
	double dt = 1.0;

	// ud1 zawiera [DA] w [m^2]
	matrix UD_DATA(1, 1);
	UD_DATA(0) = DA_m2;

	// jeœli solve_ode ma drugi parametr ud2, przeka¿ pust¹ macierz:
	matrix UD_EMPTY; // albo default constructor
	matrix* Y_sim = solve_ode(ff1S, t0, dt, tend, Y0, UD_DATA, UD_EMPTY);

	// Y_sim[0] = times, Y_sim[1] = states (N x n) po transpozycji w solve_ode
	// Chcemy kolumnê odpowiadaj¹c¹ TB — jeœli indeksacja od 0, TB jest kolumn¹ 3.
	matrix TB_history = get_col(Y_sim[1], 3);

	double max_TB = -1e9;
	int rows = get_size(TB_history)[0];
	for (int i = 0; i < rows; ++i) {
		if (TB_history(i) > max_TB) max_TB = TB_history(i);
	}

	// zwalnianie pamiêci poprawnie
	delete[] Y_sim;

	// Obliczenie wartoœci funkcji celu f(DA) = |max(TB(t)) - 50|
	double objective_value = fabs(max_TB - 50.0);
	return matrix(objective_value);
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
            ej((j+1) % n) = 1.0;
            D_next.set_col(ej, j);

        } else {
            matrix dj = vj * (1.0 / norm_vj);
            D_next.set_col(dj, j);
        }
    }

    return D_next;
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

















