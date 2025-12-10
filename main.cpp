/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include"opt_alg.h"
#include<fstream>
#include <iomanip>
#include"user_funs.h"


void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		//lab0();
		//lab1();
		//lab2();
		//lab3();
		lab4();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dok³adnoœæ
	int Nmax = 10000;										// maksymalna liczba wywo³añ funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz górne ograniczenie
		a(2, 1);											// dok³adne rozwi¹zanie optymalne
	solution opt;											// rozwi¹zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo³anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Wahadlo
	Nmax = 1000;											// dok³adnoœæ
	epsilon = 1e-2;											// maksymalna liczba wywo³añ funkcji celu
	lb = 0, ub = 5;											// dolne oraz górne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad³a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo³anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi¹zujemy równanie ró¿niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumieñ do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumieñ
	Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	Y[1].~matrix();
}

void lab1()
{

	// Testowa funkcja 100 w csv

	/*
	double x0 = 50;
	double d = 1.0;
	double alpha = 1.1;
	int Nmax_exp = 100;

	// Krok 1: Wyznaczenie przedzia³u przez ekspansjê
	solution::clear_calls();
	double* p = expansion(ff1T, x0, d, alpha, Nmax_exp);
	double a_start = p[0];
	double b_start = p[1];
	a_start = -100;
	b_start = 100;
	delete[] p;

	cout << "\n--- MINIMALIZACJA JEDNOWYMIAROWA ---" << endl;
	cout << "Przedzial startowy: [" << a_start << ", " << b_start << "]" << endl;

	// --- TEST METODY FIBONACCIEGO ---
	double epsilon_fib = 1e-4; // Dok³adnoœæ
	solution::clear_calls();
	solution opt_fib = fib(ff1T, a_start, b_start, epsilon_fib);

	cout << "\n--- Metoda Fibonacciego ---" << endl;
	cout << "Optymalne x*: " << m2d(opt_fib.x) << endl;
	cout << "Wartosc funkcji: " << m2d(opt_fib.y) << endl;
	cout << "Liczba wywolan f. celu: " << solution::f_calls << endl;

	// --- TEST METODY LAGRANGE'A ---
	double epsilon_lag = 1e-4;
	double gamma_lag = 1e-4;
	int Nmax_lag = 1000;

	solution::clear_calls();
	solution opt_lag = lag(ff1T, a_start, b_start, epsilon_lag, gamma_lag, Nmax_lag);

	cout << "\n--- Metoda Lagrange'a ---" << endl;
	cout << "Optymalne x*: " << m2d(opt_lag.x) << endl;
	cout << "Wartosc funkcji: " << m2d(opt_lag.y) << endl;
	cout << "Liczba wywolan f. celu: " << solution::f_calls << endl;

	*/


	// Model rzeczwisty



	// --- CZÊŒÆ 1: TESTOWANIE MODELU ---
    cout << "--- TEST MODELU ZBIORNIKOW ---" << endl;
    double DA_test = 50.0; // cm^2
    matrix DA_mat(DA_test);

    // POPRAWKA: U¿yto Twojej funkcji celu ff1C zamiast ff_zadanie
    matrix result_obj = ff1C(DA_mat, NAN, NAN);

    // Rêczne uruchomienie symulacji do weryfikacji max TB
    double DA_m2 = DA_test / 10000.0;
    matrix Y0(4, 1);
    Y0(0) = 5.0; Y0(1) = 95.0; Y0(2) = 1.0; Y0(3) = 20.0;
    matrix UD_DATA(DA_m2);

    matrix* Y_sim = solve_ode(ff1S, 0.0, 1.0, 2000.0, Y0, UD_DATA);

    matrix TB_history = get_col(Y_sim[1], 3);

    double max_TB = -1.0;
    int rows = get_size(TB_history)[0];
    for (int i = 0; i < rows; ++i)
    {
        if (TB_history(i) > max_TB) max_TB = TB_history(i);
    }

    Y_sim[0].~matrix();
    Y_sim[1].~matrix();

    cout << "DA = " << DA_test << " cm^2" << endl;
    cout << "Maksymalna temperatura w B (Max TB): " << max_TB << " C" << endl;
    cout << "Oczekiwana: ok. 79.14 C" << endl;
    cout << "Wartosc funkcji celu |max_TB - 50|: " << m2d(result_obj) << endl;
    cout << "---------------------------------" << endl;

    // --- CZÊŒÆ 2: OPTYMALIZACJA ---
    cout << "\n--- OPTYMALIZACJA ZBIORNIKOW ---" << endl;

    double a_start = 1.0;
    double b_start = 100.0;
    double epsilon = 1e-3;
    double gamma = 1e-3;
    int Nmax = 500;

    // 1. U¿ycie metody LAGRANGE'A
    solution::clear_calls();
    solution opt_lag = lag(ff1C, a_start, b_start, epsilon, gamma, Nmax);

    cout << "\n--- Wynik Optymalizacji (Lagrange) ---" << endl;
    cout << "Optymalne DA: " << m2d(opt_lag.x) << " cm^2" << endl;
    // Aby uzyskaæ max TB, dodajemy 50 do wartoœci funkcji celu (bo f=|maxTB-50|)
    cout << "Maksymalna uzyskana temp. w B: " << m2d(opt_lag.y) + 50.0 << " C" << endl;
    cout << "Wartosc funkcji celu f(DA): " << m2d(opt_lag.y) << endl;
    cout << "Liczba wywolan f. celu: " << solution::f_calls << endl;

    // 2. U¿ycie metody FIBONACCIEGO
    solution::clear_calls();
    // POPRAWKA: U¿yto Twojej funkcji celu ff1C zamiast ff_zadanie
    solution opt_fib = fib(ff1C, a_start, b_start, epsilon);

    cout << "\n--- Wynik Optymalizacji (Fibonacci) ---" << endl;
    cout << "Optymalne DA: " << m2d(opt_fib.x) << " cm^2" << endl;
    cout << "Maksymalna uzyskana temp. w B: " << m2d(opt_fib.y) + 50.0 << " C" << endl;
    cout << "Wartosc funkcji celu f(DA): " << m2d(opt_fib.y) << endl;
    cout << "Liczba wywolan f. celu: " << solution::f_calls << endl;



	//Symulacja z optymalnymi wynikami

	/*
	// ??????????? ???????? DA, ????????? ? ??????? (? ??^2)
	const double DA_OPT_FIB_CM2 = 18.8536;
	const double DA_OPT_LAG_CM2 = 20.138;

	// ??????????? ? ?^2 ??? ?????? ODE (ff1S)
	const double DA_OPT_FIB_M2 = DA_OPT_FIB_CM2 / 10000.0;
	const double DA_OPT_LAG_M2 = DA_OPT_LAG_CM2 / 10000.0;

	// ????? ????????? ?????????
	const double T0 = 0.0;
	const double T_END = 2000.0;
	const double DT = 1.0;
	const int STEP_SIZE = 20; // ??? ??????? ??? ????? (?????? 20 ??????)

	// ????????? ??????? Y0 = [VA0, TA0, VB0, TB0]^T
	matrix Y0(4, 1);
	Y0(0) = 5.0;
	Y0(1) = 95.0;
	Y0(2) = 1.0;
	Y0(3) = 20.0;

	// ?????? ????????? ??? DA_Fibonacci
	matrix UD_DATA_FIB(1, 1);
	UD_DATA_FIB(0) = DA_OPT_FIB_M2;
	matrix* Y_sim_FIB = solve_ode(ff1S, T0, DT, T_END, Y0, UD_DATA_FIB);

	// ?????? ????????? ??? DA_Lagrange
	matrix UD_DATA_LAG(1, 1);
	UD_DATA_LAG(0) = DA_OPT_LAG_M2;
	matrix* Y_sim_LAG = solve_ode(ff1S, T0, DT, T_END, Y0, UD_DATA_LAG);

	// Y_sim[0] = ????? t (Nx1), Y_sim[1] = ????????? [VA, TA, VB, TB] (Nx4)

	// ----------------------------------------------------------------------
	// 2. ????? ??????????? ? ???? ? ???????? ???????
	// ----------------------------------------------------------------------

	try {
	    std::ofstream file("symulacja_opt.csv");
	    if (!file.is_open()) {
	        throw std::string("Nie mozna otworzyc pliku symulacja_opt.csv");
	    }

	    // ????????? ?????????????? (6 ?????? ????? ???????)
	    file << std::fixed << std::setprecision(6);

	    // ????????? ? ????? ???????? ???????
	    file << "t;VAFibonacci;VALagrange;VBFibonacci;VBLagrange;TBFibonacci;TBLagrange\n";

	    // ????????? ???????? (?????? ???? ????????? ??? ????? ?????????)
	    int rows_count = get_size(Y_sim_FIB[0])[0];

	    // ???????? ? ????? 20
	    for (int i = 0; i < rows_count; i += STEP_SIZE) {

	        // ????? (Y_sim[0] ???????? ??????? ???????)
	        file << m2d(Y_sim_FIB[0](i));

	        // VAFibonacci (Y_sim[1] col 0)
	        file << ";" << m2d(Y_sim_FIB[1](i, 0));

	        // VALagrange (Y_sim[1] col 0)
	        file << ";" << m2d(Y_sim_LAG[1](i, 0));

	        // VBFibonacci (Y_sim[1] col 2)
	        file << ";" << m2d(Y_sim_FIB[1](i, 2));

	        // VBLagrange (Y_sim[1] col 2)
	        file << ";" << m2d(Y_sim_LAG[1](i, 2));

	        // TBFibonacci (Y_sim[1] col 3)
	        file << ";" << m2d(Y_sim_FIB[1](i, 3));

	        // TBLagrange (Y_sim[1] col 3)
	        file << ";" << m2d(Y_sim_LAG[1](i, 3));

	        file << "\n"; // ??????? ?? ????? ??????
	    }

	    file.close();
	    std::cout << "\nWyniki symulacji dla optymalnych DA zapisano do pliku: symulacja_opt.csv\n";
	    std::cout << "Probkowanie co 20 sekund.\n";
	}
	catch (std::string ex_info) {
	    std::cerr << "Blad w trakcie zapisu pliku: " << ex_info << std::endl;
	}

	// ----------------------------------------------------------------------
	// 3. ???????????? ??????
	// ----------------------------------------------------------------------
	delete[] Y_sim_FIB;
	delete[] Y_sim_LAG;
	*/
}

void lab2()
{
	try
	{
		// Funkcja testowa

		solution::clear_calls();

		int N = 2;

		double x0_data[] = { 0.4, -0.4 }; // Punkt startowy
		matrix x0(N, x0_data);

		double s0_data[] = { 0.1, 0.1 }; // Pocz¹tkowe kroki
		matrix s0(N, s0_data);

		double alpha = 3.0;     // Wspó³czynnik ekspansji
		double beta = 0.5;      // Wspó³czynnik kontrakcji
		double epsilon = 1e-4;  // Dok³adnoœæ
		int Nmax = 1000;        // Maksymalna liczba wywo³añ funkcji celu

		matrix ud1(NAN), ud2(NAN);

		std::cout << "--- Rozpoczêcie optymalizacji metod¹ Rosenbrocka ---\n";
		std::cout << "Punkt startowy x0: " << x0 << "\n";

		solution Xopt = Rosen(ff2T, x0, s0, alpha, beta, epsilon, Nmax, ud1, ud2);

		std::cout << "\n--- Wyniki optymalizacji ---\n";
		std::cout << Xopt;

		std::cout << "\nKomentarz do flagi wyjœcia:\n";
		if (Xopt.flag == 2) {
			std::cout << "Flaga 2: Algorytm zbieg³ do rozwi¹zania z zadan¹ dok³adnoœci¹ (max(|s_j|) < epsilon).\n";
		} else if (Xopt.flag == 0) {
			std::cout << "Flaga 0: Przekroczono maksymaln¹ liczbê wywo³añ funkcji celu (Nmax).\n";
		} else {
			std::cout << "Flaga " << Xopt.flag << ": Nieznany status wyjœcia.\n";
		}

	}
	catch (string ex_info)
	{
		std::cerr << "Wyst¹pi³ b³¹d w run_rosen_test():\n" << ex_info << "\n";
	}
}

void lab3()
{
    try
    {
       // std::cout << "--- LAB 3: Nelder-Mead with Penalty Functions ---" << endl;
       //
       // // Experiment parameters
       // double epsilon = 1e-3;
       // int Nmax = 5000;
       // int attempts = 100;
       //
       // std::vector<double> a_vals = { 4.0, 4.4934, 5.0 };
       // // Przechowamy punkty startowe, aby wykorzystaæ je do obu typów kary
       // std::vector<matrix> start_points(attempts);
       //
       // // Setup output file
       // ofstream out("lab3_results.csv");
       // // Nag³ówek: a;PenaltyType;x1_pocz;x2_pocz;x1_opt;x2_opt;r_from_origin;y_opt;f_calls
       // out << "a;PenaltyType;x1_pocz;x2_pocz;x1_opt;x2_opt;r_from_origin;y_opt;f_calls\n";
       //
       // // Iterate over 'a'
       // for (double a : a_vals)
       // {
       //    matrix ud2(1, 1);
       //    ud2(0) = a;
       //
       //    std::cout << "\n========================================\n";
       //    std::cout << "Starting tests for a=" << a << "\n";
       //    std::cout << "========================================\n";
       //
       //    // 1. GENEROWANIE 100 LOSOWYCH PUNKTÓW STARTOWYCH (jednokrotnie dla danego 'a')
       //    std::cout << "Generating 100 feasible start points...\n";
       //    for (int i = 0; i < attempts; ++i)
       //    {
       //       matrix x0(2, 1);
       //       bool feasible = false;
       //       while (!feasible)
       //       {
       //          // Losowanie w szerokim zakresie, potem przeskalowanie lub odrzucenie
       //          x0 = rand_mat(2, 1) * 5.0;
       //
       //          // Warunki dopuszczalnoœci: x1 >= 1, x2 >= 1, r <= a
       //          // Aby punkt nadawa³ siê do kary wewnêtrznej, musi byæ ŒCIŒLE wewnêtrzny.
       //          double r = norm(x0);
       //
       //          if (x0(0) > 1.0 + 1e-6 && x0(1) > 1.0 + 1e-6 && r < a - 1e-6)
       //          {
       //             feasible = true;
       //          }
       //       }
       //       start_points[i] = x0;
       //    }
       //
       //    // 2. PÊTLA PO TYPACH KARY (grupowanie wed³ug typu kary)
       //    for (int pen_type = 1; pen_type <= 2; ++pen_type)
       //    {
       //       std::string pen_name = (pen_type == 1 ? "External" : "Internal");
       //       std::cout << "--- Testing Penalty Type: " << pen_name << " ---\n";
       //
       //       matrix ud1(2, 1); // [c, type]
       //       ud1(1) = (double)pen_type;
       //
       //       // Stats for averaging
       //       double sum_x1 = 0, sum_x2 = 0, sum_y = 0, sum_calls = 0, sum_r = 0;
       //
       //       // PARAMETRY KARY
       //       double c_start = (pen_type == 1) ? 1.0 : 5.0;
       //       double dc = 2.0;
       //
       //       // 3. PÊTLA PO 100 ZAPISANYCH PUNKTACH STARTOWYCH
       //       for (int i = 0; i < attempts; ++i)
       //       {
       //          matrix x0 = start_points[i];
       //          solution::clear_calls();
       //
       //          double x1_pocz = m2d(x0(0));
       //          double x2_pocz = m2d(x0(1));
       //
       //          ud1(0) = c_start;
       //
       //          // Run optimization
       //          solution res = pen(ff3T, x0, c_start, dc, epsilon, Nmax, ud1, ud2);
       //
       //          // Collect data
       //          double r_res = norm(res.x);
       //
       //          // Zapis do pliku CSV w wymaganej kolejnoœci
       //          out << a << ";" << pen_type << ";"
       //             << x1_pocz << ";" << x2_pocz << ";"
       //             << m2d(res.x(0)) << ";" << m2d(res.x(1)) << ";"
       //             << r_res << ";"
       //             << m2d(res.y) << ";"
       //             << solution::f_calls << "\n";
       //
       //          sum_x1 += m2d(res.x(0));
       //          sum_x2 += m2d(res.x(1));
       //          sum_y += m2d(res.y);
       //          sum_calls += solution::f_calls;
       //          sum_r += r_res;
       //       }
       //
       //       // Print average stats to console
       //       std::cout << "Avg results for a=" << a << ", Type=" << pen_name << ":\n";
       //       std::cout << "  x1: " << sum_x1/attempts << ", x2: " << sum_x2/attempts << "\n";
       //       std::cout << "  f(x): " << sum_y/attempts << "\n";
       //       std::cout << "  calls: " << sum_calls/attempts << "\n";
       //       std::cout << "  r: " << sum_r/attempts << "\n\n";
       //    }
       // }
       // out.close();
       // std::cout << "Results saved to lab3_results.csv" << endl;



    	std::cout << "--- Optimization of Magnus Effect Ball Trajectory ---" << endl;

		// --- 1. MODEL VERIFICATION ---
		// Parameters from description: v0x = 5, omega = 10
		// Expected: x_end approx 41.41, x_50 approx 21.61
		std::cout << "\n[Verification Step]\n";
		matrix x_test(2, 1);
		x_test(0) = 5.0;
		x_test(1) = 10.0;

		// We call the objective function with c=0 to check pure simulation results
		// Note: The objective function calculates penalties but we print the raw simulation values below manually to be sure.

		// Run manual simulation for verification printing
		matrix Y0_test(4, 1);
		Y0_test(0) = 0; Y0_test(1) = 5.0; Y0_test(2) = 100; Y0_test(3) = 0;
		matrix ode_ud_test(1, 1);
		ode_ud_test(0) = 10.0; // omega
		matrix* Y_ver = solve_ode(df4, 0, 0.2, 7, Y0_test, ode_ud_test, NAN);

		int n = get_size(Y_ver[0])[0];
		double x_50_ver = 0, x_end_ver = 0;
		for (int i = 0; i < n - 1; ++i) {
			if (Y_ver[1](i, 2) >= 50 && Y_ver[1](i + 1, 2) < 50) {
				double f = (50 - Y_ver[1](i, 2)) / (Y_ver[1](i + 1, 2) - Y_ver[1](i, 2));
				x_50_ver = Y_ver[1](i, 0) + f * (Y_ver[1](i + 1, 0) - Y_ver[1](i, 0));
			}
			if (Y_ver[1](i, 2) >= 0 && Y_ver[1](i + 1, 2) < 0) {
				double f = (0 - Y_ver[1](i, 2)) / (Y_ver[1](i + 1, 2) - Y_ver[1](i, 2));
				x_end_ver = Y_ver[1](i, 0) + f * (Y_ver[1](i + 1, 0) - Y_ver[1](i, 0));
			}
		}
		std::cout << "Test parameters: v0x=5, omega=10\n";
		std::cout << "Calculated x_end: " << x_end_ver << " (Expected ~41.41)\n";
		std::cout << "Calculated x_50:  " << x_50_ver << " (Expected ~21.61)\n";
		delete[] Y_ver;


		// --- 2. OPTIMIZATION ---
		std::cout << "\n[Optimization Start]\n";

		// Initial guess (within bounds [-10, 10])
		matrix x0(2, 1);
		x0(0) = 1.0; // v0x start
		x0(1) = 1.0; // omega start

		// Penalty method parameters
		double c_start = 1.0;
		double dc = 2.0;       // Penalty increase factor
		double epsilon = 1e-5; // Tolerance
		int Nmax = 5000;       // Max function calls

		matrix ud1(1, 1); // ud1(0) will be updated by 'pen' function with current 'c'
		matrix ud2;       // Empty

		solution::clear_calls();

		// Use External Penalty function (implicit in how we wrote ff4R and using 'pen')
		// The 'pen' function calls sym_NM (Nelder-Mead) internally.
		solution opt = pen(ff4R, x0, c_start, dc, epsilon, Nmax, ud1, ud2);

		std::cout << "\n[Optimization Results]\n";
		std::cout << "Optimal v0x:   " << m2d(opt.x(0)) << " m/s\n";
		std::cout << "Optimal omega: " << m2d(opt.x(1)) << " rad/s\n";

		// We need to re-run simulation to show the physical result (since opt.y includes penalty)
		matrix final_Y0(4, 1);
		final_Y0(0) = 0; final_Y0(1) = m2d(opt.x(0)); final_Y0(2) = 100; final_Y0(3) = 0;
		matrix final_omega(1, 1);
		final_omega(0) = m2d(opt.x(1));
		matrix* Y_final = solve_ode(df4, 0, 0.2, 7, final_Y0, final_omega, NAN);

		// Find final metrics for display
		double x_end_opt = 0, x_50_opt = 0;
		n = get_size(Y_final[0])[0];
		for (int i = 0; i < n - 1; ++i) {
			if (Y_final[1](i, 2) >= 50 && Y_final[1](i + 1, 2) < 50) {
				double f = (50 - Y_final[1](i, 2)) / (Y_final[1](i + 1, 2) - Y_final[1](i, 2));
				x_50_opt = Y_final[1](i, 0) + f * (Y_final[1](i + 1, 0) - Y_final[1](i, 0));
			}
			if (Y_final[1](i, 2) >= 0 && Y_final[1](i + 1, 2) < 0) {
				double f = (0 - Y_final[1](i, 2)) / (Y_final[1](i + 1, 2) - Y_final[1](i, 2));
				x_end_opt = Y_final[1](i, 0) + f * (Y_final[1](i + 1, 0) - Y_final[1](i, 0));
			}
		}

		std::cout << "Maximized x_end: " << x_end_opt << " m\n";
		std::cout << "Constraint check (x at y=50): " << x_50_opt << " m (Target: [3, 7])\n";
		std::cout << "Function calls: " << solution::f_calls << endl;

		// --- 3. SAVE SIMULATION TO CSV ---
		ofstream Sout("symulacja_lab4.csv");
		Sout << "t;x;vx;y;vy\n"; // Header
		Sout << hcat(Y_final[0], Y_final[1]);
		Sout.close();
		std::cout << "Simulation data saved to 'symulacja_lab4.csv'\n";

		delete[] Y_final;

    }
    catch (string ex_info)
    {
       cerr << "ERROR in lab3:\n" << ex_info << endl;
    }
}

void lab4()
{
    // // Funkcja celu: f(x1, x2) = 1/6*x1^6 - 1.05*x1^4 + 2*x1^2 + x2^2 + x1*x2
    //
    // // Parametry symulacji
    // const int N_TRIALS = 100;
    // double epsilon = 1e-4;
    // int Nmax = 1000;
    // matrix ud1, ud2; // Puste dane u¿ytkownika
    //
    // // Ograniczenia dla punktu startowego: x1, x2 nale¿¹ do [-2, 2]
    // matrix lb(2, 1, -2.0), ub(2, 1, 2.0);
    // double range_x1 = ub(0) - lb(0); // 4.0
    // double range_x2 = ub(1) - lb(1); // 4.0
    //
    // // Szukane minimum g³êbokie lokalne (wartoœæ f)
    // // f* ? 0.298638 (dla x* ? [±1.7475, ?0.8737]^T)
    // const double MINIMUM_TARGET = 0.298638;
    // const double TOLERANCE = 0.001; // Tolerancja dla uznania, ¿e znaleziono minimum g³êbokie
    //
    // // Ustawienie ziarna losowoœci
    // srand(time(NULL));
    //
    // // Struktura do zbierania statystyk dla ka¿dej serii 100 prób
    // struct RunStats {
    //     // Sumy dla wszystkich prób
    //     double sum_x1_start = 0.0;
    //     double sum_x2_start = 0.0;
    //
    //     // Sumy tylko dla UDANYCH optymalizacji (flag=1 i nie NaN)
    //     double sum_x1_final = 0.0;
    //     double sum_x2_final = 0.0;
    //     double sum_y_final = 0.0;
    //     long long sum_f_calls = 0;
    //     long long sum_g_calls = 0;
    //     long long sum_H_calls = 0;
    //
    //     // Liczniki
    //     int total_runs = 0; // Zawsze 100
    //     int successful_runs = 0; // Udane konwergencje (flag = 1 i nie NaN)
    //     int global_min_count = 0; // Znalezione minimum g³êbokie (y* ? 0.298)
    // };
    //
    // // --- 1. Generacja 100 punktów startowych (jednorazowo) ---
    // std::vector<matrix> start_points(N_TRIALS);
    // for (int t = 0; t < N_TRIALS; ++t) {
    //     matrix x0(2, 1);
    //     x0(0) = lb(0) + range_x1 * ((double)rand() / RAND_MAX);
    //     x0(1) = lb(1) + range_x2 * ((double)rand() / RAND_MAX);
    //     start_points[t] = x0;
    // }
    // std::cout << "Wygenerowano " << N_TRIALS << " losowych punktów startowych w zakresie [-2, 2].\n";
    // // ------------------------------------------------------------------
    //
    // // Definicja metod i kroków
    // std::vector<std::string> methods = {"SD", "CG", "Newton"};
    // std::vector<double> steps = {0.05, 0.25, 0.0}; // 0.0 oznacza krok zmienny (Golden)
    //
    // // Otwarcie pliku do zapisu wyników
    // ofstream results_file("lab4_results.csv");
    // results_file << std::fixed << std::setprecision(6);
    //
    // // Nag³ówek w wymaganym formacie
    // results_file << "Metoda;Dlugosc_kroku;x1_0;x2_0;x1_gwiazdka;x2_gwiazdka;y_gwiazdka;f_calls;g_calls;H_calls;Minimum_Globalne[TAK/NIE]\n";
    //
    // std::cout << "--- LAB 4: Symulacja Monte Carlo (" << N_TRIALS << " prób) ---\n";
    //
    // // --- G³ówny cykl: Iteracja po d³ugoœciach kroku ---
    // for (double h : steps) {
    //
    //     // --- Wewnêtrzny cykl: Iteracja po metodach ---
    //     for (const std::string& method : methods) {
    //
    //         RunStats stats; // Struktura do zbierania statystyk dla bie¿¹cej grupy (Metoda, Krok)
    //
    //         // --- Najbardziej wewnêtrzny cykl: Iteracja po 100 punktach startowych ---
    //         for (int t = 0; t < N_TRIALS; ++t) {
    //             matrix x0 = start_points[t]; // Pobranie wczeœniej wylosowanego punktu
    //
    //             stats.total_runs++;
    //             stats.sum_x1_start += x0(0);
    //             stats.sum_x2_start += x0(1);
    //
    //             solution::clear_calls();
    //             solution opt;
    //             bool error_occurred = false;
    //
    //             try {
    //                 if (method == "SD") {
    //                     opt = SD(ff5T, gf5T, x0, h, epsilon, Nmax, ud1, ud2);
    //                 } else if (method == "CG") {
    //                     opt = CG(ff5T, gf5T, x0, h, epsilon, Nmax, ud1, ud2);
    //                 } else if (method == "Newton") {
    //                     opt = Newton(ff5T, gf5T, Hf5T, x0, h, epsilon, Nmax, ud1, ud2);
    //                 }
    //
    //                 // Sprawdzenie, czy optymalizacja by³a udana (flag=1) i wynik jest skoñczony (nie NaN)
    //                 if (opt.flag != 1 || !std::isfinite(m2d(opt.y))) {
    //                     error_occurred = true;
    //                 }
    //
    //                 // Okreœlenie, czy znaleziono minimum globalne
    //                 bool found_global = !error_occurred && std::fabs(m2d(opt.y) - MINIMUM_TARGET) < TOLERANCE;
    //
    //                 // --- Zbieranie statystyk dla udanych prób ---
    //                 if (!error_occurred) {
    //                     stats.successful_runs++;
    //                     stats.sum_x1_final += opt.x(0);
    //                     stats.sum_x2_final += opt.x(1);
    //                     stats.sum_y_final += m2d(opt.y);
    //                     stats.sum_f_calls += solution::f_calls;
    //                     stats.sum_g_calls += solution::g_calls;
    //                     stats.sum_H_calls += solution::H_calls;
    //
    //                     // Licznik znalezienia minimum globalnego
    //                     if (found_global) {
    //                         stats.global_min_count++;
    //                     }
    //                 }
    //
    //                 // Zapis wyników do pliku CSV
    //                 results_file << method << ";" << h << ";"
    //                              << x0(0) << ";" << x0(1) << ";"
    //                              << (error_occurred ? "NaN" : std::to_string(opt.x(0))) << ";"
    //                              << (error_occurred ? "NaN" : std::to_string(opt.x(1))) << ";"
    //                              << (error_occurred ? "NaN" : std::to_string(m2d(opt.y))) << ";"
    //                              << solution::f_calls << ";"
    //                              << solution::g_calls << ";"
    //                              << solution::H_calls << ";"
    //                              << (found_global ? "TAK" : "NIE") << "\n";
    //
    //             } catch (string ex) {
    //                 // Zapis b³êdu (np. b³¹d w inwersji macierzy)
    //                 results_file << method << ";" << h << ";"
    //                              << x0(0) << ";" << x0(1) << ";"
    //                              << "Error;Error;Error;"
    //                              << solution::f_calls << ";"
    //                              << solution::g_calls << ";"
    //                              << solution::H_calls << ";"
    //                              << "NIE" << "\n";
    //             }
    //         } // Koniec 100 prób dla danego (Metoda, Krok)
    //
    //         // --- WYŒWIETLENIE ZBIORCZEJ STATYSTYKI DO KONSOLI ---
    //
    //         // Formatowanie wyœwietlania kroku
    //         std::string h_str = (h == 0.0) ? "Zmienny (Z³oty Podzia³)" : std::to_string(h);
    //
    //         std::cout << std::fixed << std::setprecision(5);
    //         std::cout << "\n======================================================\n";
    //         std::cout << "PODSUMOWANIE STATYSTYCZNE\n";
    //         std::cout << "Metoda: " << method << ", D³ugoœæ kroku (h): " << h_str << "\n";
    //         std::cout << "Liczba prób (N): " << stats.total_runs << "\n";
    //         std::cout << "Liczba udanych konwergencji: " << stats.successful_runs << "\n";
    //         std::cout << "Liczba znalezionych minimów globalnych: " << stats.global_min_count << "\n";
    //         std::cout << "------------------------------------------------------\n";
    //
    //         if (stats.total_runs > 0) {
    //             std::cout << "Œrednia x1_0:         " << (stats.sum_x1_start / stats.total_runs) << "\n";
    //             std::cout << "Œrednia x2_0:         " << (stats.sum_x2_start / stats.total_runs) << "\n";
    //         }
    //
    //         if (stats.successful_runs > 0) {
    //             std::cout << "\nŒREDNIE WYNIKI DLA UDANYCH OPTYMALIZACJI:\n";
    //             std::cout << "Œrednia x1_gwiazdka:  " << (stats.sum_x1_final / stats.successful_runs) << "\n";
    //             std::cout << "Œrednia x2_gwiazdka:  " << (stats.sum_x2_final / stats.successful_runs) << "\n";
    //             std::cout << "Œrednia y_gwiazdka:   " << (stats.sum_y_final / stats.successful_runs) << "\n";
    //             std::cout << "Œrednia f_calls:      " << (double)stats.sum_f_calls / stats.successful_runs << "\n";
    //             std::cout << "Œrednia g_calls:      " << (double)stats.sum_g_calls / stats.successful_runs << "\n";
    //             std::cout << "Œrednia H_calls:      " << (double)stats.sum_H_calls / stats.successful_runs << "\n";
    //         } else {
    //             std::cout << "Brak udanych konwergencji do obliczenia œrednich wyników.\n";
    //         }
    //         std::cout << "======================================================\n";
    //         // -------------------------------------------------------------------
    //
    //     } // Koniec metod
    // } // Koniec kroków
    //
    // results_file.close();
    // std::cout << "\nSymulacja zakoñczona. Wyniki zapisano do lab4_results.csv" << endl;



	// double start_x1 = -1.1951;
 //    double start_x2 = 0.666415;
 //
 //    // Punkt startowy
 //    matrix x0(2, 1);
 //    x0(0) = start_x1;
 //    x0(1) = start_x2;
 //
 //    double epsilon = 1e-4;
 //    int Nmax = 1000;
 //    matrix ud1, ud2;
 //
 //    std::vector<std::string> methods = {"SD", "CG", "Newton"};
 //    std::vector<double> steps = {0.05, 0.25, 0.0}; // 0.0 = zmienny krok
 //
 //    // Otwarcie pliku do zapisu trajektorii
 //    ofstream trajectory_file("trajectory_all_steps.csv"); // Zmieniono nazwê pliku
 //    trajectory_file << std::fixed << std::setprecision(8);
 //
 //    // Nag³ówek pliku CSV dla trajektorii
 //    trajectory_file << "Metoda;Krok_h;Iteracja;x1;x2;y_val;f_calls;g_calls;H_calls\n";
 //
 //    std::cout << "\n--- LAB 6: Generowanie trajektorii dla x0 = [" << start_x1 << ", " << start_x2 << "] ---\n";
 //
 //    // Pêtla po metodach i krokach
 //    for (double h : steps) {
 //        for (const std::string& method : methods) {
 //
 //            solution::clear_calls();
 //            solution opt;
 //            bool error_occurred = false;
 //
 //            // --- ZAPIS PUNKTU STARTOWEGO (Iteracja 0) ---
 //            trajectory_file << method << ";" << h << ";" << 0 << ";"
 //                            << x0(0) << ";" << x0(1) << ";"
 //                            << m2d(ff5T(x0, ud1, ud2)) << ";" // Obliczenie y dla x0
 //                            << 0 << ";" << 0 << ";" << 0 << "\n";
 //            // ------------------------------------------
 //
 //            try {
 //                if (method == "SD") {
 //                    opt = SD_Logged(ff5T, gf5T, x0, h, epsilon, Nmax, ud1, ud2, trajectory_file, method, h);
 //                } else if (method == "CG") {
 //                    opt = CG_Logged(ff5T, gf5T, x0, h, epsilon, Nmax, ud1, ud2, trajectory_file, method, h);
 //                } else if (method == "Newton") {
 //                    opt = Newton_Logged(ff5T, gf5T, Hf5T, x0, h, epsilon, Nmax, ud1, ud2, trajectory_file, method, h);
 //                }
 //
 //                if (opt.flag == 1 && std::isfinite(m2d(opt.y))) {
 //                    std::cout << "  Zapisano trajektoriê (sukces) dla: " << method << ", h=" << h << std::endl;
 //                } else {
 //                    error_occurred = true;
 //                    std::cout << "  B£¥D/Niepowodzenie konwergencji dla: " << method << ", h=" << h << std::endl;
 //                }
 //
 //            } catch (string ex) {
 //                error_occurred = true;
 //                std::cerr << "  Wyst¹pi³ wyj¹tek dla " << method << ", h=" << h << ": " << ex << std::endl;
 //            }
 //        }
 //    }
 //
 //    trajectory_file.close();
 //    std::cout << "\nZakoñczono generowanie pliku trajectory_all_steps.csv.\n";
 //    std::cout << "PLIK ZAWIERA PE£NE TRAJEKTORIE DLA KA¯DEJ KONFIGURACJI.\n";



	std::cout << "\n--- Regresja Logistyczna (Tylko CG) ---\n";

    // 1. Wczytanie danych
    matrix X_data;
    matrix Y_data;
    try {
        // X: 3x100 (theta0, x1, x2)
        X_data = read_data_matrix("XData.txt", 3, 100);
        // Y: 1x100
        Y_data = read_data_matrix("YData.txt", 1, 100);

    } catch (const std::string& ex) {
        std::cerr << "B£¥D: Nie mo¿na wczytaæ plików XData.txt lub YData.txt:\n" << ex << std::endl;
        return;
    }

    // Dane dla algorytmu optymalizacyjnego
    matrix x0(3, 1, 0.0); // Wektor startowy theta(0) = [0, 0, 0]^T
    double epsilon = 1e-4;
    int Nmax = 10000;
    matrix ud_X = X_data; // Dane u¿ytkownika 1 (X)
    matrix ud_Y = Y_data; // Dane u¿ytkownika 2 (Y)

    // Konfiguracje optymalizacji (Tylko CG, wymagane kroki sta³e)
    std::vector<std::string> methods = {"CG"};
    std::vector<double> steps = {0.01, 0.001, 0.0001};

    // Nag³ówek dla Tabeli 3 (na konsolê)
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\n-----------------------------------------------------------------------------------------------------------------\n";
    std::cout << "| Metoda | Krok_h | J(theta*) | P(theta*) [%] | f_calls | g_calls | H_calls | theta0 | theta1 | theta2 |\n";
    std::cout << "-----------------------------------------------------------------------------------------------------------------\n";

    solution best_opt;
    double best_accuracy = 0.0;
    std::string best_case_info;
    int best_g_calls = 0; // Zmienna do przechowania liczby wywo³añ gradientu dla najlepszego przypadku

    for (double h : steps) {
        for (const std::string& method : methods) {

            solution::clear_calls();
            solution opt;

            try {
                // Uruchomienie metody CG
                if (method == "CG") {
                    opt = CG(ff_logistic, gf_logistic, x0, h, epsilon, Nmax, ud_X, ud_Y);
                }

                if (opt.flag == 1) {
                    double accuracy = classification_accuracy(opt.x, X_data, Y_data);

                    // Standardowy wydruk do tabeli
                    std::cout << "| " << std::setw(6) << method << " | " << std::setw(6) << h << " | "
                              << m2d(opt.y) << " | " << std::setw(11) << accuracy << " | "
                              << std::setw(7) << solution::f_calls << " | " << std::setw(7) << solution::g_calls << " | "
                              << std::setw(7) << solution::H_calls << " | "
                              << opt.x(0) << " | " << opt.x(1) << " | " << opt.x(2) << " |\n";

                    // Znajdowanie i zapisywanie najlepszego przypadku
                    if (accuracy > best_accuracy) {
                        best_accuracy = accuracy;
                        best_opt = opt;
                        best_case_info = method + ", h=" + std::to_string(h);
                        best_g_calls = solution::g_calls; // KLUCZOWE: Zapisujemy g_calls
                    }
                } else {
                     std::cout << "| " << std::setw(6) << method << " | " << std::setw(6) << h << " | BRAK KONWERGENCJI (flaga=" << opt.flag << ") |\n";
                }

            } catch (std::string ex) {
                std::cerr << "  Wyst¹pi³ wyj¹tek dla " << method << ", h=" << h << ": " << ex << std::endl;
            }
        }
    }
    std::cout << "-----------------------------------------------------------------------------------------------------------------\n";

    // Wyœwietlenie w wymaganym, skondensowanym formacie
    if (best_opt.flag == 1) {
        double t0 = best_opt.x(0);
        double t1 = best_opt.x(1);
        double t2 = best_opt.x(2);
        double cost = m2d(best_opt.y);

        // Wymagany format: theta0_*, theta1_*, theta2_*, J(theta_*), P(theta_*), g_calls
        std::cout << "\n\n" << std::fixed << std::setprecision(8)
                  << t0 << ", "
                  << t1 << ", "
                  << t2 << ", "
                  << cost << ", "
                  << best_accuracy << ", "
                  << best_g_calls << "\n";

        // Dodatkowy wydruk granicy klasyfikacji (dla sprawdzenia)
        std::cout << "\nNajlepszy przypadek: " << best_case_info << ", Dok³adnoœæ: " << best_accuracy << "%\n";
        std::cout << "Równanie granicy klasyfikacji (do naniesienia na wykres, h_theta(x) = 0.5):\n";

        if (std::abs(t2) > 1e-8) {
             std::cout << "x2 = - (" << t0 << " + " << t1 << " * x1) / " << t2 << "\n";
        } else {
             std::cout << "x2 = (NIEMO¯LIWE DO OBLICZENIA) Theta2 jest bliskie zeru. \n";
        }
    }
}

void lab5()
{

}

void lab6()
{

}
