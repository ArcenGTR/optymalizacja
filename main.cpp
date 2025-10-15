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
		lab1();
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

	/*

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

    // POPRAWKA: Pe³ne zwalnianie pamiêci, tak jak w lab0
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
    int Nmax = 500; // Zmniejszono Nmax, bo symulacje s¹ kosztowne

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

	*/

	//Symulacja z optymalnymi wynikami

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

}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
