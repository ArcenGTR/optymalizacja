/*********************************************
Kod stanowi uzupełnienie materiałów do ćwiczeń
w ramach przedmiotu metody optymalizacji.
Kod udostępniony na licencji CC BY-SA 3.0
Autor: dr inż. Łukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include"opt_alg.h"
#include <cmath> // Add this include directive at the top of the file for M_PI
#ifndef M_PI
#define M_PI 3.14159265358979323846 // Define M_PI if not already defined
#endif

void lab0();
void lab1();
void lab2();
void lab2_HJ();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		//lab2();
		lab2_HJ();
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
	double epsilon = 1e-2;									// dokładność
	int Nmax = 10000;										// maksymalna liczba wywołań funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz górne ograniczenie
		a(2, 1);											// dokładne rozwiązanie optymalne
	solution opt;											// rozwiązanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywołanie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Wahadlo
	Nmax = 1000;											// dokładność
	epsilon = 1e-2;											// maksymalna liczba wywołań funkcji celu
	lb = 0, ub = 5;											// dolne oraz górne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahadła
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywołanie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki początkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment siły działający na wahadło oraz czas działania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwiązujemy równanie różniczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumień do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumień
	Y[0].~matrix();											// usuwamy z pamięci rozwiązanie RR
	Y[1].~matrix();
}


//matrix ff(matrix x, matrix ud1, matrix ud2)
//{
//	//double val = -cos(0.1 * xval) * exp(-pow(0.1 * xval - 2 * M_PI, 2)) + 0.002 * pow(0.1 * xval, 2);
//	//return matrix(1,1,val)
//	return matrix(1, 1, 0);
//}

void lab1()
{
		//std::cout << "start lab1 " <<std::endl;
		//matrix ud1, ud2;

		//double x0 = 1.0;
		//double y0 = 1.0;
		//double d = 0.1;
		//double alpha = 1.1;
		//int Nmax = 10000;

		//double* wynik_ekspansji = expansion(ff, x0, d, alpha, Nmax, ud1, ud2);
		//std::cout << "Wynik algorytmu ekspansji: " << wynik_ekspansji[0] << ", " << wynik_ekspansji[1] << std::endl;

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

		double s0_data[] = { 0.1, 0.1 }; // Poczatkowe kroki
		matrix s0(N, s0_data);

		double alpha = 3.0;     // Wspólczynnik ekspansji
		double beta = 0.5;      // Wspólczynnik kontrakcji
		double epsilon = 1e-4;  // Dokladnosc
		int Nmax = 1000;        // Maksymalna liczba wywolan funkcji celu

		matrix ud1(NAN), ud2(NAN);

		std::cout << "--- Rozpoczecie optymalizacji metoda Rosenbrocka ---\n";
		std::cout << "Punkt startowy x0: " << x0 << "\n";

		solution Xopt = Rosen(ff2T, x0, s0, alpha, beta, epsilon, Nmax, ud1, ud2);

		std::cout << "\n--- Wyniki optymalizacji ---\n";
		std::cout << Xopt;

		std::cout << "\nKomentarz do flagi wyjscia:\n";
		if (Xopt.flag == 2) {
			std::cout << "Flaga 2: Algorytm zbiegl do rozwiazania z zadana dokladnoscia (max(|s_j|) < epsilon).\n";
		}
		else if (Xopt.flag == 0) {
			std::cout << "Flaga 0: Przekroczono maksymalna liczbe wywolan funkcji celu (Nmax).\n";
		}
		else {
			std::cout << "Flaga " << Xopt.flag << ": Nieznany status wyjscia.\n";
		}

	}
	catch (string ex_info)
	{
		std::cerr << "Wystapil blad w run_rosen_test():\n" << ex_info << "\n";
	}
}

void lab2_HJ()
{
	std::cout << "start lab2_HJ (Hooke-Jeeves)" << std::endl;

	// Punkt startowy dla funkcji 2D ff2T
	double x0_data[] = { 0.5, 0.5 }; // punkt startowy
	matrix x0(2, x0_data);
	
	double s = 0.5;        // długość kroku
	double alpha = 0.5;    // współczynnik zmniejszania kroku
	double epsilon = 1e-3; // dokładność
	int Nmax = 10000;      // maksymalna liczba wywołań funkcji celu
	matrix ud1, ud2;

	solution opt = HJ(ff2T, x0, s, alpha, epsilon, Nmax, ud1, ud2);

	std::cout << "Wynik HJ: x = [" << opt.x(0) << ", " << opt.x(1) << "], f(x) = " << opt.y(0) << ", flag = " << opt.flag << std::endl;
	std::cout << "Liczba wywolan funkcji celu: " << solution::f_calls << std::endl;
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
