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
void lab2_rosen(const std::vector<matrix>& starting_points, int num_trials, double x_min, double x_max);
void lab2_HJ(const std::vector<matrix>& starting_points, int num_trials, double x_min, double x_max);
void lab2_robot();
void lab3();
void lab4();
void lab5();
void lab6();


// Funkcja do generowania i zapisywania punktów startowych
std::vector<matrix> generate_starting_points(int count, double x_min, double x_max) {
    std::vector<matrix> points;
    for (int i = 0; i < count; ++i) {
        matrix x0 = rand_mat(2);
        for (int j = 0; j < 2; ++j) {
            x0(j) = x_min + (x_max - x_min) * x0(j);
        }
        points.push_back(x0);
    }
    return points;
}

int main()
{
    int num_trials = 100;
    double x_min = -1.0;
    double x_max = 1.0;
    std::vector<matrix> starting_points = generate_starting_points(num_trials, x_min, x_max);

    try
    {
        lab2_rosen(starting_points, num_trials, x_min, x_max);
        lab2_HJ(starting_points, num_trials, x_min, x_max);
        //lab2_robot();
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




void lab2_rosen(const std::vector<matrix>& starting_points, int num_trials, double x_min, double x_max) {
    try {
        std::cout << "\n=== Lab2 Rosenbrock - 3x100 optymalizacji ===" << std::endl;

        int N = 2;
        double alpha = 2.0;
        double beta = 0.5;
        double epsilon = 1e-6;
        int Nmax = 10000;
        matrix ud1(NAN), ud2(NAN);

        double step_sizes[3][2] = { {0.1, 0.1}, {0.5, 0.5}, {1.0, 1.0} };
        int num_steps = 3;

        ofstream results_file("wyniki_rosen_3x100.csv");
        results_file << "Nr,s1,s2,x1_start,x2_start,x1_opt,x2_opt,f_opt,f_calls,flag\n";

        std::cout << "\nPrzeprowadzanie 3x100 optymalizacji dla różnych długości kroku..." << std::endl;

        for (int step_idx = 0; step_idx < num_steps; ++step_idx) {
            double s0_data[2] = { step_sizes[step_idx][0], step_sizes[step_idx][1] };
            matrix s0(2, s0_data);

            std::cout << "\n--- Długości kroku s = [" << s0(0) << ", " << s0(1) << "] ---" << std::endl;

            int success_count = 0;
            int global_min_count = 0;
            double sum_f_calls = 0.0;
            double sum_f_opt = 0.0;

            for (int trial = 0; trial < num_trials; ++trial) {
                matrix x0 = starting_points[trial];

                solution::clear_calls();
                solution opt = Rosen(ff2T, x0, s0, alpha, beta, epsilon, Nmax, ud1, ud2);

                results_file << (trial + 1) << "," << s0(0) << "," << s0(1) << ","
                    << x0(0) << "," << x0(1) << ","
                    << opt.x(0) << "," << opt.x(1) << ","
                    << opt.y(0) << "," << solution::f_calls << ","
                    << opt.flag << "\n";

                if (opt.y(0) < 1e-3 && (opt.flag == 0 || opt.flag == 1 || opt.flag == 2)) {
                    global_min_count++;
                    sum_f_calls += solution::f_calls;
                    sum_f_opt += opt.y(0);
                }

                if (opt.flag == 0 || opt.flag == 1 || opt.flag == 2) {
                    success_count++;
                }

                if ((trial + 1) % 20 == 0)
                    std::cout << "  Ukończono " << (trial + 1) << "/" << num_trials << " optymalizacji" << std::endl;
            }

            std::cout << "Wyniki dla s = [" << s0(0) << ", " << s0(1) << "]:" << std::endl;
            std::cout << "  Udane optymalizacje: " << success_count << "/" << num_trials << std::endl;
            std::cout << "  Minimum globalne: " << global_min_count << "/" << num_trials << std::endl;
            if (global_min_count > 0) {
                std::cout << "  Średnia liczba wywołań funkcji celu: " << (sum_f_calls / global_min_count) << std::endl;
                std::cout << "  Średnia wartość funkcji celu: " << (sum_f_opt / global_min_count) << std::endl;
            }
        }

        results_file.close();
        std::cout << "\nWyniki zapisane do: wyniki_rosen_3x100.csv" << std::endl;

    }
    catch (string ex_info) {
        std::cerr << "Wystąpił błąd w lab2_rosen():\n" << ex_info << "\n";
    }
}

void lab2_HJ(const std::vector<matrix>& starting_points, int num_trials, double x_min, double x_max) {
    std::cout << "\n=== Lab2 HJ - Testowa funkcja celu (100 optymalizacji) ===" << std::endl;

    matrix ud1, ud2;
    double epsilon = 1e-4;
    int Nmax = 50000;

    double step_sizes[] = { 0.1, 0.5, 1.0 };
    int num_steps = 3;

    ofstream results_file("wyniki_100_optymalizacji.csv");
    results_file << "Nr,s,x1_start,x2_start,x1_opt,x2_opt,f_opt,f_calls,flag\n";

    std::cout << "\nPrzeprowadzanie 100 optymalizacji dla kazdej dlugosci kroku..." << std::endl;

    for (int step_idx = 0; step_idx < num_steps; ++step_idx) {
        double s = step_sizes[step_idx];
        double alpha = 0.5;

        std::cout << "\n--- Dlugosc kroku s = " << s << " ---" << std::endl;

        int success_count = 0;
        double sum_f_calls = 0.0;

        for (int trial = 0; trial < num_trials; ++trial) {
            matrix x0 = starting_points[trial];

            solution::clear_calls();
            solution opt = HJ(ff2T, x0, s, alpha, epsilon, Nmax, ud1, ud2);

            results_file << (trial + 1) << "," << s << ","
                << x0(0) << "," << x0(1) << ","
                << opt.x(0) << "," << opt.x(1) << ","
                << opt.y(0) << "," << solution::f_calls << ","
                << opt.flag << "\n";

            if (opt.y(0) < 0.01 && opt.flag == 1) {
                success_count++;
                sum_f_calls += solution::f_calls;
            }

            if ((trial + 1) % 20 == 0)
                std::cout << "  Ukończono " << (trial + 1) << "/" << num_trials << " optymalizacji" << std::endl;
        }

        std::cout << "Wyniki dla s = " << s << ":" << std::endl;
        std::cout << "  Sukces: " << success_count << "/" << num_trials << std::endl;
        if (success_count > 0) {
            std::cout << "  Srednia liczba wywolan funkcji celu: " << (sum_f_calls / success_count) << std::endl;
        }
    }

    results_file.close();
    std::cout << "\nWyniki zapisane do: wyniki_100_optymalizacji.csv" << std::endl;
}

void lab2_robot()
{
	std::cout << "\n=== Lab2 Robot - Optymalizacja ramienia robota ===" << std::endl;
	
	// Test poprawności implementacji dla k1=5, k2=5
	std::cout << "\n--- Test poprawnosci implementacji ---" << std::endl;
	double test_k[] = { 5.0, 5.0 };
	matrix x_test(2, test_k);
	
	// Wartości referencyjne
	double ref_data[] = { M_PI, 0.0 }; // alpha_ref = π, omega_ref = 0
	matrix ud1(2, ref_data);
	matrix ud2; // nie używane w ff2R
	
	matrix Q_test = ff2R(x_test, ud1, ud2);
	std::cout << "Q(k1=5, k2=5) = " << Q_test(0) << std::endl;
	std::cout << "Oczekiwana wartosc: ~775.229" << std::endl;
	
	// Optymalizacja metodą Hooke-Jeeves - pojedyncza
	std::cout << "\n--- Optymalizacja metoda Hooke-Jeeves (pojedyncza) ---" << std::endl;
	
	solution::clear_calls();
	
	// Punkt startowy z przedziału [0,20] x [0,20]
	double x0_data[] = { 10.0, 10.0 }; // k1, k2
	matrix x0(2, x0_data);
	
	double s = 1.0;        // długość kroku
	double alpha = 0.5;    // współczynnik zmniejszania kroku
	double epsilon = 1e-4; // dokładność
	int Nmax = 50000;      // maksymalna liczba wywołań funkcji celu
	
	std::cout << "Punkt startowy: k1 = " << x0(0) << ", k2 = " << x0(1) << std::endl;
	std::cout << "Parametry: s = " << s << ", alpha = " << alpha << ", epsilon = " << epsilon << std::endl;
	
	solution opt = HJ(ff2R, x0, s, alpha, epsilon, Nmax, ud1, ud2);
	
	std::cout << "\n--- Wyniki optymalizacji ---" << std::endl;
	std::cout << "Optymalne wspolczynniki wzmocnienia:" << std::endl;
	std::cout << "  k1 = " << opt.x(0) << " Nm" << std::endl;
	std::cout << "  k2 = " << opt.x(1) << " Nms" << std::endl;
	std::cout << "Wartosc funkcjonalu jakosci: Q = " << opt.y(0) << std::endl;
	std::cout << "Liczba wywolan funkcji celu: " << solution::f_calls << std::endl;
	std::cout << "Flaga: " << opt.flag << " (1 = sukces, 0 = przekroczono Nmax)" << std::endl;
	
	// Symulacja dla optymalnych parametrów
	std::cout << "\n--- Generowanie symulacji dla optymalnych parametrow ---" << std::endl;
	
	matrix Y0(2, 1);
	Y0(0) = 0.0;  // alpha(0) = 0
	Y0(1) = 0.0;  // omega(0) = 0
	
	matrix* Y = solve_ode(df2R, 0.0, 0.1, 100.0, Y0, ud1, opt.x);
	
	ofstream Sout("symulacja_robot.csv");
	Sout << "t,alpha,omega\n";
	int n = get_len(Y[0]);
	for (int i = 0; i < n; ++i)
	{
		Sout << Y[0](i) << "," << Y[1](i, 0) << "," << Y[1](i, 1) << "\n";
	}
	Sout.close();
	std::cout << "Wyniki symulacji zapisane do: symulacja_robot.csv" << std::endl;
	
	Y[0].~matrix();
	Y[1].~matrix();
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
