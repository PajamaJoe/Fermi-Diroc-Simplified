#pragma once

// constants
namespace FormulaConstants {
	const double M_PI = 3.1415926535897;
	const int N = 3000;
	const int g = 2;
}

double f(double E, double T, double miu);

double compute_p_i(double Emax, double m, double T, double miu);
double compute_n_i(double Emax, double m, double T, double miu);
double compute_e_i(double Emax, double m, double T, double miu);