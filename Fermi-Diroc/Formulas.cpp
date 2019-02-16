#include "stdafx.h"
#include "Formulas.h"
#include <cmath>

using namespace std;

// get Fermat value
double f(double E, double T, double miu) {
	double beta = 1 / T;
	return 1 / (exp(beta * (E - miu)) + 1);
}

// function that should be integrated for p_i
double p_iFunctionInsideIntegral(double E, double m, double T, double miu) {
	return pow(E * E - m * m, 1.5) * f(E, T, miu);
}

// function that should be integrated for n_i
double n_iFunctionInsideIntegral(double E, double m, double T, double miu) {
	return sqrt(E * E - m * m) * E * f(E, T, miu);
}

// function that should be integrated for e_i
double e_iFunctionInsideIntegral(double E, double m, double T, double miu) {
	return sqrt(E * E - m * m) * E * E * f(E, T, miu);
}

// Get p_i value
double compute_p_i(double Emax, double m, double T, double miu) {
	double before_intrgration = (FormulaConstants::g / (6 * FormulaConstants::M_PI * FormulaConstants::M_PI));

	double Emin = m;
	double deltaE = Emax - Emin / FormulaConstants::N - 1;
	double n = 0;
	double E = m;
	for (int j = 0; j <= FormulaConstants::N - 1; j++) {
		E += deltaE;
		double n_j = p_iFunctionInsideIntegral(E, m, T, miu) * deltaE;
		n += n_j;
	}
	return n * before_intrgration;
}

// Get n_i value
double compute_n_i(double Emax, double m, double T, double miu) {
	double before_intrgration = (FormulaConstants::g / (2 * FormulaConstants::M_PI * FormulaConstants::M_PI));

	double Emin = m;
	double deltaE = Emax - Emin / FormulaConstants::N - 1;
	double n = 0;
	double E = m;
	for (int j = 0; j <= FormulaConstants::N - 1; j++) {
		E += deltaE;
		double n_j = n_iFunctionInsideIntegral(E, m, T, miu) * deltaE;
		n += n_j;
	}
	return n * before_intrgration;
}

// Get e_i value
double compute_e_i(double Emax, double m, double T, double miu) {
	double before_intrgration = (FormulaConstants::g / (2 * FormulaConstants::M_PI * FormulaConstants::M_PI));

	double Emin = m;
	double deltaE = Emax - Emin / FormulaConstants::N - 1;
	double n = 0;
	double E = m;
	for (int j = 0; j <= FormulaConstants::N - 1; j++) {
		E += deltaE;
		double n_j = e_iFunctionInsideIntegral(E, m, T, miu) * deltaE;
		n += n_j;
	}
	return n * before_intrgration;
}