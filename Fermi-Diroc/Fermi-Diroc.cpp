// Fermi-Diroc.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Formulas.h"
#include <iostream>

using namespace std;

int main()
{
	// input
	double E, T, miu, m;

	double p_i;
	double n_i;
	double e_i;

	cout << "E = ";
	cin >> E;

	cout << "T = ";
	cin >> T;

	cout << "miu = ";
	cin >> miu;

	cout << "m = ";
	cin >> m;

	// output
	double fermi = f(E, T, miu);

	p_i = compute_p_i(E, m, T, miu);
	n_i = compute_n_i(E, m, T, miu);
	e_i = compute_e_i(E, m, T, miu);

	cout << "Fermi = " << fermi << endl;
	cout << "p_i = " << p_i << endl;
	cout << "n_i = " << n_i << endl;
	cout << "e_i = " << e_i << endl;

	return 0;
}

