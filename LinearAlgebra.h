#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include <iostream>
#include <vector>
#include <functional>
#include <string>

// ================== ALGEBRA LINIOWA ==================
namespace LinearAlgebra {

	// Operacje na macierzach
	std::vector<std::vector<double>> createIdentityMatrix(int n);
	std::vector<std::vector<double>> multiplyMatrices(
		const std::vector<std::vector<double>>& A,
		const std::vector<std::vector<double>>& B
	);
	bool areMatricesEqual(
		const std::vector<std::vector<double>>& A,
		const std::vector<std::vector<double>>& B,
		double eps = 1e-6
	);

	// Rozk³ad LU
	void luDecomposition(
		const std::vector<std::vector<double>>& A_orig,
		std::vector<std::vector<double>>& L,
		std::vector<std::vector<double>>& U,
		std::vector<std::vector<double>>& P
	);

	// Rozwi¹zywanie uk³adów równañ
	std::vector<double> solveLowerTriangular(
		const std::vector<std::vector<double>>& L,
		const std::vector<double>& b,
		const std::vector<std::vector<double>>& P
	);
	std::vector<double> solveUpperTriangular(
		const std::vector<std::vector<double>>& U,
		const std::vector<double>& y
	);
	std::vector<double> solveLU(
		std::vector<std::vector<double>> A,
		std::vector<double> b
	);

	// Eliminacja Gaussa
	void gaussianElimination(
		int N,
		std::vector<double>& b,
		std::vector<std::vector<double>>& A,
		bool verbose = false
	);
	std::vector<double> backSubstitution(
		int N,
		const std::vector<std::vector<double>>& A,
		const std::vector<double>& b
	);
	bool verifySystem(
		int N,
		const std::vector<std::vector<double>>& A,
		const std::vector<double>& x,
		const std::vector<double>& b,
		double epsilon = 1e-6
	);
	// Obliczanie ilorazów ró¿nicowych (interpolacja Newtona)
	std::vector<double> dividedDifferences(
		const std::vector<double>& x_values,
		const std::vector<double>& y_values
	);

	// Wartoœæ wielomianu Newtona w punkcie
	double evaluateNewton(
		const std::vector<double>& coefficients,
		const std::vector<double>& x_values,
		double x
	);

	// Wartoœæ wielomianu w postaci naturalnej (pomocnicze)
	double evaluateNatural(const std::vector<double>& coefficients, double x);

	// Wartoœæ wielomianu metod¹ Hornera (pomocnicze)
	double evaluateHorner(const std::vector<double>& coefficients, double x);
}
#endif // LINEAR_ALGEBRA_H#pragma once
