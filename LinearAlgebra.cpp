#include "LinearAlgebra.h"
#include <iostream>
#include <cmath>

// ================== ALGEBRA LINIOWA ==================
namespace LinearAlgebra {

	std::vector<std::vector<double>> createIdentityMatrix(int n) {
		std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
		for (int i = 0; i < n; ++i)
			I[i][i] = 1.0;
		return I;
	}

	std::vector<std::vector<double>> multiplyMatrices(
		const std::vector<std::vector<double>>& A,
		const std::vector<std::vector<double>>& B) {

		if (A.empty() || B.empty() || A[0].size() != B.size())
			throw std::invalid_argument("Nieprawid³owe rozmiary macierzy do mno¿enia");

		int n = A.size();
		int m = B[0].size();
		int k = B.size();
		std::vector<std::vector<double>> result(n, std::vector<double>(m, 0.0));
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				for (int l = 0; l < k; l++)
					result[i][j] += A[i][l] * B[l][j];
		return result;
	}

	bool areMatricesEqual(
		const std::vector<std::vector<double>>& A,
		const std::vector<std::vector<double>>& B,
		double eps) {

		int n = A.size();
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				if (std::abs(A[i][j] - B[i][j]) > eps)
					return false;
		return true;
	}

	void luDecomposition(
		const std::vector<std::vector<double>>& A_orig,
		std::vector<std::vector<double>>& L,
		std::vector<std::vector<double>>& U,
		std::vector<std::vector<double>>& P) {

		int n = A_orig.size();
		std::vector<std::vector<double>> A = A_orig;

		L = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
		U = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
		P = createIdentityMatrix(n);

		for (int k = 0; k < n; k++) {
			double max_val = 0.0;
			int max_idx = k;

			for (int i = k; i < n; i++) {
				if (std::abs(A[i][k]) > max_val) {
					max_val = std::abs(A[i][k]);
					max_idx = i;
				}
			}

			if (max_idx != k) {
				std::swap(A[k], A[max_idx]);
				std::swap(P[k], P[max_idx]);
				for (int j = 0; j < k; j++) {
					std::swap(L[k][j], L[max_idx][j]);
				}
			}

			L[k][k] = 1.0;

			for (int i = k + 1; i < n; i++) {
				L[i][k] = A[i][k] / A[k][k];
			}

			for (int j = k; j < n; j++) {
				U[k][j] = A[k][j];
			}

			for (int i = k + 1; i < n; i++) {
				for (int j = k + 1; j < n; j++) {
					A[i][j] -= L[i][k] * U[k][j];
				}
			}
		}
	}

	std::vector<double> solveLowerTriangular(
		const std::vector<std::vector<double>>& L,
		const std::vector<double>& b,
		const std::vector<std::vector<double>>& P) {

		int n = L.size();
		std::vector<double> y(n, 0.0);
		std::vector<double> Pb(n, 0.0);

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				Pb[i] += P[i][j] * b[j];

		for (int i = 0; i < n; i++) {
			double sum = 0.0;
			for (int j = 0; j < i; j++) {
				sum += L[i][j] * y[j];
			}
			y[i] = Pb[i] - sum;
		}
		return y;
	}

	std::vector<double> solveUpperTriangular(
		const std::vector<std::vector<double>>& U,
		const std::vector<double>& y) {

		int n = U.size();
		std::vector<double> x(n, 0.0);
		for (int i = n - 1; i >= 0; i--) {
			double sum = 0.0;
			for (int j = i + 1; j < n; j++) {
				sum += U[i][j] * x[j];
			}
			x[i] = (y[i] - sum) / U[i][i];
		}
		return x;
	}

	std::vector<double> solveLU(
		std::vector<std::vector<double>> A,
		std::vector<double> b) {

		int n = A.size();
		if (n == 0 || A[0].size() != n || b.size() != n) {
			return {};
		}

		std::vector<std::vector<double>> L, U, P;
		luDecomposition(A, L, U, P);
		std::vector<double> z = solveLowerTriangular(L, b, P);
		std::vector<double> x = solveUpperTriangular(U, z);
		return x;
	}

	void gaussianElimination(
		int N,
		std::vector<double>& b,
		std::vector<std::vector<double>>& A,
		bool verbose) {

		const double EPSILON = 1e-9;
		for (int i = 0; i < N; ++i) {
			int maxRow = i;
			for (int k = i + 1; k < N; ++k) {
				if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
					maxRow = k;
				}
			}

			if (std::abs(A[maxRow][i]) < EPSILON) {
				if (verbose) {
					std::cout << "Uk³ad równañ jest sprzeczny lub ma nieskoñczenie wiele rozwi¹zañ." << std::endl;
				}
				return;
			}

			std::swap(A[maxRow], A[i]);
			std::swap(b[maxRow], b[i]);

			for (int k = i + 1; k < N; ++k) {
				double c = -A[k][i] / A[i][i];
				for (int j = i; j < N; ++j) {
					A[k][j] += c * A[i][j];
				}
				b[k] += c * b[i];
			}
		}
	}

	std::vector<double> backSubstitution(
		int N,
		const std::vector<std::vector<double>>& A,
		const std::vector<double>& b) {

		std::vector<double> x(N);
		for (int i = N - 1; i >= 0; --i) {
			x[i] = b[i];
			for (int j = i + 1; j < N; ++j) {
				x[i] -= A[i][j] * x[j];
			}
			x[i] /= A[i][i];
		}
		return x;
	}

	bool verifySystem(
		int N,
		const std::vector<std::vector<double>>& A,
		const std::vector<double>& x,
		const std::vector<double>& b,
		double epsilon) {

		for (int i = 0; i < N; ++i) {
			double suma = 0;
			for (int j = 0; j < N; ++j) {
				suma += A[i][j] * x[j];
			}
			if (std::abs(suma - b[i]) > epsilon) {
				return false;
			}
		}
		return true;
	}
}