#ifndef APROXIMATION_H
#define APROXIMATION_H
#include "Integration.h"
#include <vector>
#include <functional>
#include <cmath>


// ================== APROKSYMACJA ==================
namespace Approximation {

    // Funkcje bazowe i aproksymowane
    double polynomialBasis(int n, double x);

    // Aproksymacja metod¹ najmniejszych kwadratów
    void computeSystemMatrix(
        int degree,
        double a, double b,
        std::vector<std::vector<double>>& A,
        std::vector<double>& b_vec,
        std::function<double(double)> func
    );

    double evaluateApproximation(
        const std::vector<double>& coefficients,
        double x
    );

    double calculateAverageError(
        const std::vector<double>& coefficients,
        std::function<double(double)> func,
        double a, double b,
        int num_points = 16
    );
}
#endif
