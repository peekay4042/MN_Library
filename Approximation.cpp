#include "Approximation.h"
#include <iostream>
// ================== APROKSYMACJA ==================
namespace Approximation {

    double polynomialBasis(int n, double x) {
        if (n < 0)
            throw std::invalid_argument("Stopieñ wielomianu nie mo¿e byæ ujemny");
        return std::pow(x, n);
    }

    void computeSystemMatrix(
        int degree,
        double a, double b,
        std::vector<std::vector<double>>& A,
        std::vector<double>& b_vec,
        std::function<double(double)> func) {

        auto integral = [&](std::function<double(double)> integrand) {
            return Integration::gaussLegendreIntegration(integrand, a, b, 5, 1000);
            };

        for (int i = 0; i <= degree; ++i) {
            for (int j = 0; j <= degree; ++j) {
                A[i][j] = integral([=](double x) {
                    return polynomialBasis(i, x) * polynomialBasis(j, x);
                    });
            }
            b_vec[i] = integral([=](double x) {
                return polynomialBasis(i, x) * func(x);
                });
        }
    }

    double evaluateApproximation(
        const std::vector<double>& coefficients,
        double x) {

        double wynik = 0.0;
        for (int i = 0; i < coefficients.size(); ++i) {
            wynik += coefficients[i] * polynomialBasis(i, x);
        }
        return wynik;
    }

    double calculateAverageError(
        const std::vector<double>& coeffs,
        std::function<double(double)> func,
        double a, double b, int numSamples)
    {
        if (numSamples <= 0)
            return 0.0;
        double sum = 0.0;
        for (int i = 0; i < numSamples; ++i) {
            double x = a + (b - a) * i / (numSamples - 1);
            double diff = func(x) - evaluateApproximation(coeffs, x);
            sum += std::abs(diff);
        }
        return sum / numSamples;
    }
}