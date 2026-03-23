#include "Integration.h"
#include <iostream>
#include <vector>

namespace Integration {

    double evaluatePolynomial(const std::vector<double>& coefficients, double x) {
        double wynik = 0.0;
        for (int i = coefficients.size() - 1; i >= 0; --i) {
            wynik = wynik * x + coefficients[i];
        }
        return wynik;
    }

    double rectangleRule(
        const std::vector<double>& coefficients,
        double a, double b, int n) {

        if (n <= 0) return 0.0;
        double h = (b - a) / n;
        double suma = 0.0;
        for (int i = 0; i < n; ++i) {
            double x = a + i * h;
            suma += evaluatePolynomial(coefficients, x);
        }
        return suma * h;
    }

    double trapezoidalRule(
        const std::vector<double>& coefficients,
        double a, double b, int n) {

        double h = (b - a) / n;
        double suma = 0.0;
        double pdst_a = evaluatePolynomial(coefficients, a);
        double pdst_b;
        for (int i = 1; i <= n; ++i) {
            pdst_b = evaluatePolynomial(coefficients, a + i * h);
            suma += (pdst_a + pdst_b);
            pdst_a = pdst_b;
        }
        return suma * h * 0.5;
    }

    double simpsonsRule(
        const std::vector<double>& coefficients,
        double a, double b, int n) {

        if (n % 2 != 0) return 0.0;
        if (n % 2 != 0) {
            throw std::invalid_argument("Liczba przedzia³ów n musi byæ parzysta dla metody Simpsona!");
        }
        double h = (b - a) / n;
        double suma = evaluatePolynomial(coefficients, a) + evaluatePolynomial(coefficients, b);
        for (int i = 1; i < n; ++i) {
            double x = a + i * h;
            if (i % 2 == 0) {
                suma += 2 * evaluatePolynomial(coefficients, x);
            }
            else {
                suma += 4 * evaluatePolynomial(coefficients, x);
            }
        }
        return suma * h / 3.0;
    }

    double rectangleRule(
        std::function<double(double)> func,
        double a, double b, int n) {

        if (n <= 0) return 0.0;
        double h = (b - a) / n;
        double suma = 0.0;
        for (int i = 0; i < n; ++i) {
            double x = a + i * h;
            suma += func(x);
        }
        return suma * h;
    }

    double trapezoidalRule(
        std::function<double(double)> func,
        double a, double b, int n) {

        double h = (b - a) / n;
        double suma = 0.0;
        double pdst_a = func(a);
        double pdst_b;
        for (int i = 1; i <= n; ++i) {
            double x = a + i * h;
            pdst_b = func(x);
            suma += (pdst_a + pdst_b);
            pdst_a = pdst_b;
        }
        return suma * h * 0.5;
    }

    double simpsonsRule(
        std::function<double(double)> func,
        double a, double b, int n) {

        if (n % 2 != 0) {
            return 0.0;
        }
        double h = (b - a) / n;
        double suma = func(a) + func(b);
        for (int i = 1; i < n; ++i) {
            double x = a + i * h;
            if (i % 2 == 0) {
                suma += 2 * func(x);
            }
            else {
                suma += 4 * func(x);
            }
        }
        return suma * h / 3.0;
    }

    GaussLegendreData getGaussLegendreData(int n) {
        if (n == 2) {
            return { {-0.5773502692, 0.5773502692}, {1.0, 1.0} };
        }
        else if (n == 3) {
            return { {-0.7745966692, 0.0, 0.7745966692}, {0.5555555556, 0.8888888889, 0.5555555556} };
        }
        else if (n == 4) {
            return { {-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116},
                    {0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451} };
        }
        else if (n == 5) {
            return { {-0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459},
                    {0.2369268856, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268856} };
        }
        else {
            throw std::invalid_argument("Nieobs³ugiwana liczba wêz³ów");
        }
    }

    double gaussLegendreIntegration(
        std::function<double(double)> func,
        double a, double b, int n, int m) {

        if (n <= 0) return 0.0;
        GaussLegendreData gl = getGaussLegendreData(n);
        double wynik = 0.0;
        double h = (b - a) / m;

        for (int j = 0; j < m; j++) {
            double aj = a + j * h;
            double bj = aj + h;

            for (int i = 0; i < gl.nodes.size(); i++) {
                double xi = 0.5 * (bj - aj) * gl.nodes[i] + 0.5 * (bj + aj);
                wynik += gl.weights[i] * func(xi);
            }
        }

        wynik *= 0.5 * h;
        return wynik;
    }
}