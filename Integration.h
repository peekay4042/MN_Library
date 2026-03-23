#ifndef INTGRATION_H
#define INTGRATION_H

#include <vector>
#include <functional>
#include <string>

// ================== ALGEBRA LINIOWA ==================
namespace Integration {

    // Metody podstawowe
    double rectangleRule(
        const std::vector<double>& coefficients,
        double a, double b, int n
    );
    double trapezoidalRule(
        const std::vector<double>& coefficients,
        double a, double b, int n
    );
    double simpsonsRule(
        const std::vector<double>& coefficients,
        double a, double b, int n
    );

    // Metody dla funkcji
    double rectangleRule(
        std::function<double(double)> func,
        double a, double b, int n
    );
    double trapezoidalRule(
        std::function<double(double)> func,
        double a, double b, int n
    );
    double simpsonsRule(
        std::function<double(double)> func,
        double a, double b, int n
    );

    // Kwadratura Gaussa-Legendre'a
    struct GaussLegendreData {
        std::vector<double> nodes;
        std::vector<double> weights;
    };

    GaussLegendreData getGaussLegendreData(int n);
    double gaussLegendreIntegration(
        std::function<double(double)> func,
        double a, double b, int n, int m
    );

    // Pomocnicze
    double evaluatePolynomial(const std::vector<double>& coefficients, double x);
}
#endif // INTGRATION_H