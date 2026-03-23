#ifndef NONLINEAR_ALGEBRA_H
#define NONLINEAR_ALGEBRA_H

#include <vector>
#include <functional>
#include <string>

namespace NonlinearEquations {

    using Function = double(*)(double);
    using Derivative = double(*)(double);

    // Metody numeryczne
    double bisection(Function f, double a, double b, double eps, int max_iter,
        bool log_steps = false, std::ofstream* out = nullptr, const char* fname = "");

    double newton(Function f, Derivative df, double x0, double eps, int max_iter,
        bool log_steps = false, std::ofstream* out = nullptr, const char* fname = "");

    double secant(Function f, double x0, double x1, double eps, int max_iter,
        bool log_steps = false, std::ofstream* out = nullptr, const char* fname = "");

    std::vector<std::pair<double, double>> findIntervals(Function f, double a, double b, double step);

    void addUnique(std::vector<double>& roots, double x, double eps = 1e-7);
    double minAbsError(const std::vector<double>& reference, double x);
}
#endif // NONLINEAR_ALGEBRA_H
