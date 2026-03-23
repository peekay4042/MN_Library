#include "NonlinearAlgebra.h"
#include <iostream>
#include <cmath>
#include <fstream>

namespace NonlinearEquations {

    using Function = double(*)(double);
    using Derivative = double(*)(double);

    double bisection(Function f, double a, double b, double eps, int max_iter,
        bool log_steps, std::ofstream* out, const char* fname) {
        if (f(a) * f(b) >= 0) return NAN;
        double c;
        for (int i = 0; i < max_iter; ++i) {
            c = (a + b) / 2.0;
            if (log_steps && out)
                *out << fname << ",Bisekcja," << (i + 1) << "," << a << "," << b << "," << c << "\n";

            if (std::fabs(f(c)) < eps || std::fabs(b - a) < eps) return c;
            if (f(a) * f(c) < 0) b = c;
            else a = c;
        }
        return c;
    }

    double newton(Function f, Derivative df, double x0, double eps, int max_iter,
        bool log_steps, std::ofstream* out, const char* fname) {
        double x = x0;
        for (int i = 0; i < max_iter; ++i) {
            double fx = f(x);
            double dfx = df(x);
            if (log_steps && out)
                *out << fname << ",Newton," << (i + 1) << "," << x0 << ",," << x << "\n";

            if (std::fabs(dfx) < 1e-12) return NAN;
            double x1 = x - fx / dfx;
            if (std::fabs(x1 - x) < eps) return x1;
            x = x1;
        }
        return x;
    }

    double secant(Function f, double x0, double x1, double eps, int max_iter,
        bool log_steps, std::ofstream* out, const char* fname) {
        double f0 = f(x0), f1 = f(x1);
        for (int i = 0; i < max_iter; ++i) {
            if (log_steps && out)
                *out << fname << ",Sieczne," << (i + 1) << "," << x0 << "," << x1 << "," << x1 << "\n";

            if (std::fabs(f1 - f0) < 1e-12) return NAN;
            double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
            if (std::fabs(x2 - x1) < eps) return x2;
            x0 = x1; f0 = f1; x1 = x2; f1 = f(x1);
        }
        return x1;
    }

    std::vector<std::pair<double, double>> findIntervals(Function f, double a, double b, double step) {
        std::vector<std::pair<double, double>> intervals;
        double x1 = a, x2 = a + step;
        while (x2 <= b) {
            double y1 = f(x1), y2 = f(x2);
            if (std::isfinite(y1) && std::isfinite(y2) && std::fabs(y1) < 1e6 && std::fabs(y2) < 1e6) {
                if ((y1 * y2 < 0) || (y1 == 0) || (y2 == 0))
                    intervals.emplace_back(x1, x2);
            }
            x1 = x2;
            x2 += step;
        }
        return intervals;
    }

    void addUnique(std::vector<double>& roots, double x, double eps) {
        for (double r : roots)
            if (std::fabs(x - r) < 10 * eps) return;
        roots.push_back(x);
    }

    double minAbsError(const std::vector<double>& reference, double x) {
        double min_err = std::numeric_limits<double>::max();
        for (double r : reference) {
            double err = std::fabs(x - r);
            if (err < min_err) min_err = err;
        }
        return min_err;
    }
}