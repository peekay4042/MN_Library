#include "Interpolation.h"
#include <iostream>
#include <vector>

// ================== INTERPOLACJA ==================
namespace Interpolation {

    double lagrangeInterpolation(
        const std::vector<double>& x,
        const std::vector<double>& y,
        double xi) {

        if (x.size() != y.size())
            throw std::invalid_argument("Rozmiary x i y musz¹ byæ równe");

        double suma = 0.0;
        int n = x.size();

        for (int i = 0; i < n; i++) {
            double l = 1.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    l *= (xi - x[j]) / (x[i] - x[j]);
                }
            }
            suma += l * y[i];
        }
        return suma;
    }

    double calculateMSE(
        const std::vector<double>& x_all,
        const std::vector<double>& y_all,
        const std::vector<double>& x_nodes,
        const std::vector<double>& y_nodes) {

        double suma_bl = 0.0;
        int licz = 0;

        for (int i = 0; i < x_all.size(); i++) {
            if (std::find(x_nodes.begin(), x_nodes.end(), x_all[i]) == x_nodes.end()) {
                double interpol = lagrangeInterpolation(x_nodes, y_nodes, x_all[i]);
                double bl = std::pow(y_all[i] - interpol, 2);
                suma_bl += bl;
                licz++;
            }
        }

        return licz > 0 ? suma_bl / licz : 0.0;
    }

    int findOptimalNodes(
        const std::vector<double>& x,
        const std::vector<double>& y,
        int max_nodes) {

        if (max_nodes == 0) max_nodes = x.size() / 2;

        int best_k = 0;
        double min_mse = std::numeric_limits<double>::max();

        for (int k = 2; k <= max_nodes; k++) {
            std::vector<double> x_k, y_k;
            int step = x.size() / k;
            if (step == 0) step = 1;

            for (int i = 0; i < x.size(); i += step) {
                x_k.push_back(x[i]);
                y_k.push_back(y[i]);
            }

            double mse = calculateMSE(x, y, x_k, y_k);
            if (mse < min_mse) {
                min_mse = mse;
                best_k = k;
            }
        }

        return best_k;
    }

    std::vector<double> selectEveryKthPoint(
        const std::vector<double>& data,
        int k) {

        std::vector<double> result;
        if (data.empty() || k <= 0) return result;

        result.push_back(data[0]);
        if (k >= data.size()) return result;

        for (int i = 1; i < data.size(); i++) {
            if (i % k == 0) {
                result.push_back(data[i]);
            }
        }
        return result;
    }

    std::vector<double> dividedDifferences(
        const std::vector<double>& x_values,
        const std::vector<double>& y_values
    ) {
        int n = x_values.size();
        if (n == 0 || y_values.size() == 0)
            return {};

        std::vector<std::vector<double>> table(n, std::vector<double>(n));
        std::vector<double> coefficients(n);

        for (int i = 0; i < n; ++i)
            table[i][0] = y_values[i];

        for (int j = 1; j < n; ++j)
            for (int i = 0; i < n - j; ++i)
                table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x_values[i + j] - x_values[i]);

        for (int i = 0; i < n; ++i)
            coefficients[i] = table[0][i];

        return coefficients;
    }

    // Interpolation::evaluateNewton
    double evaluateNewton(
        const std::vector<double>& coefficients,
        const std::vector<double>& x_values,
        double x
    ) {
        double result = coefficients[0];
        double product = 1.0;

        for (int i = 1; i < coefficients.size(); ++i) {
            product *= (x - x_values[i - 1]);
            result += coefficients[i] * product;
        }

        return result;
    }

    // Interpolation::evaluateNatural
    double evaluateNatural(const std::vector<double>& coefficients, double x) {
        double result = 0.0;
        double power = 1.0;

        for (double c : coefficients) {
            result += c * power;
            power *= x;
        }

        return result;
    }

    // Interpolation::evaluateHorner
    double evaluateHorner(const std::vector<double>& coefficients, double x) {
        double result = 0.0;

        for (int i = coefficients.size() - 1; i >= 0; --i)
            result = result * x + coefficients[i];

        return result;
    }
}