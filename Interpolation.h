#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <iostream>
#include <vector>
#include <cmath>

// ================== INTERPOLACJA ==================
namespace Interpolation {

    double lagrangeInterpolation(
        const std::vector<double>& x,
        const std::vector<double>& y,
        double xi
    );

    double calculateMSE(
        const std::vector<double>& x_all,
        const std::vector<double>& y_all,
        const std::vector<double>& x_nodes,
        const std::vector<double>& y_nodes
    );

    int findOptimalNodes(
        const std::vector<double>& x,
        const std::vector<double>& y,
        int max_nodes = 0
    );

    std::vector<double> selectEveryKthPoint(
        const std::vector<double>& data,
        int k
    );
    std::vector<double> dividedDifferences(const std::vector<double>& x, const std::vector<double>& y);

    double evaluateNewton(const std::vector<double>& coeffs, const std::vector<double>& x, double value);

    double evaluateNatural(const std::vector<double>& coeffs, double x);

    double evaluateHorner(const std::vector<double>& coeffs, double x);

}
#endif
