#include "DifferentialEquations.h"

namespace DifferentialEquations {
    double equation(double T) {
        const double k = 7e-12;
        return -k * std::pow(T, 4);
    }

    double euler(double T0, double dt, double t) {
        double T = T0;
        int steps = static_cast<int>(t / dt);
        for (int i = 0; i < steps; i++) {
            T += dt * equation(T);
            if (T < 0) {
                T = 0;
                break;
            }
        }
        return T;
    }

    double heun(double T0, double dt, double t) {
        double T = T0;
        int steps = static_cast<int>(t / dt);
        for (int i = 0; i < steps; i++) {
            double k1 = equation(T);
            double k2 = equation(T + dt * k1);
            T += (dt / 2.0) * (k1 + k2);
            if (T < 0) {
                T = 0;
                break;
            }
        }
        return T;
    }

    double midpoint(double T0, double dt, double t) {
        double T = T0;
        int steps = static_cast<int>(t / dt);
        for (int i = 0; i < steps; i++) {
            double k1 = equation(T);
            double k2 = equation(T + (dt / 2.0) * k1);
            T += dt * k2;
            if (T < 0) {
                T = 0;
                break;
            }
        }
        return T;
    }

    double runge_kutta(double T0, double dt, double t) {
        double T = T0;
        int steps = static_cast<int>(t / dt);
        for (int i = 0; i < steps; i++) {
            double k1 = equation(T);
            double k2 = equation(T + (dt / 2.0) * k1);
            double k3 = equation(T + (dt / 2.0) * k2);
            double k4 = equation(T + dt * k3);
            T += (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
            if (T < 0) {
                T = 0;
                break;
            }
        }
        return T;
    }

    double analytical_solution(double T) {
        return 970000 / std::cbrt(19166133 * T + 1000000000);
    }
}