#include "DifferentialEquations.h"
#include <iostream>
#include <iomanip>

using namespace DifferentialEquations;

int main1() {
    double T0 = 1200.0; // temperatura pocz¹tkowa
    double dt = 0.5;   // krok czasowy
    double t_max = 1000.0; // czas koñcowy

    std::cout << "Symulacja chlodzenia ciala (Euler vs Runge-Kutta)\n";
    std::cout << std::fixed << std::setprecision(15);
    std::cout << "Czas\tEuler\tRunge-Kutta\n";
    for (double t = 0; t <= t_max; t += 100.0) {
        double T_euler = euler(T0, dt, t);
        double T_rk = runge_kutta(T0, dt, t);
        std::cout << t << "\t" << T_euler << "\t" << T_rk << "\n";
    }
    return 0;
}