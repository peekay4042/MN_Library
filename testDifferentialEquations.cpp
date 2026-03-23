#include "DifferentialEquations.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>

using namespace DifferentialEquations;

// Testy dla equation
void test_equation() {
    // Poprawny przypadek
    double T = 100.0;
    double result = equation(T);
    assert(result < 0);

    // Przypadek brzegowy: T = 0
    assert(equation(0.0) == 0.0);
}

// Testy dla euler
void test_euler() {
    // Poprawny przypadek
    double T = euler(100.0, 0.1, 1.0);
    assert(T >= 0 && T < 100.0);

    // B³êdny przypadek: dt = 0 (powinno zwróciæ T0)
    double T0 = 50.0;
    assert(euler(T0, 0.0, 1.0) == T0);
}

// Testy dla heun
void test_heun() {
    // Poprawny przypadek
    double T = heun(100.0, 0.1, 1.0);
    assert(T >= 0 && T < 100.0);

    // B³êdny przypadek: t = 0 (powinno zwróciæ T0)
    double T0 = 80.0;
    assert(heun(T0, 0.1, 0.0) == T0);
}

// Testy dla midpoint
void test_midpoint() {
    // Poprawny przypadek
    double T = midpoint(100.0, 0.1, 1.0);
    assert(T >= 0 && T < 100.0);

    // B³êdny przypadek: dt < 0 (nie powinno zmieniaæ T0)
    double T0 = 60.0;
    assert(midpoint(T0, -0.1, 1.0) == T0);
}

// Testy dla runge_kutta
void test_runge_kutta() {
    // Poprawny przypadek
    double T = runge_kutta(100.0, 0.1, 1.0);
    assert(T >= 0 && T < 100.0);

    // B³êdny przypadek: t < 0 (nie powinno zmieniaæ T0)
    double T0 = 70.0;
    assert(runge_kutta(T0, 0.1, -1.0) == T0);
}

// Testy dla analytical_solution
void test_analytical_solution() {
    // Poprawny przypadek
    double T = 100.0;
    double result = analytical_solution(T);
    assert(result > 0);

    // B³êdny przypadek: bardzo du¿e T
    double bigT = 1e12;
    double res = analytical_solution(bigT);
    assert(res >= 0);
}

void run_tests_differential_equations() {
    test_equation();
    test_euler();
    test_heun();
    test_midpoint();
    test_runge_kutta();
    test_analytical_solution();
    std::cout << "Wszystkie testy DifferentialEquations.h zaliczone!" << std::endl;
}