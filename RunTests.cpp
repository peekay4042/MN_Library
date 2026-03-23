#include <iostream>

void run_tests_approximation();
void run_tests_integration();
void run_tests_interpolation();
void run_tests_linear_algebra();
void run_tests_nonlinear_algebra();
void run_tests_differential_equations();

int main() {
    run_tests_approximation();
    run_tests_integration();
    run_tests_interpolation();
    run_tests_linear_algebra();
    run_tests_nonlinear_algebra();
    run_tests_differential_equations();
    std::cout << "Wszystkie testy zakonczone.\n";
    return 0;
}