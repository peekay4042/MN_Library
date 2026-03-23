#include "NonlinearAlgebra.h"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace NonlinearEquations;

// Przyk³adowa funkcja: x^3 - x - 2 = 0
double f(double x) { return x * x * x - x - 2; }
double df(double x) { return 3 * x * x - 1; }

int main3() {
    std::cout << "Rozwiazywanie rownania x^3 - x - 2 = 0\n";
    std::cout << std::fixed << std::setprecision(15);

    // Metoda bisekcji
    double root_bisect = bisection(f, 1, 2, 1e-8, 100);
    std::cout << "Bisekcja: pierwiastek w przyblizeniu = " << root_bisect << "\n";

    // Metoda Newtona
    double root_newton = newton(f, df, 1.5, 1e-8, 100);
    std::cout << "Newton: pierwiastek w przyblizeniu = " << root_newton << "\n";

    // Metoda siecznych
    double root_secant = secant(f, 1, 2, 1e-8, 100);
    std::cout << "Sieczne: pierwiastek w przyblizeniu = " << root_secant << "\n";

    // Porównanie z wartoœci¹ referencyjn¹
    double ref = 1.5213797;
    std::cout << "Blad bisekcji: " << std::abs(root_bisect - ref) << "\n";
    std::cout << "Blad Newtona: " << std::abs(root_newton - ref) << "\n";
    std::cout << "Blad siecznych: " << std::abs(root_secant - ref) << "\n";

    return 0;
}