#include "NonlinearAlgebra.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>

using namespace NonlinearEquations;

double f1(double x) { return x * x - 4; }
double df1(double x) { return 2 * x; }
double f2(double x) { return std::sin(x); }
double df2(double x) { return std::cos(x); }

void test_bisection() {
    // Test 1: Pierwiastek równania x^2-4=0 w [0,3]
    double root = bisection(f1, 0, 3, 1e-8, 100, false, nullptr, "");
    assert(std::abs(root - 2.0) < 1e-7);

    // Test 2: Brak zmiany znaku w [3,5] (powinno zwróciæ NAN)
    double nan_root = bisection(f1, 3, 5, 1e-8, 100, false, nullptr, "");
    assert(std::isnan(nan_root));
}

void test_newton() {
    // Test 1: Pierwiastek x^2-4=0, start x0=3
    double root = newton(f1, df1, 3, 1e-8, 100, false, nullptr, "");
    assert(std::abs(root - 2.0) < 1e-7);

    // Test 2: Brak pochodnej (df=0), x0=0 (powinno zwróciæ NAN)
    double nan_root = newton(f1, df1, 0, 1e-8, 10, false, nullptr, "");
    assert(std::isnan(nan_root));
}

void test_secant() {
    // Test 1: Pierwiastek x^2-4=0, start x0=0, x1=3
    double root = secant(f1, 0, 3, 1e-8, 100, false, nullptr, "");
    assert(std::abs(root - 2.0) < 1e-7);

    // Test 2: Brak zmiany wartoœci (f0==f1), x0=1, x1=1 (powinno zwróciæ NAN)
    double nan_root = secant(f1, 1, 1, 1e-8, 10, false, nullptr, "");
    assert(std::isnan(nan_root));
}

void test_findIntervals() {
    // Test 1: sin(x) w [0, 4], krok 1
    auto intervals = findIntervals(f2, 0, 4, 1.0);
    assert(!intervals.empty());

    // Test 2: x^2-4 w [0, 3], krok 0.5 (powinien znaleŸæ przedzia³ z pierwiastkiem)
    auto intervals2 = findIntervals(f1, 0, 3, 0.5);
    bool found = false;
    for (auto& p : intervals2)
        if (p.first <= 2.0 && p.second >= 2.0) found = true;
    assert(found);
}

void test_addUnique() {
    // Test 1: Dodanie unikalnej wartoœci
    std::vector<double> roots = { 1.0, 2.0 };
    addUnique(roots, 3.0, 1e-8);
    assert(roots.size() == 3);

    // Test 2: Próba dodania istniej¹cej wartoœci (nie powinno dodaæ)
    addUnique(roots, 2.0 + 1e-9, 1e-8);
    assert(roots.size() == 3);
}

void test_minAbsError() {
    // Test 1: Najmniejszy b³¹d bezwzglêdny
    std::vector<double> ref = { 1.0, 2.0, 3.0 };
    double err = minAbsError(ref, 2.1);
    assert(std::abs(err - 0.1) < 1e-12);

    // Test 2: x równe jednemu z referencyjnych
    err = minAbsError(ref, 3.0);
    assert(std::abs(err) < 1e-12);
}

void test_bisection_errors() {
    // Przedzia³ o zerowej d³ugoœci
    double root = bisection(f1, 2.0, 2.0, 1e-8, 100, false, nullptr, "");
    assert(std::isnan(root));
}

void test_newton_errors() {
    // Maksymalna liczba iteracji = 0
    double root = newton(f1, df1, 3, 1e-8, 0, false, nullptr, "");
    assert(std::isnan(root) || std::abs(root - 3) < 1e-12);
}

void test_secant_errors() {
    // Maksymalna liczba iteracji = 0
    double root = secant(f1, 0, 3, 1e-8, 0, false, nullptr, "");
    assert(std::isnan(root) || std::abs(root - 3) < 1e-12);
}

void run_tests_nonlinear_algebra() {
    test_bisection();
    test_newton();
    test_secant();
    test_findIntervals();
    test_addUnique();
    test_minAbsError();
    test_bisection_errors();
    test_newton_errors();
    test_secant_errors();
    std::cout << "Wszystkie testy NonlinearAlgebra.h zaliczone!" << std::endl;
}