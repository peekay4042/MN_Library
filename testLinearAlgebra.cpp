#include "LinearAlgebra.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>

using namespace LinearAlgebra;

void test_createIdentityMatrix() {
    // Test 1: Macierz jednostkowa 3x3
    auto I3 = createIdentityMatrix(3);
    assert(I3.size() == 3 && I3[0][0] == 1 && I3[1][1] == 1 && I3[2][2] == 1);

    // Test 2: Macierz jednostkowa 1x1
    auto I1 = createIdentityMatrix(1);
    assert(I1.size() == 1 && I1[0][0] == 1);
}

void test_multiplyMatrices() {
    // Test 1: Mno¿enie macierzy jednostkowej przez macierz
    std::vector<std::vector<double>> A = { {1,2},{3,4} };
    auto I = createIdentityMatrix(2);
    auto res = multiplyMatrices(A, I);
    assert(res == A);

    // Test 2: Mno¿enie dwóch macierzy 2x2
    std::vector<std::vector<double>> B = { {2,0},{1,2} };
    auto res2 = multiplyMatrices(A, B);
    assert(res2[0][0] == 4 && res2[0][1] == 4 && res2[1][0] == 10 && res2[1][1] == 8);
}

void test_areMatricesEqual() {
    // Test 1: Identyczne macierze
    std::vector<std::vector<double>> A = { {1,2},{3,4} };
    std::vector<std::vector<double>> B = { {1,2},{3,4} };
    assert(areMatricesEqual(A, B, 1e-12));

    // Test 2: Ró¿ne macierze
    std::vector<std::vector<double>> C = { {1,2},{3,5} };
    assert(!areMatricesEqual(A, C, 1e-12));
}

void test_luDecomposition() {
    // Test 1: Prosta macierz 2x2
    std::vector<std::vector<double>> A = { {4,3},{6,3} };
    std::vector<std::vector<double>> L, U, P;
    luDecomposition(A, L, U, P);
    auto PA = multiplyMatrices(P, A);
    auto LU = multiplyMatrices(L, U);
    assert(areMatricesEqual(PA, LU, 1e-9));

    // Test 2: Macierz jednostkowa
    auto I = createIdentityMatrix(2);
    luDecomposition(I, L, U, P);
    assert(areMatricesEqual(L, I, 1e-12) && areMatricesEqual(U, I, 1e-12));
}

void test_solveLowerTriangular() {
    // Test 1: L = jednostkowa, b = [1,2]
    auto L = createIdentityMatrix(2);
    std::vector<double> b = { 1,2 };
    auto y = solveLowerTriangular(L, b, createIdentityMatrix(2));
    assert(y == b);

    // Test 2: L = [[1,0],[2,1]], b = [3,7]
    L = { {1,0},{2,1} };
    b = { 3,7 };
    auto y2 = solveLowerTriangular(L, b, createIdentityMatrix(2));
    assert(std::abs(y2[0] - 3) < 1e-12 && std::abs(y2[1] - 1) < 1e-12);
}

void test_solveUpperTriangular() {
    // Test 1: U = jednostkowa, y = [1,2]
    auto U = createIdentityMatrix(2);
    std::vector<double> y = { 1,2 };
    auto x = solveUpperTriangular(U, y);
    assert(x == y);

    // Test 2: U = [[2,1],[0,1]], y = [5,2]
    U = { {2,1},{0,1} };
    y = { 5,2 };
    auto x2 = solveUpperTriangular(U, y);
    assert(std::abs(x2[1] - 2) < 1e-12 && std::abs(x2[0] - 1.5) < 1e-12);
}

void test_solveLU() {
    // Test 1: Uk³ad 2x2
    std::vector<std::vector<double>> A = { {2,1},{1,3} };
    std::vector<double> b = { 5,10 };
    auto x = solveLU(A, b);
    assert(std::abs(x[0] - 1) < 1e-12 && std::abs(x[1] - 3) < 1e-12);

    // Test 2: Uk³ad 1x1
    A = { {4} };
    b = { 8 };
    x = solveLU(A, b);
    assert(std::abs(x[0] - 2) < 1e-12);
}

void test_gaussianElimination() {
    // Test 1: Uk³ad 2x2
    std::vector<std::vector<double>> A = { {2,1},{1,3} };
    std::vector<double> b = { 5,10 };
    gaussianElimination(2, b, A, false);
    assert(std::abs(A[1][0]) < 1e-9);

    // Test 2: Uk³ad sprzeczny
    A = { {1,1},{1,1} };
    b = { 2,3 };
    gaussianElimination(2, b, A, false);
    assert(std::abs(A[1][0]) < 1e-9);
}

void test_backSubstitution() {
    // Test 1: Uk³ad 2x2 górnotrójk¹tny
    std::vector<std::vector<double>> A = { {2,1},{0,1} };
    std::vector<double> b = { 5,2 };
    auto x = backSubstitution(2, A, b);
    assert(std::abs(x[0] - 1.5) < 1e-12 && std::abs(x[1] - 2) < 1e-12);

    // Test 2: Uk³ad 1x1
    A = { {4} };
    b = { 8 };
    x = backSubstitution(1, A, b);
    assert(std::abs(x[0] - 2) < 1e-12);
}

void test_verifySystem() {
    // Test 1: Poprawne rozwi¹zanie
    std::vector<std::vector<double>> A = { {2,1},{1,3} };
    std::vector<double> x = { 1,3 };
    std::vector<double> b = { 5,10 };
    assert(verifySystem(2, A, x, b, 1e-12));

    // Test 2: B³êdne rozwi¹zanie
    x = { 0,0 };
    assert(!verifySystem(2, A, x, b, 1e-12));
}

void test_createIdentityMatrix_errors() {
    // Rozmiar 0
    auto I = createIdentityMatrix(0);
    assert(I.empty());
}

void test_multiplyMatrices_errors() {
    // Niezgodne rozmiary
    std::vector<std::vector<double>> A = { {1,2} };
    std::vector<std::vector<double>> B = { {1,2},{3,4},{5,6} };
    try {
        auto res = multiplyMatrices(A, B);
        assert(false && "Brak wyj¹tku dla niezgodnych rozmiarów");
    }
    catch (...) {
        // Oczekiwany wyj¹tek lub b³¹d
    }
}

void test_solveLU_errors() {
    // Macierz osobliwa
    std::vector<std::vector<double>> A = { {0,0},{0,0} };
    std::vector<double> b = { 0,0 };
    auto x = solveLU(A, b);
    for (auto v : x) assert(std::isnan(v) || v == 0.0);
}

void run_tests_linear_algebra() {
    test_createIdentityMatrix();
    test_multiplyMatrices();
    test_areMatricesEqual();
    test_luDecomposition();
    test_solveLowerTriangular();
    test_solveUpperTriangular();
    test_solveLU();
    test_gaussianElimination();
    test_backSubstitution();
    test_verifySystem();
    test_createIdentityMatrix_errors();
    test_multiplyMatrices_errors();
    test_solveLU_errors();
    std::cout << "Wszystkie testy LinearAlgebra.h zaliczone!" << std::endl;
}