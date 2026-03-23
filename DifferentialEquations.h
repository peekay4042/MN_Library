#ifndef DIFFERENTIAL_EQUATIONS_H
#define DIFFERENTIAL_EQUATIONS_H

#include <vector>
#include <cmath>

namespace DifferentialEquations {
	double equation(double T);

	double euler(
		double T0,
		double dt,
		double t);

	double heun(
		double T0,
		double dt,
		double t);

	double midpoint(
		double T0,
		double dt,
		double t);

	double runge_kutta(
		double T0,
		double dt,
		double t);

	double analytical_solution(double T);
}

#endif // !DIFFERENTIAL_EQUATIONS_H

