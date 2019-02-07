#pragma once
#include "RungeKutta.h"

//redefined function for next point calculation
double RungeKutta::get_next_y(const double &x, const double &y, double(*f)(const double &x, const double &y)) {
	const double k1 = f(x, y);
	const double k2 = f(x + (this->h / 2), y + (this->h*k1 / 2));
	const double k3 = f(x + (this->h / 2), y + (this->h*k2 / 2));
	const double k4 = f(x + this->h, y + this->h*k3);
	return y + this->h*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

//initialize object for regular numeric graph
RungeKutta::RungeKutta(const double &x0, const double &y0, const double &h, const double &x_to, double(*f)(const double &x, const double &y))
	: Solver(x0, y0, h, x_to, f)
{
}

//initialize object for error to step-size graph
RungeKutta::RungeKutta(const double& x0, const double& y0, const double& x_to, const double& h_from, const double& h_to,
	const double& h_steps, double(*f)(const double& x, const double& y),
	double(*f_aly)(const double& x, const double& c), double(*ccalc)(const double& x, const double& y), const int &method)
	: Solver(x0, y0, x_to, h_from, h_to, h_steps, f, f_aly, ccalc, method)
{
}
