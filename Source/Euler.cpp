#pragma once
#include "Euler.h"

//redefined function for next point calculation
double Euler::get_next_y(const double &x, const double &y, double(*f)(const double &x, const double &y)) {
	return y + this->h * f(x, y);
}

//initialize object for regular numeric graph
Euler::Euler(const double &x0, const double &y0, const double &h, const double &x_to, double(*f)(const double &x, const double &y))
	: Solver(x0, y0, h, x_to, f) 
{
}

//initialize object for error to step-size graph
Euler::Euler(const double& x0, const double& y0, const double& x_to, const double& h_from, const double& h_to,
	const double& h_steps, double(*f)(const double& x, const double& y),
	double(*f_aly)(const double& x, const double& c), double(*ccalc)(const double& x, const double& y), const int &method)
	: Solver(x0, y0, x_to, h_from, h_to, h_steps, f, f_aly, ccalc, method)
{
}
