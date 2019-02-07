#pragma once
#include "Solver.h"

class Euler_Improved : public Solver {
private:
	double get_next_y(const double &x, const double &y, double(*f)(const double &x, const double &y));
	Euler_Improved();
public:
	explicit Euler_Improved(const double &x0, const double &y0, const double &h, const double &x_to, double(*f)(const double &x, const double &y));
	explicit Euler_Improved(const double& x0, const double& y0, const double& x_to, const double& h_from, const double& h_to,
		const double& h_steps, double(*f)(const double& x, const double& y),
		double(*f_aly)(const double& x, const double& c), double(*ccalc)(const double& x, const double& y), const int &method);
};