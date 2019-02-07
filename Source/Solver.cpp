#pragma once
#include "Solver.h"
#include "SimplePlot.h"

//common function for getting next point, redefined in child classes
double Solver::get_next_y(const double& x, const double& y, double (*f)(const double& x, const double& y))
{
	return double();
}

//initialize object for regular numeric graph
Solver::Solver(const double &x0, const double &y0, const double &h, const double &x_to, double(*f)(const double &x, const double &y))
	:	x0(x0),
		y0(y0),
		h(h),
		x_to(x_to)
{
	this->f = f;
}

//initialize object for error to step-size graph
Solver::Solver(const double& x0, const double& y0, const double& x_to, const double& h_from, const double& h_to,
	const double& h_steps, double(* f)(const double& x, const double& y),
	double(* f_aly)(const double& x, const double& c), double(* ccalc)(const double& x, const double& y), const int &method)
	:	x0(x0),
		y0(y0),
		x_to(x_to),
		h(h_from),
		h_from(h_from),
		h_to(h_to),
		h_steps(h_steps),
		err_method(method)
{
	this->f = f;
	this->f_aly = f_aly;
	this->ccalc = ccalc;
}

Solver::~Solver()
{
	//nothing allocated on heap
}

//Get numerical solution of constructed numerical method object
//x_to must be larger than x0
//step h must be smaller than (x_to - x0)
PlotDataObject Solver::solve() {
	const int NUM_OF_POINTS = ceil((this->x_to - this->x0) / this->h);
	QVector<QVector<double>> result(2, QVector<double>(NUM_OF_POINTS + 1));
	result[0][0] = this->x0;
	result[1][0] = this->y0;
	for (int i = 1; i <= NUM_OF_POINTS; i++) {
		result[0][i] = result[0][i - 1] + this->h;
		result[1][i] = get_next_y(result[0][i - 1], result[1][i - 1], f);
	}
	return PlotDataObject(result[0], result[1]);
}

void Solver::process()
{
	qInfo("stating solving");
	emit result(solve());
}

//calculate Mean Squared Error from numerical graph with respect to local analytical function
double Solver::getMSE(const PlotDataObject &numerical) const
{
	const double c = ccalc(x0, y0);
	double err = 0;
	//Mean Squared Error Method
	for (int i = 0; i < numerical.y.size(); i++)
	{
		err += std::pow(numerical.y[i] - f_aly(numerical.x[i], c), 2.0);
	}
	err /= numerical.y.size();
	return err;
}

//calculate Mean Absolute Error from numerical graph with respect to local analytical function
double Solver::getMAE(const PlotDataObject &numerical) const
{
	const double c = ccalc(x0, y0);
	double err = 0;
	//Mean Absolute Error Method
	for (int i = 0; i < numerical.y.size(); i++)
	{
		err += std::abs(numerical.y[i] - f_aly(numerical.x[i], c));
	}
	err /= numerical.y.size();
	return err;
}

//calculate Maximum Error from numerical graph with respect to local analytical function
double Solver::getMXE(const PlotDataObject &numerical) const
{
	const double c = ccalc(x0, y0);
	double max_err = 0;
	double curr_err;
	for (int j = 0; j < numerical.y.size(); ++j)
	{
		//find largest error
		curr_err = std::abs(numerical.y[j] - f_aly(numerical.x[j], c));
		if (curr_err > max_err)
			max_err = curr_err;
	}
	return max_err;
}

//launch process with error to step-size graph building
void Solver::err_proc()
{	
	//initialize lists for x and y
	QVector<double> err_list_x, err_list_y;
	//get outer step size for sampling
	const double step_size = (h_to - h_from) / h_steps;
	
	//for every outer step
	for(int i = 0; i < h_steps; ++i)
	{	
		//get inner step size, i.e. step size for actual numeric method
		double c_step_size = (x_to - x0) / static_cast<int>(ceil(h_from));
		//assign it to variable
		this->h = c_step_size;
		
		//make record of that step size as x axis
		err_list_x.push_back(c_step_size);

		//solve for current numeric method
		PlotDataObject sol_num = solve();

		//get error from resulting data with respect to chosen error-calculation method
		switch(err_method)
		{
			case 0:	err_list_y.push_back(getMSE(sol_num)); break;
			case 1: err_list_y.push_back(getMAE(sol_num)); break;
			case 2: err_list_y.push_back(getMXE(sol_num)); break;
		}
		
		//increment outer step size
		h_from+=step_size;

		//set new outer step size as current
		this->h = h_from;
	}
	//return resulting data
	emit result(PlotDataObject(err_list_x, err_list_y));
}
