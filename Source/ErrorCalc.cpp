#include "ErrorCalc.h"

//initialize object for single value calc
ErrorCalc::ErrorCalc(const QVector<double>& y_num, const QVector<double>& y_aly, const int &type)
	:	y_num(y_num),
		y_aly(y_aly),
		type(type)
{
}

//initialize object for simple graph calc
ErrorCalc::ErrorCalc(const QVector<double>& y_num, const QVector<double>& y_aly, const QVector<double>& x_num, const int &type)
	:	y_num(y_num),
		y_aly(y_aly),
		x_num(x_num),
		type(type)
{
}

ErrorCalc::~ErrorCalc()
{
	//nothing allocated on heap
}

//Launch process of calculation, data is taken from initialized object
void ErrorCalc::process()
{
	//main variable for storing final error value
	double err = 0;
	if(type == 0)
	{
		//Mean Squared Error Method
		for (int i = 0; i < y_num.size(); i++)
		{
			err += std::pow(y_num[i] - y_aly[i], 2.0);
		}
		err /= y_num.size();
		emit result(err);
	}else if(type == 1)
	{
		//Mean Absolute Error Method
		for (int i = 0; i < y_num.size(); i++)
		{
			err += std::abs(y_num[i] - y_aly[i]);
		}
		err /= y_num.size();
		emit result(err);
	}else if(type == 2)
	{
		//Maximal Error Method
		err = std::abs(y_num[0] - y_aly[0]);
		double curr_err;
		for (int i = 1; i < y_num.size(); i++)
		{
			curr_err = std::abs(y_num[i] - y_aly[i]);
			if (curr_err > err)
				err = curr_err;
		}
		emit result(err);
	}else if(type == 3)
	{
		//Plotting graph for Absoulte Error across numerical and analytical graphs
		QVector<double> error_vector;
		for (int i = 0; i < y_num.size(); i++)
		{
			error_vector.push_back(std::abs(y_num[i] - y_aly[i]));
			
		}
		emit result(PlotDataObject(x_num, error_vector));
	}
}

