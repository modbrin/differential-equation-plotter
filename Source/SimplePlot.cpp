#include "SimplePlot.h"

SimplePlot::SimplePlot(const double& x0, const double& y0, const double& h, const double& x_to, double (* f)(const double& x, const double &c), double (*ccalc)(const double& x, const double &y))
	:	x0(x0),
		y0(y0),
		h(h),
		x_to(x_to)
{
	this->f = f;
	this->ccalc = ccalc;
}
SimplePlot::~SimplePlot()
{
	//nothing allocated on heap
}

PlotDataObject SimplePlot::solve()
{
	QVector<QVector<double>> result(2, QVector<double>());
	double current_x = this->x0;
	const double c = ccalc(this->x0, this->y0);
	while (current_x < this->x_to)
	{
		result[0].push_back(current_x);
		result[1].push_back(this->f(current_x, c));
		current_x += this->h;
	}
	result[0].push_back(this->x_to);
	result[1].push_back(this->f(this->x_to, c));
	return PlotDataObject(result[0], result[1]);
}

void SimplePlot::setH(const double& h)
{
	this->h = h;
}

void SimplePlot::process()
{
	qInfo("start simple plot in thread");
	emit result(solve());
}

