#pragma once
#include <QObject>
#include "PlotDataObject.h"

//class for creating plot dataset from analytical solution
class SimplePlot : public QObject
{
	Q_OBJECT
private:
	double x0, y0, h, x_to, c;
	double(*f)(const double &x, const double &c);
	double(*ccalc)(const double &x, const double &y);
public:
	SimplePlot(const double &x0, const double &y0, const double &h, const double &x_to, double(*f)(const double &x, const double &c), double(*ccalc)(const double &x, const double &y));
	~SimplePlot();
	PlotDataObject solve();
	void setH(const double &h);
public slots :
	void process();
signals:
	void finished();
	void error(QString err);
	void result(PlotDataObject obj);
};
