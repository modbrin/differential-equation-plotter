#pragma once
#include <QVector>

//container for holding plot data, i.e. x and y arrays
class PlotDataObject 
{
public:
	QVector<double> x, y;
	PlotDataObject(QVector<double> x, QVector<double> y);
	PlotDataObject();
	void clear();
};