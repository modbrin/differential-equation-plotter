#include "PlotDataObject.h"

//object initialization
PlotDataObject::PlotDataObject(QVector<double> x, QVector<double> y)
	:	x(x),
		y(y)
{
}

PlotDataObject::PlotDataObject()
{
}

//release memory used by container's data
void PlotDataObject::clear()
{
	this->x.clear();
	this->y.clear();
}
