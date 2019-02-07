#pragma once
#include <QObject>
#include "PlotDataObject.h"

//first version of Solver was a class template, but since it's incompatible with QObject it had to be specialized into double
class Solver : public QObject{
	Q_OBJECT
protected:
	double x0, y0, h, x_to, h_from, h_to, h_steps;
	double(*f)(const double &x, const double &y);
	double(*f_aly)(const double &x, const double &c);
	double(*ccalc)(const double &x, const double &y);
	virtual double get_next_y(const double &x, const double &y, double(*f)(const double &x, const double &y));
	Solver();
	int err_method;
	double getMSE(const PlotDataObject &numerical) const;
	double getMAE(const PlotDataObject &numerical) const;
	double getMXE(const PlotDataObject &numerical) const;
public:
	Solver(const double &x0, const double &y0, const double &h, const double &x_to, double(*f)(const double &x, const double &y));
	Solver(const double &x0, const double &y0, const double &x_to, const double &h_from, const double &h_to, const double &h_steps, double(*f)(const double &x, const double &y), double(*f_aly)(const double &x, const double &c), double(*ccalc)(const double &x, const double &y), const int &method);
	~Solver();
	PlotDataObject solve();
public slots:
	void process();
	void err_proc();
signals:
	void finished();
	void error(QString err);
	void result(PlotDataObject obj);
};