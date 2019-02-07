#pragma once
#include <QObject>
#include <QVector>
#include "PlotDataObject.h"

class ErrorCalc : public QObject
{
	Q_OBJECT
private:
	QVector<double> y_num, y_aly, x_num;
	int type;
public:
	ErrorCalc(const QVector<double> &y_num, const QVector<double> &y_aly, const int &type);
	ErrorCalc(const QVector<double> &y_num, const QVector<double> &y_aly, const QVector<double> &x_num, const int &type);
	~ErrorCalc();
public slots:
	void process();
signals:
	void finished();
	void error(QString err);
	void result(double err);
	void result(PlotDataObject plt);
};
