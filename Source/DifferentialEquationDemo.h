#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_DifferentialEquationDemo.h"
#include "PlotDataObject.h"

class DifferentialEquationDemo : public QMainWindow
{
	Q_OBJECT

public:
	DifferentialEquationDemo(QWidget *parent = Q_NULLPTR);

private:
	Ui::DifferentialEquationDemoClass ui;
	void initInputFields();
	void enableInputGeneral(const bool &isEnabled);
	bool emptyCheckMain();
	bool emptyCheckMain_E();
	bool emptyCheckErr();
	bool dataCheckMain(const double &x_init, const double &y_init, const double &h, const double &x_to);
	bool dataCheckMain_E(const double &x_init, const double &y_init, const double &x_to);
	bool dataCheckErr(const double &e_from, const double &e_to, const double &e_num);
	inline int getStepSize(const double &from, const double &to, const int &num);
	QVector<double> x_num, y_num, x_aly, y_aly;
	inline int DTOI(const double &val);
public slots:
	void launchMain();
	void catchMain(PlotDataObject);
	void launchSimple();
	void catchSimple(PlotDataObject);
	void setDefaults();
	void legendSwitch();
	void launchError();
	void errorWarning(QString warning);
	void launchSE(const int &type);
	void launchSE3();
	void catchSE(double err);
	void catchSE(PlotDataObject plt);
	void xplus();
	void yplus();
	void xyplus();
	void liCalc();
	void launchErrGraph();
	void catchErrGraph(PlotDataObject plt);
	void mSizeHandler();
	void mNumbHandler();
};
