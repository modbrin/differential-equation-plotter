#include "DifferentialEquationDemo.h"
#include "PlotDataObject.h"
#include "Euler.h"
#include "EulerImp.h"
#include "RungeKutta.h"
#include "Solver.h"
#include "SimplePlot.h"
#include "ErrorCalc.h"
#include <qmath.h>

#ifndef THREAD_WRAPPER(target)
	#define THREAD_WRAPPER(target) connect(target, SIGNAL(error(QString)), this, SLOT(errorWarning(QString)));\
	connect(thread, SIGNAL(started()), target, SLOT(process()));\
	connect(target, SIGNAL(finished()), thread, SLOT(quit()));\
	connect(target, SIGNAL(finished()), target, SLOT(deleteLater()));\
	connect(thread, SIGNAL(finished()), thread, SLOT(deleteLater()));
#endif

//margin for graph inside the window
const double BORDER_MARGIN = 0.01;

//list that contains default values for input fields
const QList<double> DEF_LS = { 0.0, 1.0, 9.0, 100, 10, 100, 90};

//That is differential equation of form dy/dx = f(x, y)
inline double e_19(const double &x, const double &y)
{
	return 2.0 * sqrt(y) + 2.0 * y;
}

//That is analytical solution of form y = f(x, C)
double e_19_sol(const double &x, const double &c) {
	return std::pow(std::pow(M_E, (c + x)) - 1, 2.0);
}

/*
 *That is expression of const in analytical solution of form C = f(x0, y0)
 *it actually doesn't mean that this c is the only one, that is taken 
 *as a single real root from family of equations with imaginary part as a component
 */
double e_19_const(const double &x, const double &y)
{
	return std::log(std::sqrt(y) + 1) - x;
}

//Main initialization function, binding button events and setting default values, plotting graph with default settings
DifferentialEquationDemo::DifferentialEquationDemo(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	//set default functionaliy for the UI
	initInputFields();

	//create signal mapping from UI buttons to function slots
	connect(ui.plot_main, SIGNAL(released()), this, SLOT(launchMain()));
	connect(ui.s_plot, SIGNAL(released()), this, SLOT(launchError()));
	connect(ui.reset_main, SIGNAL(released()), this, SLOT(setDefaults()));
	connect(ui.x_nav, SIGNAL(released()), this, SLOT(xplus()));
	connect(ui.y_nav, SIGNAL(released()), this, SLOT(yplus()));
	connect(ui.both_nav, SIGNAL(released()), this, SLOT(xyplus()));
	connect(ui.li_show, SIGNAL(released()), this, SLOT(liCalc()));
	connect(ui.s_cr_gr, SIGNAL(released()), this, SLOT(launchSE3()));
	connect(ui.m_number, SIGNAL(released()), this, SLOT(mNumbHandler()));
	connect(ui.m_size, SIGNAL(released()), this, SLOT(mSizeHandler()));
	connect(ui.e_plt, SIGNAL(released()), this, SLOT(launchErrGraph()));
	connect(ui.legend, SIGNAL(stateChanged(int)), this, SLOT(legendSwitch()));
	//register custom container so it's recognisable by singal-slot mechanism
	qRegisterMetaType<PlotDataObject>("PlotDataObject");
	
	//create header in window
	ui.customPlot->plotLayout()->insertRow(0);
	ui.customPlot->plotLayout()->addElement(0, 0, new QCPTextElement(ui.customPlot, "Equation 19: dy/dx = 2*y^0.5 + 2*y", QFont("Roboto Light", 20)));
	
	//set default values for input fields
	setDefaults();
	
	//make starter plot
	launchMain();
}

//set restriction on zoom for x axis only
void DifferentialEquationDemo::xplus()
{
	ui.customPlot->axisRect()->setRangeZoom(Qt::Horizontal);
	ui.y_nav->setChecked(false);
	ui.both_nav->setChecked(false);
}

//set restriction on zoom for y axis only
void DifferentialEquationDemo::yplus()
{
	ui.customPlot->axisRect()->setRangeZoom(Qt::Vertical);
	ui.x_nav->setChecked(false);
	ui.both_nav->setChecked(false);
}

//set no restrictios on zoom
void DifferentialEquationDemo::xyplus()
{
	ui.customPlot->axisRect()->setRangeZoom(Qt::Vertical | Qt::Horizontal);
	ui.x_nav->setChecked(false);
	ui.y_nav->setChecked(false);
}

//General error message function
void DifferentialEquationDemo::errorWarning(QString message)
{
	QMessageBox::warning(this, "Runtime Error", message);
}

//Set decimal regular expression for all input fields
void DifferentialEquationDemo::initInputFields()
{
	QFont pfont("Roboto");
	pfont.setPointSize(8);
	ui.customPlot->xAxis->setTickLabelFont(pfont);
	ui.customPlot->yAxis->setTickLabelFont(pfont);
	//ui.customPlot->xAxis->setAutoTickStep(false);   // <-- disable to use your own value
	//ui.customPlot->xAxis->setTickStep(7200);
	
	const QRegExp DECIMAL_ONLY("^-?[0-9]+([.][0-9]+)?$");
	const QRegExp POSITIVE_INT("^[1-9]\d*$");
	QRegExpValidator *validator_d = new QRegExpValidator(DECIMAL_ONLY, this);
	QRegExpValidator *validator_i = new QRegExpValidator(DECIMAL_ONLY, this);
	for (auto o : ui.Settings->findChildren<QLineEdit *>())
		o->setValidator(validator_d);
	ui.e_from->setValidator(validator_i);
	ui.e_to->setValidator(validator_i);
	ui.e_num->setValidator(validator_i);
}

//function for making input (not) availavle for entire application
void DifferentialEquationDemo::enableInputGeneral(const bool &isEnabled)
{
	for (auto o : ui.Settings->findChildren<QWidget *>())
		o->setEnabled(isEnabled);
	ui.customPlot->setInteraction(QCP::iRangeDrag, isEnabled);
	ui.customPlot->setInteraction(QCP::iRangeZoom, isEnabled);
}

//set default values for main input fields
void DifferentialEquationDemo::setDefaults()
{
	const QList<QLineEdit *> MN_LIST = { ui.x_init, ui.y_init, ui.x_final, ui.step_size, ui.e_from, ui.e_to, ui.e_num};
	for (int i = 0; i < DEF_LS.size(); i++)
		MN_LIST[i]->setText(QString::number(DEF_LS.at(i)));
	ui.customPlot->legend->setFont(QFont("Roboto", 9));
	ui.customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft | Qt::AlignTop);
}

//check main input fields to be non-empty
bool DifferentialEquationDemo::emptyCheckMain()
{
	if (ui.x_init->text() == "" || ui.y_init->text() == "" || ui.x_final->text() == "" || ui.step_size->text() == "")
	{
		QMessageBox::warning(this, "Invalid input", "Empty fields present");
		enableInputGeneral(true);
		return false;
	}
	return true;
}

//check main input fields to be non-empty, but in this case ignore step field since error graph has own step settings
bool DifferentialEquationDemo::emptyCheckMain_E()
{
	if (ui.x_init->text() == "" || ui.y_init->text() == "" || ui.x_final->text() == "")
	{
		QMessageBox::warning(this, "Invalid input", "Empty fields present");
		enableInputGeneral(true);
		return false;
	}
	return true;
}

//check main input fields's values to be valid in order to use them in algorithms
bool DifferentialEquationDemo::dataCheckMain(const double &x_init, const double &y_init, const double &h, const double &x_to)
{
	//restriction on initial values: right bound must be bigger than left bound, step size should be smaller than chosen interval
	if (x_init < x_to && (x_to - x_init) > h)
	{
		//for concrete equation number 19 y value should be bigger than 0 becuase of square root
		if(y_init>0.0)
		{
			return true;
		}else
		{
			QMessageBox::warning(this, "Invalid input", "Y must be grater than 0");
			enableInputGeneral(true);
			return false;
		}
	}
	else
	{
		QMessageBox::warning(this, "Invalid input", "Can not plot graph with such inputs");
		enableInputGeneral(true);
		return false;
	}
}

//check main input fields's values to be valid in order to use them in algorithms
bool DifferentialEquationDemo::dataCheckMain_E(const double &x_init, const double &y_init, const double &x_to)
{
	//restriction on initial values: right bound must be bigger than left bound, step size should be smaller than chosen interval
	if (x_init < x_to)
	{
		//for concrete equation number 19 y value should be bigger than 0 becuase of square root
		if (y_init>0.0)
		{
			return true;
		}
		else
		{
			QMessageBox::warning(this, "Invalid input", "Y must be grater than 0");
			enableInputGeneral(true);
			return false;
		}
	}
	else
	{
		QMessageBox::warning(this, "Invalid input", "Can not plot graph with such inputs");
		enableInputGeneral(true);
		return false;
	}
}

//simple search for min value in qvector
double minQV(const QVector<double> &vec)
{
	double curr = vec[0];
	for (double o : vec)
	{
		if (o < curr)
			curr = o;
	}
	return curr;
}

//simple search for max value in qvector
double maxQV(const QVector<double> &vec)
{
	double curr = vec[0];
	for (double o : vec)
	{
		if (o > curr)
			curr = o;
	}
	return curr;
}

//wrapper for info message generation
QString numInfo(double num)
{
	return QString("Graph  (generated ") + QString::number(num, 'g', 10) + QString(" points)");
}

//starter slot for actual numerical calculations
void DifferentialEquationDemo::launchMain()
{
	//Disable user input
	enableInputGeneral(false);

	//Prepare object containers for new thread
	QThread* thread = new QThread;
	Solver* solver = nullptr;

	//Check input fields to be non-empty
	if (!emptyCheckMain()) return;

	//Get data from input fields
	const double x_init = ui.x_init->text().toDouble();
	const double y_init = ui.y_init->text().toDouble();
	double h = ui.step_size->text().toDouble();
	const double x_to = ui.x_final->text().toDouble();
	if(ui.m_number->isChecked())
	{
		h = DTOI(h);
		h = (x_to - x_init) / h;
	}

	//Check data to be valid
	if (!dataCheckMain(x_init, y_init, h, x_to)) return;

	//Decide what numerical method to use and assign to solver pointer
	if (ui.choice_euler->isChecked())
	{
		solver = new Euler(x_init, y_init, h, x_to, e_19);
	}
	else if (ui.choice_eulerimp->isChecked())
	{
		solver = new Euler_Improved(x_init, y_init, h, x_to, e_19);
	}
	else if (ui.choice_rungek->isChecked())
	{
		solver = new RungeKutta(x_init, y_init, h, x_to, e_19);
	}
	else
	{
		qDebug() << "No choice for plot-type-main taken.";
		return;
	}

	//Now we have everything checked, so we can throw all data into new thread
	solver->moveToThread(thread);
	THREAD_WRAPPER(solver)
	connect(solver, SIGNAL(result(PlotDataObject)), this, SLOT(catchMain(PlotDataObject)));
	thread->start();
	//thus current main thread has nothing to do with calculations
}

//slot that will be activated on callback from numerical-calculations thread
void DifferentialEquationDemo::catchMain(PlotDataObject obj)
{
	// Clear previous graph data from viewer and local storage
	ui.customPlot->clearGraphs();
	this->x_num.clear();
	this->y_num.clear();
	this->x_aly.clear();
	this->y_aly.clear();
	
	// Create new graph
	ui.customPlot->addGraph();
	ui.customPlot->graph(0)->setData(obj.x, obj.y);
	ui.customPlot->xAxis->setLabel("X");
	ui.customPlot->yAxis->setLabel("Y");
	
	// Apply design for graph
	QPen bluePen;
	bluePen.setColor(QColor(30, 40, 255, 150));
	bluePen.setWidthF(5);
	ui.customPlot->graph(0)->setPen(bluePen);
	ui.customPlot->graph(0)->setName("Numeric Plot");

	// Set axes ranges, so we see all data conveniently:
	const double range = obj.x[obj.x.size() - 1] - obj.x[0];
	ui.customPlot->xAxis->setRange(obj.x[0] - BORDER_MARGIN * range, obj.x[obj.x.size() - 1] + BORDER_MARGIN * range);
	// See if analytical plot is needed or not
	if(ui.isAnalytical->isChecked())
	{
		//if it is we will need to save plotting data because on this branch it is possible to work with error data
		this->x_num = obj.x;
		this->y_num = obj.y;
		//request starting analytical calculation
		launchSimple();
	}else
	{
		//if analytical plot is not needed then we can already set the data on graph viewer
		const double minW = minQV(obj.y);
		const double maxW = maxQV(obj.y);
		const double range = maxW - minW;
		ui.customPlot->yAxis->setRange(minW - BORDER_MARGIN * range, maxW + BORDER_MARGIN * range);
		enableInputGeneral(true);
		ui.customPlot->replot();
	}
	//update information header about number of generated points
	ui.Graph->setTitle(numInfo(obj.x.size()));
}

//starter function for analytical calculations
void DifferentialEquationDemo::launchSimple()
{
	//get the values from input fields again
	const double x_init = ui.x_init->text().toDouble();
	const double y_init = ui.y_init->text().toDouble();
	const double x_to = ui.x_final->text().toDouble();
	double h = ui.step_size->text().toDouble();
	if (ui.m_number->isChecked())
	{
		h = DTOI(h);
		h = (x_to - x_init) / h;
	}

	//prepare object containers for new thread
	QThread* thread = new QThread;
	SimplePlot* simplot = new SimplePlot(x_init, y_init, h, x_to, e_19_sol, e_19_const);
	
	//create signal bindings
	simplot->moveToThread(thread);
	THREAD_WRAPPER(simplot)
	connect(simplot, SIGNAL(result(PlotDataObject)), this, SLOT(catchSimple(PlotDataObject)));
	//start prepared thread
	thread->start();
}

//slot that will be activated on callback from analytical-calculations thread
void DifferentialEquationDemo::catchSimple(PlotDataObject obj)
{
	//put resulting data into viewer
	this->x_aly = obj.x;
	this->y_aly = obj.y;
	ui.customPlot->addGraph();
	ui.customPlot->graph(1)->setData(obj.x, obj.y);
	
	//set design for chart
	QPen redDotPen;
	redDotPen.setColor(QColor(255, 0, 0, 150));
	redDotPen.setStyle(Qt::DotLine);
	redDotPen.setWidthF(4);
	ui.customPlot->graph(1)->setPen(redDotPen);
	ui.customPlot->graph(1)->setName("Analytic Plot");
	
	// set axes ranges, so we see all data conveniently:
	const double minW = minQV(obj.y);
	const double maxW = maxQV(obj.y);
	const double range = maxW - minW;
	ui.customPlot->yAxis->setRange(minW - BORDER_MARGIN * range, maxW + BORDER_MARGIN * range);
	
	//re-enable the input and update window
	enableInputGeneral(true);
	ui.customPlot->replot();
}

//starter function for single error calculation from currently drawn graph 
void DifferentialEquationDemo::launchError()
{
	//check storing arrays to contain values
	if(x_num.size()>0&&y_num.size()>0&&x_aly.size()>0&&y_aly.size())
	{
	}else
	{
		QMessageBox::warning(this, "Insufficient information", "Please plot graph with analytical solution to work with error");
		return;
	}

	//decide what error graph is needed to be plotted
	//launch according function
	ui.s_result->setText(QString("..."));
	if (ui.s_err_mse->isChecked())
	{
		ui.s_label->setText(QString("MSE:"));
		launchSE(0);
	}
	else if (ui.s_err_mae->isChecked())
	{
		ui.s_label->setText(QString("MAE:"));
		launchSE(1);
	}
	else if (ui.s_err_mxe->isChecked())
	{
		ui.s_label->setText(QString("MXE:"));
		launchSE(2);
	}else
	{
		qDebug() << "No choice for plot-type-error taken.";
	}
}

//starter function for simple error calculation, i.e. single value from current graph
void DifferentialEquationDemo::launchSE(const int &type)
{
	//check staring arrays to contain values
	if (x_num.size()>0 && y_num.size()>0 && x_aly.size()>0 && y_aly.size())
	{
	}
	else
	{
		QMessageBox::warning(this, "Insufficient information", "Please plot graph with analytical solution to work with error");
		return;
	}
	//create new object for the thread
	QThread* thread = new QThread;
	ErrorCalc* calc;
	switch(type)
	{
		case 0: 
		case 1: 
		case 2: calc = new ErrorCalc(y_num, y_aly, type); break;
		case 3: calc = new ErrorCalc(y_num, y_aly, x_num, type); break;
		default:qDebug("Non existant error calculation type"); return;
	}
	calc->moveToThread(thread);
	THREAD_WRAPPER(calc)
	switch (type)
	{
		case 0:
		case 1:
		case 2: connect(calc, SIGNAL(result(double)), this, SLOT(catchSE(double))); break;
		case 3: connect(calc, SIGNAL(result(PlotDataObject)), this, SLOT(catchSE(PlotDataObject))); break;
	}
	thread->start();
}

//slot that will be called when error calculating thread will return value
void DifferentialEquationDemo::catchSE(double val)
{
	ui.s_result->setText(QString::number(val));
}

//slot that will be called when error calculating thread will return error graph
void DifferentialEquationDemo::catchSE(PlotDataObject obj)
{
	////clear previous plot data
	ui.customPlot->clearGraphs();
	this->x_num.clear();
	this->y_num.clear();
	this->x_aly.clear();
	this->y_aly.clear();

	//create new graph
	ui.customPlot->addGraph();
	//assign values
	ui.customPlot->graph(0)->setData(obj.x, obj.y);
	
	//create style for that graph
	QPen redPen;
	redPen.setColor(QColor(255, 0, 0, 150));
	redPen.setStyle(Qt::SolidLine);
	redPen.setWidthF(4);
	//assign style
	ui.customPlot->graph(0)->setPen(redPen);
	
	//create brush that will fill space under the graph
	QBrush crossedFill(Qt::CrossPattern);
	crossedFill.setColor(redPen.color());
	//assign the brush for graph
	ui.customPlot->graph(0)->setBrush(crossedFill);
	
	//set labels for axes and name the graph
	ui.customPlot->graph(0)->setName("Error");
	ui.customPlot->xAxis->setLabel("X");
	ui.customPlot->yAxis->setLabel("Error Value");
	
	// set axes ranges, so graphs are visible
	const double rangeX = obj.x[obj.x.size() - 1] - obj.x[0];
	ui.customPlot->xAxis->setRange(obj.x[0] - BORDER_MARGIN * rangeX, obj.x[obj.x.size() - 1] + BORDER_MARGIN * rangeX);
	const double minW = minQV(obj.y);
	const double maxW = maxQV(obj.y);
	const double rangeY = maxW - minW;
	ui.customPlot->yAxis->setRange(minW - BORDER_MARGIN * rangeY, maxW + BORDER_MARGIN * rangeY);
	
	//re-enable the input and update window
	enableInputGeneral(true);
	ui.customPlot->replot();
}

//draw or hide on-graph legend that shows name of each graph in window
void DifferentialEquationDemo::legendSwitch()
{
	ui.customPlot->legend->setVisible(ui.legend->isChecked());
	ui.customPlot->replot();
}

//special wrapper for simple error graph building, it is forwarded to function with according value
//written separately to accessible by button signal
void DifferentialEquationDemo::launchSE3()
{
	launchSE(3);
}

//check if given value is in bounds of given list
double inRange(const QVector<double> &list, const double &val)
{
	return val >= list[0] && val <= list[list.size()-1];
}

//find index of value that is just below given in a list
double findLowerIndex(const QVector<double> &list, const double &val)
{
	for (int i = 1; i < list.size(); i++)
	{
		if(val>=list[i-1] && val<=list[i])
		{
			return i - 1;
		}
	}
	return -1;
}

//given two points x0, x2 with values y0,y2 find y1 that corresponds to given x1 by linearly interpolating on that line
double interpolateLinearY(double x_left, double y_left, double x_right, double y_right, double x_get) {
	double L = x_right - x_left;
	double H = y_right - y_left;
	double L1 = x_get - x_left;
	double H1 = H*L1 / L;
	return y_left + H1;
}

//starter function to launch calculations of value y1 from x1 given by user
void DifferentialEquationDemo::liCalc()
{
	//check that there is a value
	if (!ui.li_val->text().isEmpty())
	{
		enableInputGeneral(false);
		//load the value into double
		double val = ui.li_val->text().toDouble();
		int index = -1;
		double result;

		//check that numerical function containers are non-empty
		if (x_num.size() > 1 && y_num.size() > 1)
		{
			//check if given value is in range
			if (inRange(x_num, val))
			{
				//find the value by given x1
				index = findLowerIndex(x_num, val);
				result = interpolateLinearY(x_num[index], y_num[index], x_num[index + 1], y_num[index + 1], val);
				ui.li_res_num->setText(QString::number(result));
			}
			else
			{
				QMessageBox::information(this, "Interpolation Stopped", "X is out of range for Numerical graph");
			}
		}
		else
		{
			QMessageBox::information(this, "Interpolation Stopped", "No graph for Numerical interpolation");
		}
		
		//check that analytical function containers are non-empty
		if(x_aly.size()>1&&y_aly.size()>1)
		{
			//check if given value is in range
			if (inRange(x_aly, val))
			{
				//find the value by given x1
				index = findLowerIndex(x_aly, val);
				result = interpolateLinearY(x_aly[index], y_aly[index], x_aly[index + 1], y_aly[index + 1], val);
				ui.li_res_aly->setText(QString::number(result));
			}
			else
			{
				QMessageBox::information(this, "Interpolation Stopped", "X is out of range for Analytical graph");
			}
		}
		else
		{
			QMessageBox::information(this, "Interpolation Stopped", "No graph for Analytical interpolation");
		}
		enableInputGeneral(true);
	}
	else
	{
		QMessageBox::warning(this, "Error", "Empty Interpolation Value");
	}
}

//complementary handler for UI
void DifferentialEquationDemo::mSizeHandler()
{
	ui.m_number->setChecked(false);
	ui.m_size->setChecked(true);
}

//complementary handler for UI
void DifferentialEquationDemo::mNumbHandler()
{
	ui.m_size->setChecked(false);
	ui.m_number->setChecked(true);
}

//check that error plotting settings 
bool DifferentialEquationDemo::emptyCheckErr()
{
	if (ui.e_from->text() == "" || ui.e_to->text() == "" || ui.e_num->text() == "")
	{
		QMessageBox::warning(this, "Invalid input", "Empty fields present (Error graph)");
		enableInputGeneral(true);
		return false;
	}
	return true;
}

//check values to be valid for calculations
bool DifferentialEquationDemo::dataCheckErr(const double& e_from, const double& e_to, const double& e_num)
{
	//number of samples should be greater than 0
	if(DTOI(e_num)>0)
	{
		return true;
	}
	else
	{
		QMessageBox::warning(this, "Invalid input", "Can not plot graph with such inputs");
		enableInputGeneral(true);
		return false;
	}

}

//starter function for plotting error graph (error to step-size)
void DifferentialEquationDemo::launchErrGraph()
{
	//Disable user input
	enableInputGeneral(false);

	//Prepare object containers for new thread
	QThread* thread = new QThread;
	Solver* solver = nullptr;

	//Check input fields to be non-empty
	if (!emptyCheckErr()) return;

	//Get data from input fields
	const double e_from = ui.e_from->text().toDouble();
	const double e_to = ui.e_to->text().toDouble();
	const double e_num = ui.e_num->text().toDouble();

	//Check data to be valid
	if (!dataCheckErr(e_from, e_to, e_num)) return;

	//Check input fields to be non-empty
	if (!emptyCheckMain_E()) return;

	//Get data from input fields
	const double x_init = ui.x_init->text().toDouble();
	const double y_init = ui.y_init->text().toDouble();
	const double x_to = ui.x_final->text().toDouble();

	//Check data to be valid
	if (!dataCheckMain_E(x_init, y_init, x_to)) return;
	
	int method;
	//Decide what method for error calculation to use
	if (ui.s_err_mse->isChecked()) //for mean-squared-error method
	{
		method = 0;
	}else if(ui.s_err_mae->isChecked()) //for mean-absolute-error method
	{
		method = 1;
	}else if(ui.s_err_mxe->isChecked()) //for maximum-error method
	{
		method = 2;
	}else
	{
		qDebug() << "No choice for plot-type-error-method taken.";
		return;
	}

	//Decide what numerical method to use
	if (ui.choice_euler->isChecked())
	{
		solver = new Euler(x_init, y_init, x_to, e_from, e_to, e_num, e_19, e_19_sol, e_19_const, method);
	}
	else if (ui.choice_eulerimp->isChecked())
	{
		solver = new Euler_Improved(x_init, y_init, x_to, e_from, e_to, e_num, e_19, e_19_sol, e_19_const, method);
	}
	else if (ui.choice_rungek->isChecked())
	{
		solver = new RungeKutta(x_init, y_init, x_to, e_from, e_to, e_num, e_19, e_19_sol, e_19_const, method);
	}
	else
	{
		qDebug() << "No choice for plot-type-error-graph taken.";
		return;
	}

	//Now we have everything checked, so we can throw all data into new thread
	solver->moveToThread(thread);
	connect(solver, SIGNAL(error(QString)), this, SLOT(errorWarning(QString)));
	connect(thread, SIGNAL(started()), solver, SLOT(err_proc()));
	connect(solver, SIGNAL(finished()), thread, SLOT(quit()));
	connect(solver, SIGNAL(finished()), solver, SLOT(deleteLater()));
	connect(thread, SIGNAL(finished()), thread, SLOT(deleteLater()));
	connect(solver, SIGNAL(result(PlotDataObject)), this, SLOT(catchErrGraph(PlotDataObject)));
	thread->start();
	//thus current main thread has nothing to do with calculations
}

//slot that will be called when error calculating thread returns value
void DifferentialEquationDemo::catchErrGraph(PlotDataObject obj)
{
	//clear previous plot data
	ui.customPlot->clearGraphs();
	this->x_num.clear();
	this->y_num.clear();
	this->x_aly.clear();
	this->y_aly.clear();

	//create new graph
	ui.customPlot->addGraph();
	//assign values
	ui.customPlot->graph(0)->setData(obj.x, obj.y);
	
	//create style for that graph
	QPen redPen;
	redPen.setColor(QColor(255, 0, 0, 150));
	redPen.setStyle(Qt::SolidLine);
	redPen.setWidthF(4);
	//assign the style
	ui.customPlot->graph(0)->setPen(redPen);
	//create brush that will fill space under the graph
	QBrush crossedFill(Qt::CrossPattern);
	crossedFill.setColor(redPen.color());
	//assign the brush for graph
	ui.customPlot->graph(0)->setBrush(crossedFill);
	
	//set labels for axes and name the graph
	ui.customPlot->graph(0)->setName("Error-to-Step-Size");
	ui.customPlot->xAxis->setLabel("Step Size");
	ui.customPlot->yAxis->setLabel("Error Value");
	
	// set axes ranges, so graphs are visible
	const double rangeX = obj.x[obj.x.size() - 1] - obj.x[0];
	ui.customPlot->xAxis->setRange(obj.x[0] - BORDER_MARGIN * rangeX, obj.x[obj.x.size() - 1] + BORDER_MARGIN * rangeX);
	const double minW = minQV(obj.y);
	const double maxW = maxQV(obj.y);
	const double rangeY = maxW - minW;
	ui.customPlot->yAxis->setRange(minW - BORDER_MARGIN * rangeY, maxW + BORDER_MARGIN * rangeY);
	
	//re-enable the input and update window
	enableInputGeneral(true);
	ui.customPlot->replot();
}

//separate function to convert double to int, in current context it is used to enforce local conversion policy
inline int DifferentialEquationDemo::DTOI(const double &val)
{
	return static_cast<int>(std::floor(val));
}
