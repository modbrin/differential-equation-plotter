#include "DifferentialEquationDemo.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	DifferentialEquationDemo w;
	w.show();
	return a.exec();
}
