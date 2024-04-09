//#define BOOST_PARAMETER_MAX_ARITY 12
#include "medialaxissimplification3d.h"
#include <QtGui>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	MedialAxisSimplification3D w;
	w.showMaximized(); // ×î´ó»¯
	// this line just for testing commiting
	return a.exec();
}
