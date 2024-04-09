#ifndef MEDIALAXISSIMPLIFICATION3D_H
#define MEDIALAXISSIMPLIFICATION3D_H

#include <QtGui>

#include <QtGui>
#include <QtOpenGL>

#include "ThreeDimensionalShape.h"
#include "GLWidget.h"
#include "ui_medial_axis.h"
#include "PsRender/PsRender.h"

class MedialAxisSimplification3D : public QMainWindow
{
	Q_OBJECT

public:
	MedialAxisSimplification3D(QWidget *parent = 0, Qt::WindowFlags flags = 0);
	~MedialAxisSimplification3D();

private slots:

	void openmeshfile();
	void openmeshfile2();
	void showface(bool _val);
	void showedge(bool _val);
	void showvertex(bool _val);
	void showfacetnormal(bool _val);
	void showdt(bool _val);
	void showvoronoi(bool _val);
	void showmat(bool _val);
	void showsphere(bool _val);
	void showenvelop(bool _val);

	void determineVertexType();

	bool importMA(QString name);
	void importVP(QString viewpoint);

	void importMA();
	void exportMA();
	void computeMA();

	void exportSimplifiedMA();

	// slab simplification
	void initialslab();
	void showslabvertex(bool _val);
	void showslabedge(bool _val);
	void showslabface(bool _val);
	void showslabenvelop(bool _val);
	void showslabsphere(bool _val);
	void simplifySlabByOne();
	void simplifySlab();
	void showmeansquareerror(bool _val);

	// slab method parameters
	void initialize();
	void zoomma(bool _val);

	void showstabilityratio(bool _val);

	void setratiofactor();
	void loadconfigurefile();

	void saveEPS();
	void loadViewPoint();
	void saveViewPoint();
	void setCheckReverseOrientation(bool _val);
	void setCheckHeadBar(bool _val);
	void setCheckColorBar(bool _val);
	void setSmoothRendering(bool _val);

public slots:
	void updateviews();

private:
	
	void createActions();

	QWidget * m_pCentralWidget;
	QGridLayout * m_pCentralLayout;
	QScrollArea * m_pCentralArea;
	GLWidget * m_pGLWidget;
	bool slab_initial;
	ThreeDimensionalShape * m_pThreeDimensionalShape;
	double k;

	//
	QString m_offFileName;
	bool m_isSimplified;

private:
	Ui::medial_axisClass ui;
};

#endif // MEDIALAXISSIMPLIFICATION3D_H
