#ifndef _GLWIDGET_H_
#define _GLWIDGET_H_

#include <CGAL/glu.h>
#include <QtOpenGL>
#include <QGLWidget>

#include "ColorRamp/ColorRamp.h"
#include "ThreeDimensionalShape.h"
#include "LinearAlgebra/Wm4Vector.h"
#include "LinearAlgebra/Wm4Matrix.h"

#define BUFFER_LENGTH 512
enum materialtype
{
	SILVER, 
	GOLD, 
	JADE, 
	LIGHT_BLUE, 
	EMERALD, 
	POLISHED_SILVER, 
	CHROME, 
	COPPER, 
	POLISHED_GOLD, 
	PEWTER, 
	OBSIDIAN, 
	BLACK_PLASTIC, 
	POLISHED_BRONZE, 
	POLISHED_COPPER, 
	PEARL, 
	RUBY, 
	TURQUOISE, 
	BRASS
};


class GLWidget : public QGLWidget
{
	//tobedone

	Q_OBJECT
public:
	GLWidget(QWidget * parent = 0);
	~GLWidget();

	void init();
	void set3DShape(ThreeDimensionalShape * pThreeDimensionalShape);
	//Zhiyang Modified
	void set3DSelected_Pole(vector<vector<double>> & pselect_poles);

	void updateRender();

	int xRotation() const {return xRot;}
	int yRotation() const {return yRot;}
	int zRotation() const {return zRot;}

public slots:
	void setXYRotationChanged(int xAngle, int yAngle);
	void setXZRotationChanged(int xAngle, int zAngle);
	void setRotationChanged(int xAngle, int yAngle, int zAngle);
	void setXRotation(int angle);
	void setYRotation(int angle);
	void setZRotation(int angle);
	void setZoom(double factor);
	void setTranslate(double dx, double dy);

private slots:

signals:
	void xRotationChanged(int angle);
	void yRotationChanged(int angle);
	void zRotationChanged(int angle);
	void ZoomChanged(double factor);
	void TranslateChanged(double dx, double dy);

protected:
	void initializeGL();
	void paintGL();
	void resizeGL(int width, int height);
	void mousePressEvent(QMouseEvent * event);
	void mouseMoveEvent(QMouseEvent * event);
	void mouseReleaseEvent(QMouseEvent * event);
	void wheelEvent(QWheelEvent * event);
	void mouseDoubleClickEvent(QMouseEvent * event);

private:
	void normalizeAngle(int * angle);
	void setMaterial(const int str);
	void process_select(int xPos, int yPos);

	// Rendering
	void drawScene(GLdouble dx = 0, GLdouble dy = 0, GLdouble dz = 0, GLdouble angle = 0);
	Vector3d HSVToRGB(Vector3d in);

	// Rendering Mesh
	void gl_draw_faces(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_edges(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_vertices(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_facet_normals(ThreeDimensionalShape * pThreeDimensionalShape);

	// Rendering slab simplification

	void gl_draw_slab_vertex(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_slab_edge(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_slab_face(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_slab_envelop(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_slab_sphere(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_mean_square_error(ThreeDimensionalShape * pThreeDimensionalShape);

	void gl_draw_fake_boundary(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_boundary(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_nonmanifold(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_saved_vertex(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_stability_ratio(ThreeDimensionalShape * pThreeDimensionalShape);

	void processLeftButton(int xPos, int yPos);
	void gl_draw_selected_vertices();

	// rendered by Zhiyang ...
	void gl_draw_select_poles(vector<vector<double>> select_poles);




	// Rendering Envelope Mesh
	//void gl_draw_envelope_faces(ThreeDimensionalShape * pThreeDimensionalShape);
	//void gl_draw_envelope_edges(ThreeDimensionalShape * pThreeDimensionalShape);
	//void gl_draw_envelope_vertices(ThreeDimensionalShape * pThreeDimensionalShape);

	// Rendering DT and Voronoi
	void gl_draw_input_dt(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_input_voronoi(ThreeDimensionalShape * pThreeDimensionalShape);

	// Rendering Medial Axis
	void gl_draw_ma_vertices();
	void gl_draw_ma_vertices_bplist();
	void gl_draw_ma_edges();
	void gl_draw_ma_faces();
	void gl_draw_ma_footpoints();
	void gl_draw_ma_spheres();
	void gl_draw_ma_envelopes();
	
	void gl_adjust_view(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_adjust_ma(ThreeDimensionalShape * pThreeDimensionalShape);
	void gl_draw_facet(Facet_handle pFacet, double factor = 1.);
	void gl_draw_headbar();
	void gl_draw_colorbar();

public:
	void SaveViewpoint(const char * fout);
	void LoadViewpoint(const char * fin);

public:
	CColorRamp m_ColorRamp;

	int xRot;
	int yRot;
	int zRot;
	double xTrans, yTrans;

	int red;
	int green;
	int blue;

	QPoint lastPos;
	QMenu * m_PopMenu;

	bool m_smoothshading;
	bool m_antialiasing;
	bool m_cullface;
	bool m_setcullface;
	
	double scalefactor;
	double oldscale;
	double oldtranslation[3];

	bool m_firstview;
	
	bool m_display_faces;
	bool m_display_edges;
	bool m_display_vertices;
	bool m_display_select_pole;
	bool m_display_facet_normals;

	bool m_display_envelope_faces;
	bool m_display_envelope_edges;
	bool m_display_envelope_vertices;

	bool m_display_input_dt;
	bool m_display_input_voronoi;

	bool m_display_ma_vertices;
	bool m_display_ma_vertices_bplist;
	bool m_display_ma_edges;
	bool m_display_ma_faces;
	bool m_display_ma_footpoints;
	bool m_display_ma_spheres;
	bool m_display_ma_envelopes;
	
	bool m_display_headbar;
	bool m_display_axes;
	bool m_display_colorbar;

	bool m_display_reverse_orientation;

	bool m_display_logscale;

	// trigger for qem mesh
	bool m_display_qem_vertex;
	bool m_display_qem_edge;
	bool m_display_qem_face;

	// trigger for ma qem mesh
	bool m_display_ma_qem_vertex;
	bool m_display_ma_qem_edge;
	bool m_display_ma_qem_envelop;
	bool m_display_ma_qem_sphere;
	bool m_display_centroid;

	bool m_display_fake_boundary;
	unsigned m_display_boundary;
	bool m_display_nonmanifold;
	bool m_display_saved_vertex;

	bool m_display_stability_ratio;

	// trigger for slab simplification
	bool m_display_slab_vertex;
	bool m_display_slab_edge;
	bool m_display_slab_faces;
	bool m_display_slab_envelop;
	bool m_display_slab_sphere;
	bool m_display_zoom_ma;
	bool m_display_mean_squre_error;

	// trigger for sphere mesh
	bool m_display_sphere_vertex;
	bool m_display_sphere_edge;
	bool m_display_sphere_mesh_sphere;
	bool m_display_sphere_mesh_envelop;

	bool m_select_pick_mode;

	unsigned int m_logscale_highvalue; // low value is 0. Thus, mapping from 1(10^lv) ~ 10^hv

	double m_Colorbar_MaxValue;
	double m_Colorbar_MinValue;

	double m_offsetfactor;

	GLUquadricObj *spherequadric;
	GLUquadricObj * conequadric;

	ThreeDimensionalShape * m_pThreeDimensionalShape;
	vector<vector<double>> select_poles;


	//选择缓冲区的空间
	GLuint selectBuffer[BUFFER_LENGTH];
	GLint hits;
};
#endif //_GLWIDGET_H_