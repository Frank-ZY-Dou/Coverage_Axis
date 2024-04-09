#include "GLWidget.h"

GLWidget::GLWidget(QWidget * parent) : QGLWidget(parent)
{
	hits = 0;
	xRot = yRot = zRot = 0;
	xTrans = yTrans = 0;
	m_smoothshading = false;
	m_antialiasing = true;
	m_cullface = false;
	m_setcullface = false;
	scalefactor = 1.;
	oldscale = 1.;
	oldtranslation[0] = oldtranslation[1] = oldtranslation[2] = 0.;
	m_display_faces = false;
	m_display_edges = false;
	m_display_vertices = false;
	m_display_select_pole = false;
	m_display_envelope_faces = false;
	m_display_envelope_edges = false;
	m_display_envelope_vertices = false;

	m_display_input_dt = false;
	m_display_input_voronoi = false;

	m_display_ma_vertices = false;
	m_display_ma_vertices_bplist = false;
	m_display_ma_edges = false;
	m_display_ma_faces = false;
	m_display_ma_footpoints = false;
	m_display_ma_spheres = false;
	m_display_ma_envelopes = false;
	
	m_display_headbar = false;
	m_display_axes = false;
	m_display_colorbar = false;

	m_display_facet_normals = false;

	// triggers for qem mesh
	m_display_qem_vertex = false;
	m_display_qem_edge = false;
	m_display_qem_face = false;

	// triggers for ma qem mesh
	m_display_ma_qem_vertex = false;
	m_display_ma_qem_edge = false;
	m_display_ma_qem_envelop = false;
	m_display_ma_qem_sphere = false;
	m_display_centroid = false;

	// triggers for slab simplification
	m_display_slab_vertex = false;
	m_display_slab_edge = false;
	m_display_slab_faces = false;
	m_display_slab_envelop = false;
	m_display_slab_sphere = false;
	m_display_zoom_ma = false;
	m_display_mean_squre_error = false;

	m_select_pick_mode = false;

	m_display_fake_boundary = false;
	m_display_boundary = 0;
	m_display_nonmanifold = false;
	m_display_saved_vertex = false;
	m_display_stability_ratio = false;

	// triggers for sphere mesh
	m_display_sphere_vertex = false;
	m_display_sphere_edge = false;
	m_display_sphere_mesh_envelop = false;
	m_display_sphere_mesh_sphere = false;

	m_display_reverse_orientation = false;

	m_display_logscale = false;
	m_logscale_highvalue = 3;

	red = 159;
	green = 121;
	blue = 238;

	m_Colorbar_MaxValue = 1.;
	m_Colorbar_MinValue = 1e-8;

	m_pThreeDimensionalShape = NULL;
	m_ColorRamp.BuildRainbow();

	m_firstview = true;

	m_offsetfactor = 0.9999999;
	
	spherequadric = gluNewQuadric();
	gluQuadricDrawStyle(spherequadric, GLU_FILL );

	conequadric = gluNewQuadric();
	gluQuadricDrawStyle(conequadric, GLU_FILL );
}

GLWidget::~GLWidget()
{
	gluDeleteQuadric(spherequadric);
	gluDeleteQuadric(conequadric);
	makeCurrent();
}

void GLWidget::init()
{
	scalefactor = 1.;
}

void GLWidget::set3DShape(ThreeDimensionalShape * pThreeDimensionalShape)
{
	m_pThreeDimensionalShape = pThreeDimensionalShape;
	updateGL();
}
void GLWidget::set3DSelected_Pole(vector<vector<double>> & pselect_poles) {
	select_poles = pselect_poles;
	cout << "updating..." << endl;
	updateGL();
}

void GLWidget::setXYRotationChanged(int xAngle, int yAngle)
{
	normalizeAngle(&xAngle);
	normalizeAngle(&yAngle);

	if(xAngle != xRot)
	{
		xRot = xAngle;
		emit xRotationChanged(xAngle);
	}
	if(yAngle != yRot)
	{
		yRot = yAngle;
		emit yRotationChanged(yAngle);
	}
	updateGL();
}

void GLWidget::setXZRotationChanged(int xAngle, int zAngle)
{
	normalizeAngle(&xAngle);
	normalizeAngle(&zAngle);

	if(xAngle != xRot)
	{
		xRot = xAngle;
		emit xRotationChanged(xAngle);
	}
	if(zAngle != zRot)
	{
		zRot = zAngle;
		emit zRotationChanged(zAngle);
	}
	updateGL();
}

void GLWidget::setRotationChanged(int xAngle, int yAngle, int zAngle)
{
	normalizeAngle(&xAngle);
	normalizeAngle(&yAngle);
	normalizeAngle(&zAngle);

	if(xAngle != xRot)
	{
		xRot = xAngle;
		emit xRotationChanged(xAngle);
	}
	if(yAngle != yRot)
	{
		yRot = yAngle;
		emit yRotationChanged(yAngle);
	}
	if(zAngle != zRot)
	{
		zRot = zAngle;
		emit zRotationChanged(zAngle);
	}
	updateGL();
}

void GLWidget::setXRotation(int angle)
{
	normalizeAngle(&angle);
	if(angle != xRot)
	{
		xRot = angle;
		emit xRotationChanged(angle);
		updateGL();
	}
}

void GLWidget::setYRotation(int angle)
{
	normalizeAngle(&angle);
	if(angle != yRot)
	{
		yRot = angle;
		emit yRotationChanged(angle);
		updateGL();
	}
}

void GLWidget::setZRotation(int angle)
{
	normalizeAngle(&angle);
	if(angle != zRot)
	{
		zRot = angle;
		emit zRotationChanged(angle);
		updateGL();
	}
}

void GLWidget::setZoom(double factor)
{
	scalefactor += factor / 10.;
	emit ZoomChanged(factor);
	updateGL();
}

void GLWidget::setTranslate(double dx, double dy)
{
	xTrans += dx;
	yTrans += dy;
	emit TranslateChanged(dx,dy);
	updateGL();
}

void GLWidget::initializeGL()
{
	static const GLfloat lightPos0[4] = {5.0f, 5.0f, 10.0f, 1.0f};
	//static const GLfloat lightPos0[4] = {-3.0,-3.0,3.0,1.0};
	static const GLfloat lightPos1[4] = {-5.0f, -5.0f, -10.0f, 1.0f};
	static const GLfloat lightPos2[4] = {5.0f, -5.0f, 0.0f, 1.0f};
	static const GLfloat lightPos3[4] = {-5.0f, 5.0f, 0.0f, 1.0f};
	glEnable(GL_DEPTH_TEST);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
	glClearColor(1.,1.,1.,1.);

	GLfloat light0_ambient[]={0.3,0.3,0.3,1.0}; 
	GLfloat light0_diffuse[]={0.5,0.5,0.5,1.0};
	GLfloat light0_specular[]={0.5,0.5,0.5,1.0}; 
	GLfloat light0_direction[]={1.0,1.0,-1.0}; 
	
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHT2);
	glEnable(GL_LIGHT3);

	glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
	//glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
	//glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
	//glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);

	glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);
	glLightfv(GL_LIGHT2, GL_POSITION, lightPos2);
	glLightfv(GL_LIGHT3, GL_POSITION, lightPos3);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);
	float MatAmbientBack[] = {0.0f, 1.0f, 0.0f, 1.0f};
	glMaterialfv(GL_BACK, GL_AMBIENT, MatAmbientBack);
	setMaterial(COPPER);
	glEnable(GL_NORMALIZE);
}

void GLWidget::paintGL()
{

	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glLoadIdentity();
	glTranslated( xTrans, yTrans, -20.0 );

	glRotated(xRot / 16.0, 1.0, 0.0, 0.0);
	glRotated(yRot / 16.0, 0.0, 1.0, 0.0);
	glRotated(zRot / 16.0, 0.0, 0.0, 1.0);

	drawScene();
}

// Replaces gluPerspective. Sets the frustum to perspective mode.
// fovY     - Field of vision in degrees in the y direction
// aspect   - Aspect ratio of the viewport
// zNear    - The near clipping distance
// zFar     - The far clipping distance
 
void perspectiveGL( GLdouble fovY, GLdouble aspect, GLdouble zNear, GLdouble zFar )
{
const GLdouble pi = 3.1415926535897932384626433832795;
GLdouble fW, fH;
fH = tan( fovY / 360 * pi ) * zNear;
fW = fH * aspect;
glFrustum( -fW, fW, -fH, fH, zNear, zFar );
}

void GLWidget::resizeGL(int width, int height)
{
	int side = qMin(width, height);
	if(side > 0)
	{
		glViewport(0,0,width,height);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		//gluPerspective(45.0, (GLfloat) width/height, 1.0, 100.0);
		perspectiveGL(45.0, (GLfloat) width/height, 1.0, 100.0);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
}

void GLWidget::mousePressEvent(QMouseEvent * event)
{
	makeCurrent();
	lastPos = event->pos();

	if (m_select_pick_mode)
	{
		processLeftButton(lastPos.x(), lastPos.y());
	}
}

void GLWidget::mouseMoveEvent(QMouseEvent * event)
{
	int dx = event->x() - lastPos.x();
	int dy = event->y() - lastPos.y();

	//// !Left !Right Mid --> Translate
	//if( !(event->buttons() & Qt::LeftButton) && !(event->buttons() & Qt::RightButton) && (event->buttons() & Qt::MidButton))
	//	setTranslate(dx / 16.0, -dy / 16.0);
	//// Left !Right !Mid --> X Rotation
	//if( (event->buttons() & Qt::LeftButton) && !(event->buttons() & Qt::RightButton) && !(event->buttons() & Qt::MidButton) )
	//	//setXRotation(xRot + 8*(dx));
	//	setRotationChanged(xRot + 8*(dx), yRot + 8*(dy), zRot + 8*(dx));
	//// Left !Right !Mid --> Y Rotation
	//if( (event->buttons() & Qt::LeftButton) && !(event->buttons() & Qt::RightButton) && !(event->buttons() & Qt::MidButton) )
	//	setYRotation(yRot + 8*(dy));

	//if (event->buttons() & Qt::LeftButton)
	//	setXYRotationChanged(xRot + 8*(dy), yRot + 8*(dx));
	//else if (event->buttons() & Qt::RightButton)
	//	setXZRotationChanged(xRot + 8*(dy), zRot + 8*(dx));
	

	// !Left !Right Mid --> Translate
	if( !(event->buttons() & Qt::LeftButton) && !(event->buttons() & Qt::RightButton) && (event->buttons() & Qt::MidButton))
		setTranslate(dx / 16.0, -dy / 16.0);
	// Left !Right !Mid --> X Rotation
	if( (event->buttons() & Qt::LeftButton) && !(event->buttons() & Qt::RightButton) && !(event->buttons() & Qt::MidButton) )
		setXRotation(xRot + 8*(dy+dx));
	// !Left Right !Mid --> Y Rotation
	if( !(event->buttons() & Qt::LeftButton) && (event->buttons() & Qt::RightButton) && !(event->buttons() & Qt::MidButton) )
		setYRotation(yRot + 8*(dy+dx));
	// Left Right !Mid  --> Z Rotation
	if( (event->buttons() & Qt::LeftButton) && (event->buttons() & Qt::RightButton) && !(event->buttons() & Qt::MidButton))
		setZRotation(zRot + 8*(dy+dx));
	// Left Right Mid --> Scaling
	if( (event->buttons() & Qt::LeftButton) && (event->buttons() & Qt::RightButton) && (event->buttons() & Qt::MidButton))
		setZoom((dy+dx) / 100.0);


	lastPos = event->pos();
	makeCurrent();
}

void GLWidget::mouseReleaseEvent(QMouseEvent * event)
{
	makeCurrent();
}

void GLWidget::wheelEvent(QWheelEvent * event)
{
	if( event->orientation() == Qt::Horizontal)
		setZoom(-event->delta() / 100.0);
	else
		setZoom(event->delta() / 100.0);
	makeCurrent();
}

void GLWidget::mouseDoubleClickEvent(QMouseEvent * event)
{
}

void GLWidget::normalizeAngle(int * angle)
{
	while(*angle < 0)
		*angle += 360*16;
	while(*angle > 360*16)
		*angle -= 360*16;
}

void GLWidget::setMaterial(const int str)
{
	float	ambient[] = {0.0f, 0.0f, 0.0f, 1.0f};
	float	diffuse[] = {0.0f, 0.0f, 0.0f, 1.0f};
	float	specular[] = {0.0f, 0.0f, 0.0f, 1.0f};
	float	emission[] = {0.3f, 0.3f, 0.3f, 1.0f};
	float shininess[] = {0.0f};

	if(str == SILVER)
	{
		// Ambient
		ambient[0] = 0.19225f;
		ambient[1] = 0.19225f;
		ambient[2] = 0.19225f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.50754f;
		diffuse[1] = 0.50754f;
		diffuse[2] = 0.50754f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.508273f;
		specular[1] = 0.508273f;
		specular[2] = 0.508273f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 51.2f;
	}

	else if(str == GOLD)
	{
		// Ambient
		ambient[0] = 0.24725f;
		ambient[1] = 0.1995f;
		ambient[2] = 0.0745f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.75164f;
		diffuse[1] = 0.60648f;
		diffuse[2] = 0.22648f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.928281f;
		specular[1] = 0.855802f;
		specular[2] = 0.666065f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 51.2f;
	}

	else if(str == JADE)
	{
		// Ambient
		ambient[0] = 0.135f;
		ambient[1] = 0.2225f;
		ambient[2] = 0.1575f;
		ambient[3] = 0.95f;
		// Diffuse
		diffuse[0] = 0.54f;
		diffuse[1] = 0.89f;
		diffuse[2] = 0.63f;
		diffuse[3] = 0.95f;
		// Specular
		specular[0] = 0.316228f;
		specular[1] = 0.316228f;
		specular[2] = 0.316228f;
		specular[3] = 0.95f;
		// Shininess
		shininess[0] = 12.8f;
	}

	else if(str == LIGHT_BLUE)
	{
		// Ambient
		ambient[0] = 0.0f;
		ambient[1] = 0.5f;
		ambient[2] = 0.75f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.0f;
		diffuse[1] = 0.5f;
		diffuse[2] = 1.0f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.75f;
		specular[1] = 0.75f;
		specular[2] = 0.75f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 64.0f;
	}

	else if(str == EMERALD)
	{
		// Ambient
		ambient[0] = 0.0215f;
		ambient[1] = 0.1745f;
		ambient[2] = 0.0215f;
		ambient[3] = 0.55f;
		// Diffuse
		diffuse[0] = 0.07568f;
		diffuse[1] = 0.61424f;
		diffuse[2] = 0.07568f;
		diffuse[3] = 0.55f;
		// Specular
		specular[0] = 0.633f;
		specular[1] = 0.727811f;
		specular[2] = 0.633f;
		specular[3] = 0.55f;
		// Shininess
		shininess[0] = 76.8f;
	}

	else if(str == POLISHED_SILVER)
	{
		// Ambient
		ambient[0] = 0.23125f;
		ambient[1] = 0.23125f;
		ambient[2] = 0.23125f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.2775f;
		diffuse[1] = 0.2775f;
		diffuse[2] = 0.2775f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.773911f;
		specular[1] = 0.773911f;
		specular[2] = 0.773911f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 89.6f;
	}

	else if(str == CHROME)
	{
		// Ambient
		ambient[0] = 0.25f;
		ambient[1] = 0.25f;
		ambient[2] = 0.25f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.4f;
		diffuse[1] = 0.4f;
		diffuse[2] = 0.4f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.774597f;
		specular[1] = 0.774597f;
		specular[2] = 0.774597f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 76.8f;
	}

	else if(str == COPPER)
	{
		// Ambient
		ambient[0] = 0.19125f;
		ambient[1] = 0.0735f;
		ambient[2] = 0.0225f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.7038f;
		diffuse[1] = 0.27048f;
		diffuse[2] = 0.0828f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.256777f;
		specular[1] = 0.137622f;
		specular[2] = 0.086014f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 12.8f;
	}

	else if(str == POLISHED_GOLD)
	{
		// Ambient
		ambient[0] = 0.24725f;
		ambient[1] = 0.2245f;
		ambient[2] = 0.0645f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.34615f;
		diffuse[1] = 0.3143f;
		diffuse[2] = 0.0903f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.797357f;
		specular[1] = 0.723991f;
		specular[2] = 0.208006f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 83.2f;
	}

	else if(str == PEWTER)
	{
		// Ambient
		ambient[0] = 0.105882f;
		ambient[1] = 0.058824f;
		ambient[2] = 0.113725f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.427451f;
		diffuse[1] = 0.470588f;
		diffuse[2] = 0.541176f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.333333f;
		specular[1] = 0.333333f;
		specular[2] = 0.521569f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 9.84615f;
	}

	else if(str == OBSIDIAN)
	{
		// Ambient
		ambient[0] = 0.05375f;
		ambient[1] = 0.05f;
		ambient[2] = 0.06625f;
		ambient[3] = 0.82f;
		// Diffuse
		diffuse[0] = 0.18275f;
		diffuse[1] = 0.17f;
		diffuse[2] = 0.22525f;
		diffuse[3] = 0.82f;
		// Specular
		specular[0] = 0.332741f;
		specular[1] = 0.328634f;
		specular[2] = 0.346435f;
		specular[3] = 0.82f;
		// Shininess
		shininess[0] = 38.4f;
	}

	else if(str == BLACK_PLASTIC)
	{
		// Ambient
		ambient[0] = 0.0f;
		ambient[1] = 0.0f;
		ambient[2] = 0.0f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.01f;
		diffuse[1] = 0.01f;
		diffuse[2] = 0.01f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.5f;
		specular[1] = 0.5f;
		specular[2] = 0.5f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 32.0f;
	}

	else if(str == POLISHED_BRONZE)
	{
		// Ambient
		ambient[0] = 0.25f;
		ambient[1] = 0.148f;
		ambient[2] = 0.006475f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.4f;
		diffuse[1] = 0.2368f;
		diffuse[2] = 0.1036f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.774597f;
		specular[1] = 0.458561f;
		specular[2] = 0.200621f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 76.8f;
	}


	else if(str == POLISHED_COPPER)
	{
		// Ambient
		ambient[0] = 0.2295f;
		ambient[1] = 0.08825f;
		ambient[2] = 0.0275f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.5508f;
		diffuse[1] = 0.2118f;
		diffuse[2] = 0.066f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.580594f;
		specular[1] = 0.223257f;
		specular[2] = 0.0695701f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 51.2f;
	}

	else if(str == PEARL)
	{
		// Ambient
		ambient[0] = 0.25f;
		ambient[1] = 0.20725f;
		ambient[2] = 0.20725f;
		ambient[3] = 0.922f;
		// Diffuse
		diffuse[0] = 1.0f;
		diffuse[1] = 0.829f;
		diffuse[2] = 0.829f;
		diffuse[3] = 0.922f;
		// Specular
		specular[0] = 0.296648f;
		specular[1] = 0.296648f;
		specular[2] = 0.296648f;
		specular[3] = 0.922f;
		// Shininess
		shininess[0] = 11.264f;
	}

	else if(str == RUBY)
	{
		// Ambient
		ambient[0] = 0.1745f;
		ambient[1] = 0.01175f;
		ambient[2] = 0.01175f;
		ambient[3] = 0.55f;
		// Diffuse
		diffuse[0] = 0.61424f;
		diffuse[1] = 0.04136f;
		diffuse[2] = 0.04136f;
		diffuse[3] = 0.55f;
		// Specular
		specular[0] = 0.727811f;
		specular[1] = 0.626959f;
		specular[2] = 0.626959f;
		specular[3] = 0.55f;
		// Shininess
		shininess[0] = 76.8f;
	}

	else if(str == TURQUOISE)
	{
		// Ambient
		ambient[0] = 0.1f;
		ambient[1] = 0.18725f;
		ambient[2] = 0.1745f;
		ambient[3] = 0.8f;
		// Diffuse
		diffuse[0] = 0.396f;
		diffuse[1] = 0.74151f;
		diffuse[2] = 0.69102f;
		diffuse[3] = 0.8f;
		// Specular
		specular[0] = 0.297254f;
		specular[1] = 0.30829f;
		specular[2] = 0.306678f;
		specular[3] = 0.8f;
		// Shininess
		shininess[0] = 12.8f;
	}

	else if(str == BRASS)
	{
		// Ambient
		ambient[0] = 0.329412f;
		ambient[1] = 0.223529f;
		ambient[2] = 0.027451f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.780392f;
		diffuse[1] = 0.268627f;
		diffuse[2] = 0.113725f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.992157f;
		specular[1] = 0.741176f;
		specular[2] = 0.807843f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 27.8974f;
	}

	else
	{
		;
	}
	// apply
	glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess);
	glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION, emission);
}

void GLWidget::drawScene(GLdouble dx /* = 0 */, GLdouble dy /* = 0 */, GLdouble dz /* = 0 */, GLdouble angle /* = 0 */)
{
	glTranslated(dx,dy,dz);
	glRotated(angle, 0.0, 0.0, 1.0);
	glRotated(-45.0, 1.0, 0.0, 0.0);
	glRotated(-45.0, 0.0, 0.0, 1.0);
	
	if(m_smoothshading)
		glShadeModel(GL_SMOOTH);
	else
		glShadeModel(GL_FLAT);

	if(m_setcullface)
		glCullFace(GL_FRONT);
	else
		glCullFace(GL_BACK);
	
	if(m_cullface)
		glEnable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
	
	gl_adjust_view(m_pThreeDimensionalShape);
	if (m_display_zoom_ma)
		gl_adjust_ma(m_pThreeDimensionalShape);
	
	setMaterial(COPPER);
	
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	
	if(m_antialiasing)
	{
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glLineWidth(1.0f);
	}
	else
	{
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
		glLineWidth(1.0f);
	}

	glPushMatrix();
	glEnable(GL_LIGHTING);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(3.0f, 1.0f);
	glDisable(GL_COLOR_MATERIAL);
	gl_draw_faces(m_pThreeDimensionalShape);	
	gl_draw_slab_envelop(m_pThreeDimensionalShape);
	gl_draw_ma_envelopes();

	gl_draw_input_dt(m_pThreeDimensionalShape);
	gl_draw_input_voronoi(m_pThreeDimensionalShape);

	gl_draw_slab_face(m_pThreeDimensionalShape);
	gl_draw_ma_faces();
	//glShadeModel(GL_SMOOTH);
	gl_draw_ma_spheres();
	gl_draw_slab_sphere(m_pThreeDimensionalShape);
	/*
	if(m_smoothshading)
		glShadeModel(GL_SMOOTH);
	else
		glShadeModel(GL_FLAT);
		*/
	
	
	glDisable(GL_LIGHTING);

	glColor3f(0.0f, 0.0f, 0.0f);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	gl_draw_edges(m_pThreeDimensionalShape);
	gl_draw_slab_edge(m_pThreeDimensionalShape);



	gl_draw_ma_footpoints();
	gl_draw_ma_edges();
	gl_draw_ma_vertices_bplist();

	glDisable(GL_POLYGON_OFFSET_FILL);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_LIGHTING);
	glPopMatrix();

	glDisable(GL_LIGHTING);
	gl_draw_vertices(m_pThreeDimensionalShape);
		gl_draw_facet_normals(m_pThreeDimensionalShape);
	gl_draw_ma_vertices();
	gl_draw_select_poles(select_poles);


	gl_draw_fake_boundary(m_pThreeDimensionalShape);
	gl_draw_boundary(m_pThreeDimensionalShape);
	gl_draw_nonmanifold(m_pThreeDimensionalShape);
	gl_draw_saved_vertex(m_pThreeDimensionalShape);
	gl_draw_stability_ratio(m_pThreeDimensionalShape);


	//gl_draw_initial_edge(m_pThreeDimensionalShape);
	//glEnable(GL_LIGHTING);
	//gl_draw_initial_face(m_pThreeDimensionalShape);
	//glDisable(GL_LIGHTING);

	gl_draw_slab_vertex(m_pThreeDimensionalShape);
	//gl_draw_sphere_edge(m_pThreeDimensionalShape);

	// 画出所选择的ma点
	gl_draw_selected_vertices();
	gl_draw_mean_square_error(m_pThreeDimensionalShape);

	//glEnable(GL_LIGHTING);
	//gl_draw_sphere_face(m_pThreeDimensionalShape);
	//glDisable(GL_LIGHTING);

	gl_draw_headbar();
	gl_draw_colorbar();

}

void GLWidget::gl_draw_mean_square_error(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if(m_display_mean_squre_error == false)
		return;
	if(m_pThreeDimensionalShape == NULL)
		return;
	double per;

	//if (m_display_slab_edge)
	//{
	SlabMesh *smesh = &(pThreeDimensionalShape->slab_mesh);
	glLineWidth(1.0);
	glEnable(GL_COLOR_MATERIAL);

	// find the max hyperbolic distance
	double bpmax = DBL_MIN;
	double bpmin = DBL_MAX;
	for (unsigned i = 0; i < smesh->edges.size(); i++)
	{
		if (!smesh->edges[i].first)	continue;
		bpmax = max(bpmax, smesh->edges[i].second->hyperbolic_weight);
		bpmin = min(bpmin, smesh->edges[i].second->hyperbolic_weight);
	}

	m_Colorbar_MinValue = bpmin;
	m_Colorbar_MaxValue = bpmax;

	for(unsigned i = 0; i < smesh->edges.size(); i ++)
	{
		if(!smesh->edges[i].first)
			continue;

		float mc[3];
		m_ColorRamp.RedGreenBlue((smesh->edges[i].second->hyperbolic_weight)*255., mc);
		glColor3fv(mc);

		Vector3d pos0, pos1;
		pos0 = smesh->vertices[smesh->edges[i].second->vertices_.first].second->sphere.center;
		pos1 = smesh->vertices[smesh->edges[i].second->vertices_.second].second->sphere.center;

		glBegin(GL_LINES);	
		glVertex3dv(pos0);
		glVertex3dv(pos1);
		glEnd();
	}
	return;
}

void GLWidget::gl_draw_faces(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if(m_display_faces == false)
		return;

	if(pThreeDimensionalShape == NULL)
		return;
	Mesh * pmesh = &(pThreeDimensionalShape->input);

	if(pmesh == NULL)
		return;

	glEnable(GL_COLOR_MATERIAL);
	//glDisable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	//glColor3f(.3f, .5f, .8f);
	glColor3f(136.0 / 255.0, 134.0 / 255.0, 119.0 / 255.0);
	Facet_iterator pFacet = pmesh->facets_begin();
	for(; pFacet != pmesh->facets_end(); pFacet ++)
	{
		glBegin(GL_POLYGON);
		gl_draw_facet(pFacet,1.);
		glEnd();
	}
	glFlush();

	glDisable(GL_COLOR_MATERIAL);

}

void GLWidget::gl_draw_edges(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if(m_display_edges == false)
		return;

	if(pThreeDimensionalShape == NULL)
		return;
	Mesh * pmesh = &(pThreeDimensionalShape->input);

	if(pmesh == NULL)
		return;
	
	
	glColor3f(0.5f, 0.5f, 0.5f);
	
	glLineWidth(1.);
	glBegin(GL_LINES);
	for(Edge_iterator ei = pmesh->edges_begin(); ei != pmesh->edges_end(); ei ++)
	{
		glVertex3d(ei->prev()->vertex()->point()[0] / pmesh->bb_diagonal_length, ei->prev()->vertex()->point()[1] / pmesh->bb_diagonal_length, ei->prev()->vertex()->point()[2] / pmesh->bb_diagonal_length);
		glVertex3d(ei->vertex()->point()[0] / pmesh->bb_diagonal_length, ei->vertex()->point()[1] / pmesh->bb_diagonal_length, ei->vertex()->point()[2] / pmesh->bb_diagonal_length);
	}
	glEnd();
	
	glFlush();
}
void GLWidget::gl_draw_select_poles(vector<vector<double>> select_poles) {
	if (m_display_select_pole == false)
		return;


	glEnable(GL_POINT_SMOOTH);
	glColor3f(.5f, .5f, .5f);
	glPointSize(10.0);  // 6~10 is ok. If it is too big, can not see the original point...

	int count = 0;
	for (int i = 0; i < select_poles.size(); i++) {
		glLoadName(count);
		glBegin(GL_POINTS);

		glColor3f(1.0, 0.0, 0.0); //设置颜色为红 RGB（255,0,0）
		glVertex3d(select_poles[i][0] , select_poles[i][1] , select_poles[i][2] );
		glEnd();
	}
	glDisable(GL_POINT_SMOOTH);
}


void GLWidget::gl_draw_vertices(ThreeDimensionalShape * pThreeDimensionalShape)
{

	// This is for surface point visualization...
	if(m_display_vertices == false)
		return;

	if(pThreeDimensionalShape == NULL)
		return;
	Mesh * pmesh = &(pThreeDimensionalShape->input);


	if(pmesh == NULL)
		return;

	glEnable(GL_POINT_SMOOTH);
	glColor3f(.5f, .5f, .5f);
	glPointSize(3.0);

	int count = 0;
	//glBegin(GL_POINTS);
	for(Vertex_iterator pVertex = pmesh->vertices_begin(); pVertex != pmesh->vertices_end(); pVertex ++, count++)
	{
		//if (count == 2 || count == 1750)
			glLoadName(count);
			glBegin(GL_POINTS);
			glColor3f(1.0, 0.0, 0.0); //设置颜色为红 RGB（255,0,0）
			//glVertex3d(pVertex->point()[0], pVertex->point()[1], pVertex->point()[2]);
				glVertex3d(pVertex->point()[0] / pmesh->bb_diagonal_length, pVertex->point()[1] / pmesh->bb_diagonal_length, pVertex->point()[2] / pmesh->bb_diagonal_length);
			glEnd();
			
	}
	//glEnd();


	glDisable(GL_POINT_SMOOTH);
}

void GLWidget::gl_draw_facet_normals(ThreeDimensionalShape * pThreeDimensionalShape){
	if(m_display_facet_normals == false)
		return;

	if(pThreeDimensionalShape == NULL)
		return;

	//Mesh * pmesh = &(pThreeDimensionalShape->input);
	//if(pmesh == NULL)
	//	return;

	//glColor3f(.5f, .5f, .5f);
	//glLineWidth(1.);

	//glBegin(GL_LINES);
	//for (int i = 0; i < pmesh->pFaceList.size(); i ++)
	//{
	//	MPMesh::Point p;
	//	pmesh->compute_facet_center(pmesh->pFaceList[i], p);
	//	glVertex3d(p.x(), p.y(), p.z());
	//	glVertex3d(p.x() + pmesh->pFaceList[i]->normal.x()/5,
	//		p.y() + pmesh->pFaceList[i]->normal.y()/5,
	//		p.z() + pmesh->pFaceList[i]->normal.z()/5);
	//}
	//glEnd();

	// draw the normal added to the boundary vertex
	SlabMesh *smesh = &(pThreeDimensionalShape->slab_mesh);
	if (smesh == NULL)
		return;

	glColor3f(1.0f, .0f, .0f);
	glLineWidth(1.);

	//unsigned cou = 0;
	glBegin(GL_LINES);
	for (int i = 0; i < smesh->vertices.size(); i ++)
	{
		if (smesh->vertices[i].first && smesh->vertices[i].second->fake_boundary_vertex)
			//&& smesh->vertices[i].second->boundary_edge_vec.size() >= 2)
			//&&(i == 18593 || i == 18867))
		{
			Vector3d p = smesh->vertices[i].second->sphere.center;
			Vector3d nor = smesh->vertices[i].second->boundVec;
			glVertex3d(p.X(), p.Y(), p.Z());
			glVertex3d(p.X() + nor.X() / 50,
				p.Y() + nor.Y() / 50,
				p.Z() + nor.Z() / 50);
		}
		//if (++cou >= m_display_boundary)	break;
	}
	glEnd();
}

void GLWidget::gl_draw_fake_boundary(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if (m_display_fake_boundary == false)
		return;
	if (m_pThreeDimensionalShape == NULL)
		return;
	SlabMesh *smesh = &(pThreeDimensionalShape->slab_mesh);
	if(smesh == NULL)
		return;

	glEnable(GL_POINT_SMOOTH);
	glColor3f(.0f, .0f, 0.0f);
	glPointSize(5.0);
	glBegin(GL_POINTS);

	for (unsigned i = 0; i < smesh->vertices.size(); i ++)
	{
		if (smesh->vertices[i].first 
			//&& smesh->vertices[i].second->boundary_vertex)
			&& smesh->vertices[i].second->boundary_edge_vec.size() >= 2)
		{
			glVertex3d(smesh->vertices[i].second->sphere.center.X(), smesh->vertices[i].second->sphere.center.Y(), smesh->vertices[i].second->sphere.center.Z());
		}
	}
	glEnd();
	glDisable(GL_POINT_SMOOTH);

	//glEnable(GL_LINE_SMOOTH);
	//QFont serifFont("Arial", 12, QFont::Normal);
	//for (unsigned i = 0; i < smesh->edges.size(); i ++)
	//{
	//	if(!smesh->edges[i].first)
	//		continue;

	//	Vector3d pos0, pos1;
	//	//if (smesh->edges[i].second->fake_boundary_edge)
	//	if (smesh->edges[i].second->faces_.size() == 1)
	//	{
	//		glColor3f(.0f, .0f, 1.0f);
	//		glLineWidth(2.0);
	//		pos0 = smesh->vertices[smesh->edges[i].second->vertices_.first].second->sphere.center;
	//		pos1 = smesh->vertices[smesh->edges[i].second->vertices_.second].second->sphere.center;
	//		glBegin(GL_LINES);	
	//		glVertex3dv(pos0);
	//		glVertex3dv(pos1);
	//		glEnd();
	//	}
	//}

	//unsigned cou = 0;
	//glEnable(GL_LINE_SMOOTH);
	//QFont serifFont("Arial", 12, QFont::Normal);
	//for (unsigned i = 0; i < smesh->vertices.size(); i ++)
	//{
	//	if (smesh->vertices[i].first && smesh->vertices[i].second->fake_boundary_vertex
	//		&& smesh->vertices[i].second->boundary_edge_vec.size() >= 2)
	//	{
	//		Vector3d pos0, pos1;
	//		for(auto si = smesh->vertices[i].second->boundary_edge_vec.begin(); si < smesh->vertices[i].second->boundary_edge_vec.end(); si++)
	//		{				
	//			if(smesh->edges[*si].first && smesh->edges[*si].second->fake_boundary_edge)
	//			{
	//				pos0 = smesh->vertices[smesh->edges[*si].second->vertices_.first].second->sphere.center;
	//				pos1 = smesh->vertices[smesh->edges[*si].second->vertices_.second].second->sphere.center;
	//				glEnable(GL_POINT_SMOOTH);
	//				glColor3f(.0f, .0f, 0.0f);
	//				glPointSize(5.0);
	//				glBegin(GL_POINTS);
	//				glVertex3dv(pos0);
	//				glVertex3dv(pos1);
	//				glEnd();
	//				glDisable(GL_POINT_SMOOTH);

	//				glColor3f(.0f, .0f, 1.0f);
	//				glLineWidth(2.0);
	//				glBegin(GL_LINES);	
	//				glVertex3dv(pos0);
	//				glVertex3dv(pos1);
	//				glEnd();
	//			}
	//		}
	//		if (++cou >= m_display_boundary)	
	//			break;
	//	}
	//}
}

void GLWidget::gl_draw_boundary(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if (m_display_boundary == 0)
		return;
	if (m_pThreeDimensionalShape == NULL)
		return;
	SlabMesh *smesh = &(pThreeDimensionalShape->slab_mesh);
	if(smesh == NULL)
		return;

}

void GLWidget::gl_draw_nonmanifold(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if (m_display_nonmanifold == false)
		return;
	if (m_pThreeDimensionalShape == NULL)
		return;
	SlabMesh *smesh = &(pThreeDimensionalShape->slab_mesh);
	if(smesh == NULL)
		return;

	//glEnable(GL_LIGHTING);
	//glEnable(GL_COLOR_MATERIAL);
	//glColor3f(.3f, .5f, .8f);

	//unsigned cou = 0;
	//for (unsigned i = 0; i < smesh->vertices.size(); i ++)
	//{
	//	if (smesh->vertices[i].first && smesh->vertices[i].second->fake_boundary_vertex
	//		&& smesh->vertices[i].second->boundary_edge_vec.size() >= 2)
	//	{
	//		glEnable(GL_LINE_SMOOTH);
	//		glColor3f(.3f, .5f, .8f);



	//		Cone * pcone;
	//		for(auto si = smesh->vertices[i].second->boundary_edge_vec.begin(); si < smesh->vertices[i].second->boundary_edge_vec.end(); si++)
	//		{
	//			if(smesh->edges[*si].first && smesh->edges[*si].second->fake_boundary_edge)
	//			{
	//				//// 启动混合并设置混合因子     
	//				//glEnable(GL_BLEND);      
	//				//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//				//glDepthMask(GL_FALSE);

	//				glColor3f(.3f, .5f, .8f);
	//				Vector3d pos = smesh->vertices[smesh->edges[*si].second->vertices_.first].second->sphere.center;
	//				double radius = smesh->vertices[smesh->edges[*si].second->vertices_.first].second->sphere.radius;
	//				glTranslated(pos[0],pos[1],pos[2]);
	//				gluSphere(spherequadric,radius,50,50);
	//				glTranslated(-pos[0],-pos[1],-pos[2]);

	//				glColor3f(.8f, .5f, .1f);
	//				pos = smesh->vertices[smesh->edges[*si].second->vertices_.second].second->sphere.center;
	//				radius = smesh->vertices[smesh->edges[*si].second->vertices_.second].second->sphere.radius;
	//				glTranslated(pos[0],pos[1],pos[2]);
	//				gluSphere(spherequadric,radius,50,50);
	//				glTranslated(-pos[0],-pos[1],-pos[2]);

	//				//glDepthMask(GL_TRUE);
	//				//pcone = &(smesh->edges[*si].second->cone);
	//				//glPushMatrix();
	//				//glTranslated(pcone->apex[0],pcone->apex[1],pcone->apex[2]);
	//				//glRotated(-pcone->rot_angle,pcone->rot_axis[0],pcone->rot_axis[1],pcone->rot_axis[2]);
	//				//gluCylinder(conequadric,pcone->base,pcone->top,pcone->height,120,1);
	//				//glPopMatrix();
	//			}
	//		}
	//		if (++cou >= m_display_boundary)	break;
	//	}
	//}

	//glDisable(GL_LIGHTING);
	//glDisable(GL_COLOR_MATERIAL);

	//glEnable(GL_POINT_SMOOTH);
	//glColor3f(1.0f, 0.0f, 0.0f);
	//glPointSize(5.0);
	//glBegin(GL_POINTS);
	//for (unsigned i = 0; i < smesh->vertices.size(); i ++)
	//{
	//	//if (smesh->vertices[i].first && smesh->vertices[i].second->boundary_vertex)
	//	//&& smesh->vertices[i].second->boundary_edge_vec.size() >= 2)
	//	if (smesh->vertices[i].first && smesh->vertices[i].second->fake_boundary_vertex) 
	//		//&& !smesh->vertices[i].second->boundary_vertex)
	//		//&& smesh->vertices[i].second->boundary_edge_vec.size() > 2)
	//	{
	//		glVertex3d(smesh->vertices[i].second->sphere.center.X(), smesh->vertices[i].second->sphere.center.Y(), smesh->vertices[i].second->sphere.center.Z());
	//	}
	//}
	//glEnd();
	//glDisable(GL_POINT_SMOOTH);

	glEnable(GL_LINE_SMOOTH);
	QFont serifFont("Arial", 12, QFont::Normal);
	for (unsigned i = 0; i < smesh->edges.size(); i ++)
	{
		if(!smesh->edges[i].first)
			continue;

		Vector3d pos0, pos1;
		if (smesh->edges[i].second->non_manifold_edge)
		{
			glColor3f(1.0f, 0.0f, 0.0f);
			glLineWidth(2.0);
			pos0 = smesh->vertices[smesh->edges[i].second->vertices_.first].second->sphere.center;
			pos1 = smesh->vertices[smesh->edges[i].second->vertices_.second].second->sphere.center;
			glBegin(GL_LINES);	
			glVertex3dv(pos0);
			glVertex3dv(pos1);
			glEnd();
		}
	}
}

void GLWidget::gl_draw_saved_vertex(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if (m_display_saved_vertex == false)
		return;
	if (m_pThreeDimensionalShape == NULL)
		return;

	if (m_display_slab_edge)
	{
		SlabMesh *smesh = &(pThreeDimensionalShape->slab_mesh);

		if(smesh == NULL)
			return;
		glEnable(GL_POINT_SMOOTH);
		glColor3f(.0f, .0f, 0.0f);
		glPointSize(5.0);
		glBegin(GL_POINTS);

		for (unsigned i = 0; i < smesh->vertices.size(); i ++)
		{
			if (smesh->vertices[i].first && smesh->vertices[i].second->saved_vertex)
			{
				glVertex3d(smesh->vertices[i].second->sphere.center.X(), smesh->vertices[i].second->sphere.center.Y(), smesh->vertices[i].second->sphere.center.Z());
			}
		}
		glEnd();
		glDisable(GL_POINT_SMOOTH);
	}else
	{

	}

	

	//for (unsigned i = 0; i < smesh->edges.size(); i++)
	//{
	//	if (smesh->edges[i].first && smesh->edges[i].second->faces_.size() == 0)
	//	{
	//		Vector3d ver;
	//		if (smesh->vertices[smesh->edges[i].second->vertices_.first].second->edges_.size() == 1)
	//			ver = smesh->vertices[smesh->edges[i].second->vertices_.first].second->sphere.center;
	//			
	//		if (smesh->vertices[smesh->edges[i].second->vertices_.second].second->edges_.size() == 1)
	//			ver = smesh->vertices[smesh->edges[i].second->vertices_.second].second->sphere.center;
	//		
	//		glVertex3d(ver.X(), ver.Y(), ver.Z());
	//	}
	//}


}

void GLWidget::gl_draw_stability_ratio(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if(m_display_stability_ratio == false)
		return;
	if(m_pThreeDimensionalShape == NULL)
		return;
	double per;

	//if (m_display_slab_edge)
	//{
		SlabMesh *smesh = &(pThreeDimensionalShape->slab_mesh);
		glLineWidth(1.0);
		glEnable(GL_COLOR_MATERIAL);

		// find the max hyperbolic distance
		double bpmax = DBL_MIN;
		double bpmin = DBL_MAX;
		for (unsigned i = 0; i < smesh->edges.size(); i++)
		{
			if (!smesh->edges[i].first)	continue;
			bpmax = max(bpmax, smesh->edges[i].second->hyperbolic_weight);
			bpmin = min(bpmin, smesh->edges[i].second->hyperbolic_weight);
		}

		m_Colorbar_MinValue = bpmin;
		m_Colorbar_MaxValue = bpmax;

		for(unsigned i = 0; i < smesh->edges.size(); i ++)
		{
			//if(!smesh->faces[i].first || smesh->faces[i].second->hyperbolic_weight > 0.0)
			if(!smesh->edges[i].first)
				continue;
			//glColor3f(1.0f, 0.0f, 0.0f); 
			 
			//per = smesh->edges[i].second->hyperbolic_weight / bpmax;
			//per = smesh->edges[i].second->hyperbolic_weight / bpmax;
			//Vector3d rgb;
			//if (per <= 0.1)
			//	rgb = HSVToRGB(Vector3d(240 * (1 - per / 0.3) + 90 * (per / 0.3), 1.0, 1.0));
			//else
			//	rgb = HSVToRGB(Vector3d(90 - 90 * per, 1.0, 1.0));
			//glColor3f(rgb[0], rgb[1], rgb[2]);

			float mc[3];
			m_ColorRamp.RedGreenBlue((smesh->edges[i].second->hyperbolic_weight)*255., mc);
			glColor3fv(mc);

			Vector3d pos0, pos1;
			pos0 = smesh->vertices[smesh->edges[i].second->vertices_.first].second->sphere.center;
			pos1 = smesh->vertices[smesh->edges[i].second->vertices_.second].second->sphere.center;

			glBegin(GL_LINES);	
			glVertex3dv(pos0);
			glVertex3dv(pos1);
			glEnd();

		}

		//// find the max hyperbolic distance
		//double bpmax = DBL_MIN;
		//double bpmin = DBL_MAX;
		//for (unsigned i = 0; i < smesh->faces.size(); i++)
		//{
		//	if (!smesh->faces[i].first)	continue;
		//	bpmax = max(bpmax, smesh->faces[i].second->hyperbolic_weight);
		//	bpmin = min(bpmin, smesh->faces[i].second->hyperbolic_weight);
		//}

		//for(unsigned i = 0; i < smesh->faces.size(); i ++)
		//{
		//	//if(!smesh->faces[i].first || smesh->faces[i].second->hyperbolic_weight > 0.0)
		//	if(!smesh->faces[i].first)
		//		continue;
		//	//glColor3f(1.0f, 0.0f, 0.0f); 
		//	 
		//	per = smesh->faces[i].second->hyperbolic_weight / bpmax;
		//	Vector3d rgb;
		//	if (per <= 0.1)
		//		rgb = HSVToRGB(Vector3d(240 * (1 - per / 0.3) + 90 * (per / 0.3), 1.0, 1.0));
		//	else
		//		rgb = HSVToRGB(Vector3d(90 - 90 * per, 1.0, 1.0));
		//	glColor3f(rgb[0], rgb[1], rgb[2]);

		//	Wm4::Vector3d normal = smesh->faces[i].second->normal;
		//	glBegin(GL_POLYGON);
		//	glNormal3dv(normal);
		//	for(std::set<unsigned>::iterator si = smesh->faces[i].second->vertices_.begin(); si != smesh->faces[i].second->vertices_.end(); si ++)
		//	{
		//		//per = smesh->vertices[*si].second->mean_square_error / smesh->max_mean_squre_error;
		//		////per = smesh->vertices[*si].second->bplist.size() * 1.0 / bpmax ;
		//		//Vector3d rgb;
		//		//if (per <= 0.1)
		//		//{
		//		//	rgb = HSVToRGB(Vector3d(240 * (1 - per / 0.3) + 90 * (per / 0.3), 1.0, 1.0));
		//		//} 
		//		//else
		//		//{
		//		//	rgb = HSVToRGB(Vector3d(90 - 90 * per, 1.0, 1.0));
		//		//}

		//		//glColor3f(rgb[0], rgb[1], rgb[2]);
		//		glVertex3dv((smesh->vertices[*si].second->sphere.center));
		//	}
		//	glEnd();

		//	glBegin(GL_POLYGON);
		//	glNormal3dv(-normal);
		//	for(std::set<unsigned>::iterator si = smesh->faces[i].second->vertices_.end(); si != smesh->faces[i].second->vertices_.begin();)
		//	{
		//		//per = smesh->vertices[*--si].second->mean_square_error / smesh->max_mean_squre_error;
		//		////per = smesh->vertices[*--si].second->bplist.size() * 1.0 / bpmax;
		//		//Vector3d rgb;
		//		//if (per <= 0.1)
		//		//{
		//		//	rgb = HSVToRGB(Vector3d(240 * (1 - per / 0.3) + 90 * (per / 0.3), 1.0, 1.0));
		//		//} 
		//		//else
		//		//{
		//		//	rgb = HSVToRGB(Vector3d(90 - 90 * per, 1.0, 1.0));
		//		//}

		//		//glColor3f(rgb[0], rgb[1], rgb[2]);
		//		glVertex3dv((smesh->vertices[*--si].second->sphere.center));
		//	}
		//	glEnd();
		//}
	//}
	//else
	//{
	//	VQEMMesh *smesh = &(pThreeDimensionalShape->ma_qem_mesh);
	//	glLineWidth(0.5);
	//	glEnable(GL_COLOR_MATERIAL);

	//	for(unsigned i = 0; i < smesh->faces.size(); i ++)
	//	{
	//		if(!smesh->faces[i].first)
	//			continue;

	//		Wm4::Vector3d normal = smesh->faces[i].second->normal;

	//		glBegin(GL_POLYGON);
	//		glNormal3dv(normal);
	//		for(std::set<unsigned>::iterator si = smesh->faces[i].second->vertices_.begin(); si != smesh->faces[i].second->vertices_.end(); si ++)
	//		{
	//			per = smesh->vertices[*si].second->mean_square_error / smesh->max_mean_squre_error;
	//			Vector3d rgb;
	//			if (per <= 0.1)
	//			{
	//				rgb = HSVToRGB(Vector3d(240 * (1 - per / 0.3) + 90 * (per / 0.3), 1.0, 1.0));
	//			} 
	//			else
	//			{
	//				rgb = HSVToRGB(Vector3d(90 - 90 * per, 1.0, 1.0));
	//			}

	//			glColor3f(rgb[0], rgb[1], rgb[2]);
	//			glVertex3dv((smesh->vertices[*si].second->sphere.center));
	//		}
	//		glEnd();

	//		glBegin(GL_POLYGON);
	//		glNormal3dv(-normal);
	//		for(std::set<unsigned>::iterator si = smesh->faces[i].second->vertices_.end(); si != smesh->faces[i].second->vertices_.begin();)
	//		{
	//			per = smesh->vertices[*--si].second->mean_square_error / smesh->max_mean_squre_error;
	//			Vector3d rgb;
	//			if (per <= 0.1)
	//			{
	//				rgb = HSVToRGB(Vector3d(240 * (1 - per / 0.3) + 90 * (per / 0.3), 1.0, 1.0));
	//			} 
	//			else
	//			{
	//				rgb = HSVToRGB(Vector3d(90 - 90 * per, 1.0, 1.0));
	//			}

	//			glColor3f(rgb[0], rgb[1], rgb[2]);
	//			glVertex3dv((smesh->vertices[*si].second->sphere.center));
	//		}
	//		glEnd();
	//	}
	//}

	return;
}

void GLWidget::gl_draw_slab_vertex(ThreeDimensionalShape * pThreeDimensionalShape){
	if (m_display_slab_vertex == false)
		return;
	if (m_pThreeDimensionalShape == NULL)
		return;
	SlabMesh *smesh = &(pThreeDimensionalShape->slab_mesh);
	if(smesh == NULL)
		return;
	glEnable(GL_POINT_SMOOTH);
	glPointSize(5.0);
	for (unsigned i = 0; i < smesh->vertices.size(); i ++)
	{
		//if (smesh->vertices[i].first && smesh->vertices[i].second->bplist.size() != 0)
		//{
		//	glColor3f(.0f, 1.0f, .0f);
		//	glVertex3d(smesh->vertices[i].second->sphere.center.X(), smesh->vertices[i].second->sphere.center.Y(), smesh->vertices[i].second->sphere.center.Z());
		//}
		//else 
		//if (smesh->vertices[i].first && (i == 1750 || i == 2118))
		if (smesh->vertices[i].first)
				//&& smesh->vertices[i].second->bplist.size() == 0)
		{
			glColor3f(0.0f, 0.0f, 0.0f);
			glLoadName(i);
			glBegin(GL_POINTS);
			glVertex3d(smesh->vertices[i].second->sphere.center.X(), smesh->vertices[i].second->sphere.center.Y(), smesh->vertices[i].second->sphere.center.Z());
			glEnd();
		}
	}

	glDisable(GL_POINT_SMOOTH);
}

void GLWidget::gl_draw_selected_vertices()
{
	if (hits > 0)
	{
		int i = selectBuffer[3];
		glEnable(GL_POINT_SMOOTH);
		glPointSize(8.0);
		glColor3f(1.0f, 0.0f, 0.0f);
		glBegin(GL_POINTS);
		glVertex3d(m_pThreeDimensionalShape->slab_mesh.vertices[i].second->sphere.center.X(), 
			m_pThreeDimensionalShape->slab_mesh.vertices[i].second->sphere.center.Y(), 
			m_pThreeDimensionalShape->slab_mesh.vertices[i].second->sphere.center.Z());
		glEnd();
		glDisable(GL_POINT_SMOOTH);


		glEnable(GL_LIGHTING);
		glEnable(GL_COLOR_MATERIAL);
		glColor3f(.3f, .5f, .8f);
		Vector3d pos = m_pThreeDimensionalShape->slab_mesh.vertices[i].second->sphere.center;
		double radius = m_pThreeDimensionalShape->slab_mesh.vertices[i].second->sphere.radius;
		glTranslated(pos[0],pos[1],pos[2]);
		gluSphere(spherequadric,radius,40,40);
		glTranslated(-pos[0],-pos[1],-pos[2]);

		double min_coll = DBL_MAX;
		int index = 0;
		set<unsigned> eset = m_pThreeDimensionalShape->slab_mesh.vertices[i].second->edges_;
		for (set<unsigned>::iterator si = eset.begin(); si != eset.end(); si++)
		{
			if (m_pThreeDimensionalShape->slab_mesh.edges[*si].second->collapse_cost < min_coll)
			{
				index = *si;
				min_coll = m_pThreeDimensionalShape->slab_mesh.edges[*si].second->collapse_cost;
			}
		}

		int v1 = m_pThreeDimensionalShape->slab_mesh.edges[index].second->vertices_.first;
		int v2 = m_pThreeDimensionalShape->slab_mesh.edges[index].second->vertices_.second;
		int ind = v1 == i ? v2 : v1;
		pos = m_pThreeDimensionalShape->slab_mesh.vertices[ind].second->sphere.center;
		radius = m_pThreeDimensionalShape->slab_mesh.vertices[ind].second->sphere.radius;
		glTranslated(pos[0],pos[1],pos[2]);
		gluSphere(spherequadric,radius,40,40);
		glTranslated(-pos[0],-pos[1],-pos[2]);
		glDisable(GL_LIGHTING);
		glDisable(GL_COLOR_MATERIAL);


	}
}

// 调medial axis边颜色
void GLWidget::gl_draw_slab_edge(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if(!m_display_slab_edge)
		return;
	if (m_pThreeDimensionalShape == NULL)
		return;
	SlabMesh *smesh = &(pThreeDimensionalShape->slab_mesh);
	if(smesh == NULL)
		return;

	glEnable(GL_LINE_SMOOTH);
	QFont serifFont("Arial", 12, QFont::Normal);
	for (unsigned i = 0; i < smesh->edges.size(); i ++)
	{
		if(!smesh->edges[i].first)
			continue;

		Vector3d pos0, pos1;

		// 边的颜色
		glColor3f(0.5f, 0.5f, 0.5f);
		glLineWidth(0.8);

		pos0 = smesh->vertices[smesh->edges[i].second->vertices_.first].second->sphere.center;
		pos1 = smesh->vertices[smesh->edges[i].second->vertices_.second].second->sphere.center;

		glBegin(GL_LINES);	
		glVertex3dv(pos0);
		glVertex3dv(pos1);
		glEnd();
	}
	return;
}

// 调medial axis面颜色
void GLWidget::gl_draw_slab_face(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if(m_display_slab_faces == false)
		return;
	if(m_pThreeDimensionalShape == NULL)
		return;

	glLineWidth(0.5);
	glColor3f(0.5f, 0.5f, 0.5f);

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_COLOR_MATERIAL);
	//glDisable(GL_LIGHTING);
	//glColor3f(.3f, .5f, .8f);

	for(unsigned i = 0; i < m_pThreeDimensionalShape->slab_mesh.faces.size(); i ++)
	{
		if(!m_pThreeDimensionalShape->slab_mesh.faces[i].first)
			continue;
		Wm4::Vector3d normal = m_pThreeDimensionalShape->slab_mesh.faces[i].second->normal;

		//if((m_pThreeDimensionalShape->input_nmm.faces[i].second->vertices_.find(5509) == m_pThreeDimensionalShape->input_nmm.faces[i].second->vertices_.end())
		//	&& (m_pThreeDimensionalShape->input_nmm.faces[i].second->vertices_.find(5522) == m_pThreeDimensionalShape->input_nmm.faces[i].second->vertices_.end()))
		//	continue;

		float mc[3];
		glColor3fv(mc);

		//glColor3f(82.0 / 255, 134.0 / 255.0, 139.0 / 255);
		//glColor3f(0.0f, 100.0 / 255.0, 0.0f);

		// DeepSkyBlue1
		//glColor3f(0.0f, 191.0 / 255.0, 1.0f);

		//grey51
		//glColor3f(130.0 / 255.0, 130.0 / 255.0, 130.0 / 255.0);

		//Tomato1
		//glColor3f(255.0 / 255.0, 99.0 / 255.0, 71.0 / 255.0);

		//Firebrick1
		//glColor3f(255.0 / 255.0, 48.0 / 255.0, 48.0 / 255.0);

		//MediumPurple2
		//glColor3f(159 / 255.0, 121 / 255.0, 238 / 255.0);

		glColor3f(red / 255.0, green / 255.0, blue / 255.0);


		glBegin(GL_POLYGON);
		glNormal3dv(normal);
		for(std::set<unsigned>::iterator si = m_pThreeDimensionalShape->slab_mesh.faces[i].second->vertices_.begin(); si != m_pThreeDimensionalShape->slab_mesh.faces[i].second->vertices_.end(); si ++)
			glVertex3dv(m_pThreeDimensionalShape->slab_mesh.vertices[*si].second->sphere.center);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3dv(-normal);
		for(std::set<unsigned>::iterator si = m_pThreeDimensionalShape->slab_mesh.faces[i].second->vertices_.end(); si != m_pThreeDimensionalShape->slab_mesh.faces[i].second->vertices_.begin();)
			glVertex3dv(m_pThreeDimensionalShape->slab_mesh.vertices[*(--si)].second->sphere.center);
		glEnd();

	}
	glEnable(GL_LIGHTING);
	return;
}


void GLWidget::gl_draw_slab_envelop(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if(m_display_slab_envelop == false)
		return;

	if(m_pThreeDimensionalShape == NULL)
		return;
	SlabMesh *smesh = &(pThreeDimensionalShape->slab_mesh);

	// here
	cout << "print<< I am here" << endl;
	if(smesh == NULL)
		return;

	glLineWidth(0.5);

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_COLOR_MATERIAL);

	// 白金色
	glColor3f(136.0 / 255.0, 134.0 / 255.0, 119.0 / 255.0);
	
	Cone * pcone;

	for(unsigned i = 0; i < smesh->edges.size(); i ++)
	{
		if(smesh->edges[i].first)
		{
			pcone = &(smesh->edges[i].second->cone);
			glPushMatrix();
			glTranslated(pcone->apex[0],pcone->apex[1],pcone->apex[2]);
			glRotated(-pcone->rot_angle,pcone->rot_axis[0],pcone->rot_axis[1],pcone->rot_axis[2]);
			gluCylinder(conequadric,pcone->base,pcone->top,pcone->height,120,1);
			glPopMatrix();
		}
	}

	for(unsigned i = 0; i < smesh->faces.size(); i ++)
	{
		if(!smesh->faces[i].first)
			continue;
		if(!smesh->faces[i].second->valid_st)
			continue;
		for(unsigned k = 0; k < 2; k ++)
		{
			SimpleTriangle * st;
			st = &(smesh->faces[i].second->st[k]);

			glBegin(GL_POLYGON);
			glNormal3dv(st->normal);
			glVertex3dv(st->v[0]);
			glVertex3dv(st->v[1]);
			glVertex3dv(st->v[2]);
			glEnd();
			glBegin(GL_POLYGON);
			glNormal3dv(st->normal);
			glVertex3dv(st->v[1]);
			glVertex3dv(st->v[0]);
			glVertex3dv(st->v[2]);
			glEnd();
		}
	}
	glDisable(GL_AUTO_NORMAL);

	glDisable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);
	return;
}

void GLWidget::gl_draw_slab_sphere(ThreeDimensionalShape * pThreeDimensionalShape){
	if (m_display_slab_sphere == false)
		return;
	if(pThreeDimensionalShape == NULL)
		return;
	SlabMesh *smesh = &(pThreeDimensionalShape->slab_mesh);
	if(smesh == NULL)
		return;
	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);

	// 白金色
	glColor3f(136.0 / 255.0, 134.0 / 255.0, 119.0 / 255.0);

	for(unsigned i = 0; i < pThreeDimensionalShape->slab_mesh.vertices.size(); i ++)
	{
		if(pThreeDimensionalShape->slab_mesh.vertices[i].first)
		{
			Vector3d pos = pThreeDimensionalShape->slab_mesh.vertices[i].second->sphere.center;
			double radius = pThreeDimensionalShape->slab_mesh.vertices[i].second->sphere.radius;
			glTranslated(pos[0],pos[1],pos[2]);
			gluSphere(spherequadric,radius,40,40);
			glTranslated(-pos[0],-pos[1],-pos[2]);
		}
	}
	glDisable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);
	return;
}

void GLWidget::process_select(int xPos, int yPos)
{
}

void GLWidget::gl_draw_input_dt(ThreeDimensionalShape * pThreeDimensionalShape)
{	
	if(m_display_input_dt == false)
		return;
	if(pThreeDimensionalShape == NULL)
		return;

	glEnable(GL_COLOR_MATERIAL);

	glLineWidth(2.);
	glColor3f(0.2f, 0.2f, 0.2f);

	glEnable(GL_POINT_SMOOTH);

	glPointSize(5.);
	glColor3f(1., 0., 0.);

	for(Finite_cells_iterator_t fci = pThreeDimensionalShape->input.dt.finite_cells_begin(); 
		fci != pThreeDimensionalShape->input.dt.finite_cells_end();
		fci ++)
	{
		Point_t cent = CGAL::circumcenter(pThreeDimensionalShape->input.dt.tetrahedron(fci));
		Wm4::Vector3d centroidp(cent[0], cent[1], cent[2]);

		if(fci->info().inside == false)
			continue;

		glBegin(GL_LINES);
			glVertex3d(fci->vertex(0)->point()[0],fci->vertex(0)->point()[1],fci->vertex(0)->point()[2]);
			glVertex3d(fci->vertex(1)->point()[0],fci->vertex(1)->point()[1],fci->vertex(1)->point()[2]);
			glVertex3d(fci->vertex(0)->point()[0],fci->vertex(0)->point()[1],fci->vertex(0)->point()[2]);
			glVertex3d(fci->vertex(2)->point()[0],fci->vertex(2)->point()[1],fci->vertex(2)->point()[2]);
			glVertex3d(fci->vertex(0)->point()[0],fci->vertex(0)->point()[1],fci->vertex(0)->point()[2]);
			glVertex3d(fci->vertex(3)->point()[0],fci->vertex(3)->point()[1],fci->vertex(3)->point()[2]);
			glVertex3d(fci->vertex(1)->point()[0],fci->vertex(1)->point()[1],fci->vertex(1)->point()[2]);
			glVertex3d(fci->vertex(2)->point()[0],fci->vertex(2)->point()[1],fci->vertex(2)->point()[2]);
			glVertex3d(fci->vertex(1)->point()[0],fci->vertex(1)->point()[1],fci->vertex(1)->point()[2]);
			glVertex3d(fci->vertex(3)->point()[0],fci->vertex(3)->point()[1],fci->vertex(3)->point()[2]);
			glVertex3d(fci->vertex(2)->point()[0],fci->vertex(2)->point()[1],fci->vertex(2)->point()[2]);
			glVertex3d(fci->vertex(3)->point()[0],fci->vertex(3)->point()[1],fci->vertex(3)->point()[2]);
		glEnd();
	}
}

void GLWidget::gl_draw_input_voronoi(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if(m_display_input_voronoi == false)
		return;
	if(pThreeDimensionalShape == NULL)
		return;

	glEnable(GL_COLOR_MATERIAL);

	glEnable(GL_POINT_SMOOTH);

	glPointSize(1.);
	glColor3f(1., 0., 0.);

	Triangulation * pt;
	pt = &(pThreeDimensionalShape->input.dt);

	glDisable(GL_LIGHTING);
	for(Finite_cells_iterator_t fci = pt->finite_cells_begin(); 
		fci != pt->finite_cells_end();
		fci ++)
	{
		Point_t cent = CGAL::circumcenter(pt->tetrahedron(fci));
		Wm4::Vector3d centroidp(cent[0], cent[1], cent[2]);
		if(fci->info().inside == false)
			continue;

		if(fci->info().is_pole)
		{
			glPointSize(10.);
			glColor3f(1., 0., 0.);
		}
		else
		{
			glPointSize(5.);
			glColor3f(0.,1.,0.);
		}
		glPointSize(1.);
		glColor3f(0.8,0.8,0.8);
		glBegin(GL_POINTS);
			glVertex3d(cent[0], cent[1], cent[2]);
		glEnd();
	}
	glEnable(GL_LIGHTING);



	glLineWidth(2.);
	glColor3f(0.5f, 0.5f, 0.5f);

	glEnable(GL_LINE_SMOOTH);

	glColor3f(0.5f, 0.5f, 0.5f);
	glLineWidth(1.);

	for(Finite_facets_iterator_t ffi = pt->finite_facets_begin(); 
		ffi != pt->finite_facets_end();
		ffi ++)
	{
		Triangulation::Object o = pt->dual(*ffi);
		if(const Triangulation::Segment *s = CGAL::object_cast<Triangulation::Segment>(&o))
		{
			if( (ffi->first->info().inside == false) || (pt->mirror_facet(*ffi).first->info().inside == false) )
				continue;
			glBegin(GL_LINES);
			glVertex3d(s->vertex(0)[0],s->vertex(0)[1],s->vertex(0)[2]);
			glVertex3d(s->vertex(1)[0],s->vertex(1)[1],s->vertex(1)[2]);
			glEnd();
		}
		else if(const Triangulation::Ray *r = CGAL::object_cast<Triangulation::Ray>(&o))
		{
			// does not handle infinite cells at this time
				continue;
			glBegin(GL_LINES);
			glVertex3d(s->point(0)[0],s->point(0)[1],s->point(0)[2]);
			glVertex3d(s->point(1)[0],s->point(1)[1],s->point(1)[2]);
			glEnd();
		}
	}

	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_AUTO_NORMAL);

	glColor3f(.3f, .5f, .8f);
	for(Finite_edges_iterator_t fei = pt->finite_edges_begin();
		fei != pt->finite_edges_end();
		fei ++)
	{
		// get all cells incident to the edge
		std::vector<Cell_handle_t> vec_ch;
		Cell_circulator_t cc = pt->incident_cells(*fei);
		Cell_circulator_t cc_begin = cc;
		do
		{
			vec_ch.push_back(cc);
			cc++;
		}while(cc!= cc_begin);
		
		bool all_finite_inside = true;
		for(unsigned int k = 0; k < vec_ch.size(); k ++)
		{
			if(pt->is_infinite(vec_ch[k]))
				all_finite_inside = false;
			else if(vec_ch[k]->info().inside == false)
				all_finite_inside = false;
		}
		if(all_finite_inside)
		{
			Wm4::Vector3d v0,v1,v2;
			v0 = to_wm4(CGAL::circumcenter(pt->tetrahedron(vec_ch[0])));
			v1 = to_wm4(CGAL::circumcenter(pt->tetrahedron(vec_ch[1])));
			v2 = to_wm4(CGAL::circumcenter(pt->tetrahedron(vec_ch[2])));
			Wm4::Vector3d fnorm;
			fnorm = (v1-v0).Cross(v2-v0);
			fnorm.Normalize();
			for(unsigned int k = 1; k < vec_ch.size() - 1; k ++)
			{
				glBegin(GL_POLYGON);
				glNormal3d(fnorm[0], fnorm[1],fnorm[2]);
				glVertex3d(
					CGAL::circumcenter(pt->tetrahedron(vec_ch[0]))[0],
					CGAL::circumcenter(pt->tetrahedron(vec_ch[0]))[1],
					CGAL::circumcenter(pt->tetrahedron(vec_ch[0]))[2]);
				glVertex3d(
					CGAL::circumcenter(pt->tetrahedron(vec_ch[k]))[0],
					CGAL::circumcenter(pt->tetrahedron(vec_ch[k]))[1],
					CGAL::circumcenter(pt->tetrahedron(vec_ch[k]))[2]);
				glVertex3d(
					CGAL::circumcenter(pt->tetrahedron(vec_ch[k+1]))[0],
					CGAL::circumcenter(pt->tetrahedron(vec_ch[k+1]))[1],
					CGAL::circumcenter(pt->tetrahedron(vec_ch[k+1]))[2]);
				glEnd();
				glBegin(GL_POLYGON);
				glNormal3d(-fnorm[0], -fnorm[1],-fnorm[2]);
				glVertex3d(
					CGAL::circumcenter(pt->tetrahedron(vec_ch[0]))[0],
					CGAL::circumcenter(pt->tetrahedron(vec_ch[0]))[1],
					CGAL::circumcenter(pt->tetrahedron(vec_ch[0]))[2]);
				glVertex3d(
					CGAL::circumcenter(pt->tetrahedron(vec_ch[k+1]))[0],
					CGAL::circumcenter(pt->tetrahedron(vec_ch[k+1]))[1],
					CGAL::circumcenter(pt->tetrahedron(vec_ch[k+1]))[2]);
				glVertex3d(
					CGAL::circumcenter(pt->tetrahedron(vec_ch[k]))[0],
					CGAL::circumcenter(pt->tetrahedron(vec_ch[k]))[1],
					CGAL::circumcenter(pt->tetrahedron(vec_ch[k]))[2]);
				glEnd();
			}
		}
	}
}

void GLWidget::gl_draw_ma_vertices()
{
	if(m_display_ma_vertices == false)
		return;
	if(m_pThreeDimensionalShape == NULL)
		return;

	glEnable(GL_COLOR_MATERIAL);

	glEnable(GL_POINT_SMOOTH);

	glPointSize(5.);
	glPointSize(10.);
	glPointSize(5.);
	glColor3f(1., 0., 0.);
	QFont serifFont("Arial", 12, QFont::Normal);

	glPointSize(10.);
	glColor3f(1., 0., 0.);

	for(unsigned i = 0; i < m_pThreeDimensionalShape->input_nmm.vertices.size(); i ++)
	{
		if(!m_pThreeDimensionalShape->input_nmm.vertices[i].first)
			continue;
		if(m_pThreeDimensionalShape->input_nmm.vertices[i].second->is_pole)
		{
			glPointSize(10.);
			glColor3f(1., 0., 0.);
		}
		else
		{
			glPointSize(5.);
			glColor3f(0.,1.,0.);
		}
		glPointSize(5.);
		glColor3f(0.,1.,0.);
		glBegin(GL_POINTS);	
		glVertex3dv(m_pThreeDimensionalShape->input_nmm.vertices[i].second->sphere.center);
		glEnd();

		glColor3f((GLfloat).3,(GLfloat).3,(GLfloat).3);
		QString str;
	}
	return;
	
}

void GLWidget::gl_draw_ma_vertices_bplist()
{
	if(m_display_ma_vertices_bplist == false)
		return;
	if(m_pThreeDimensionalShape == NULL)
		return;

	glEnable(GL_COLOR_MATERIAL);

	glEnable(GL_LINE_SMOOTH);

	glPointSize(2.);

	glColor3f(.1,.1,.1);

	for(unsigned i = 0; i < m_pThreeDimensionalShape->input_nmm.vertices.size(); i ++)
	{
		if(!m_pThreeDimensionalShape->input_nmm.vertices[i].first)
			continue;
		for(std::set<unsigned>::iterator si = m_pThreeDimensionalShape->input_nmm.vertices[i].second->bplist.begin();
			si != m_pThreeDimensionalShape->input_nmm.vertices[i].second->bplist.end(); si ++)
		{
			glBegin(GL_LINES);
				glVertex3dv(m_pThreeDimensionalShape->input_nmm.vertices[i].second->sphere.center);
				glVertex3dv(m_pThreeDimensionalShape->input_nmm.BoundaryPoints[*si]);
			glEnd();
		}
	}
	return;

}

void GLWidget::gl_draw_ma_edges()
{
	if(m_display_ma_edges == false)
		return;
	if(m_pThreeDimensionalShape == NULL)
		return;

	glLineWidth(2.);
	glLineWidth(5.);
	glLineWidth(1.);
	glColor3f(0.5f, 0.5f, 0.5f);

	glEnable(GL_LINE_SMOOTH);
	QFont serifFont("Arial", 12, QFont::Normal);
	for(unsigned i = 0; i < m_pThreeDimensionalShape->input_nmm.numEdges; i ++)
	{
		if(!m_pThreeDimensionalShape->input_nmm.edges[i].first)
			continue;
		if(m_pThreeDimensionalShape->input_nmm.edges[i].second->is_dangling)
		{
			glColor3f(.0f, 0.f, 0.8f);
			glLineWidth(2.);
		}
		else if(m_pThreeDimensionalShape->input_nmm.edges[i].second->is_non_manifold)
		{
			glColor3f(.8f, 0.f, 0.f);
			glLineWidth(4.);
		}
		else
		{
			glColor3f(0.5f, 0.5f, 0.5f);
			glLineWidth(1.);
		}
		glColor3f(0.5f, 0.5f, 0.5f);
		glLineWidth(1.);
		
		Vector3d pos0 = m_pThreeDimensionalShape->input_nmm.vertices[m_pThreeDimensionalShape->input_nmm.edges[i].second->vertices_.first].second->sphere.center;
		Vector3d pos1 = m_pThreeDimensionalShape->input_nmm.vertices[m_pThreeDimensionalShape->input_nmm.edges[i].second->vertices_.second].second->sphere.center;
		
		Vector3d mp;
		mp = 0.5 * (pos0 + pos1);

		glBegin(GL_LINES);	
			glVertex3dv(pos0);
			glVertex3dv(pos1);
		glEnd();

		glColor3f((GLfloat).3,(GLfloat).3,(GLfloat).3);

	}
	return;

}

void GLWidget::gl_draw_ma_faces()
{
	if(m_display_ma_faces == false)
		return;
	if(m_pThreeDimensionalShape == NULL)
		return;

	glLineWidth(0.5);
	glColor3f(0.5f, 0.5f, 0.5f);

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_COLOR_MATERIAL);

	//glColor3f(.3f, .5f, .8f);
	glColor3f(0.8f, .5f, .25f);

	for(unsigned i = 0; i < m_pThreeDimensionalShape->input_nmm.numFaces; i ++)
	{
		if(!m_pThreeDimensionalShape->input_nmm.faces[i].first)
			continue;
		Wm4::Vector3d normal = m_pThreeDimensionalShape->input_nmm.faces[i].second->normal;

		glBegin(GL_POLYGON);
		glNormal3dv(normal);
		for(std::set<unsigned>::iterator si = m_pThreeDimensionalShape->input_nmm.faces[i].second->vertices_.begin(); si != m_pThreeDimensionalShape->input_nmm.faces[i].second->vertices_.end(); si ++)
			glVertex3dv(m_offsetfactor*m_pThreeDimensionalShape->input_nmm.vertices[*si].second->sphere.center+(1.-m_offsetfactor)*m_pThreeDimensionalShape->input_nmm.faces[i].second->centroid);
		glEnd();

		glBegin(GL_POLYGON);
		glNormal3dv(-normal);
		for(std::set<unsigned>::iterator si = m_pThreeDimensionalShape->input_nmm.faces[i].second->vertices_.end(); si != m_pThreeDimensionalShape->input_nmm.faces[i].second->vertices_.begin();)
			glVertex3dv(m_offsetfactor*m_pThreeDimensionalShape->input_nmm.vertices[*(--si)].second->sphere.center+(1.-m_offsetfactor)*m_pThreeDimensionalShape->input_nmm.faces[i].second->centroid);
		glEnd();

	}
	return;

}

void GLWidget::gl_draw_ma_footpoints()
{
	if(m_display_ma_footpoints == false)
		return;
	if(m_pThreeDimensionalShape == NULL)
		return;


	glColor3f(.9f,.2f,.2f);
	glPointSize(2.0);
	glEnable(GL_POINT_SMOOTH);
	
	for(unsigned int i = 0; i < m_pThreeDimensionalShape->input_nmm.BoundaryPoints.size(); i ++)
	{
		if(m_pThreeDimensionalShape->input_nmm.max_boundarypoint_footpoint_idx == i)
			glPointSize(10.);
		else
			glPointSize(2.);
		glBegin(GL_POINTS);
			glVertex3d(m_pThreeDimensionalShape->input_nmm.BoundaryPoints[i][0],m_pThreeDimensionalShape->input_nmm.BoundaryPoints[i][1],m_pThreeDimensionalShape->input_nmm.BoundaryPoints[i][2]);
			glVertex3d(m_pThreeDimensionalShape->input_nmm.BoundaryPoints[i].fp[0],m_pThreeDimensionalShape->input_nmm.BoundaryPoints[i].fp[1],m_pThreeDimensionalShape->input_nmm.BoundaryPoints[i].fp[2]);
		glEnd();
	}
	

	glEnable(GL_LINE_SMOOTH);
	glColor3f(.8f,.8f,.8f);
	glLineWidth(2.);
	glBegin(GL_LINES);
	for(unsigned int i = 0; i < m_pThreeDimensionalShape->input_nmm.BoundaryPoints.size(); i ++)
	{
		glVertex3d(m_pThreeDimensionalShape->input_nmm.BoundaryPoints[i][0],m_pThreeDimensionalShape->input_nmm.BoundaryPoints[i][1],m_pThreeDimensionalShape->input_nmm.BoundaryPoints[i][2]);
		glVertex3d(m_pThreeDimensionalShape->input_nmm.BoundaryPoints[i].fp[0],m_pThreeDimensionalShape->input_nmm.BoundaryPoints[i].fp[1],m_pThreeDimensionalShape->input_nmm.BoundaryPoints[i].fp[2]);
	}
	glEnd();
	glDisable(GL_LINE_SMOOTH);

	return;
}

void GLWidget::gl_draw_ma_spheres()
{
	if(m_display_ma_spheres == false)
		return;
	if(m_pThreeDimensionalShape == NULL)
		return;

	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glColor3f(.3f, .5f, .8f);

	for(unsigned i = 0; i < m_pThreeDimensionalShape->input_nmm.vertices.size(); i ++)
	{
		if(m_pThreeDimensionalShape->input_nmm.vertices[i].first)
		{
			Vector3d pos = m_pThreeDimensionalShape->input_nmm.vertices[i].second->sphere.center;
			double radius = m_pThreeDimensionalShape->input_nmm.vertices[i].second->sphere.radius;
			glTranslated(pos[0],pos[1],pos[2]);
			gluSphere(spherequadric,radius,50,50);
			glTranslated(-pos[0],-pos[1],-pos[2]);
		}
	}
	glDisable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);
	return;
}

void GLWidget::gl_draw_ma_envelopes()
{
	if(m_display_ma_envelopes == false)
		return;
	if(m_pThreeDimensionalShape == NULL)
		return;

	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glColor3f(.3f, .5f, .8f);
	Cone * pcone;
	
	for(unsigned i = 0; i < m_pThreeDimensionalShape->input_nmm.edges.size(); i ++)
	{
		if(m_pThreeDimensionalShape->input_nmm.edges[i].first)
		{
			if(!m_pThreeDimensionalShape->input_nmm.edges[i].second->validenvelope)
				continue;
			pcone = &(m_pThreeDimensionalShape->input_nmm.edges[i].second->cone);
			glPushMatrix();
			glTranslated(pcone->apex[0],pcone->apex[1],pcone->apex[2]);
			glRotated(-pcone->rot_angle,pcone->rot_axis[0],pcone->rot_axis[1],pcone->rot_axis[2]);
			gluCylinder(conequadric,pcone->base,pcone->top,pcone->height,120,1);
			glPopMatrix();
		}
	}

	glColor3f(.15f, .25f, .7f);
	glColor3f(.15f, .25f, .4f);
	glColor3f(.3f, .5f, .8f);

	for(unsigned i = 0; i < m_pThreeDimensionalShape->input_nmm.faces.size(); i ++)
	{
		if(!m_pThreeDimensionalShape->input_nmm.faces[i].first)
			continue;
		if(!m_pThreeDimensionalShape->input_nmm.faces[i].second->valid_st)
			continue;
		if(!m_pThreeDimensionalShape->input_nmm.faces[i].second->has_footpoints)
			continue;
		for(unsigned k = 0; k < 2; k ++)
		{
			SimpleTriangle * st;
			st = &(m_pThreeDimensionalShape->input_nmm.faces[i].second->st[k]);

			glBegin(GL_POLYGON);
				glNormal3dv(st->normal);
				glVertex3dv(st->v[0]);
				glVertex3dv(st->v[1]);
				glVertex3dv(st->v[2]);
			glEnd();
			glBegin(GL_POLYGON);
				glNormal3dv(st->normal);
				glVertex3dv(st->v[1]);
				glVertex3dv(st->v[0]);
				glVertex3dv(st->v[2]);
			glEnd();
		}
	}
	glDisable(GL_AUTO_NORMAL);

	glDisable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);
	return;
}

Vector3d GLWidget::HSVToRGB(Vector3d in)
{
	double      hh, p, q, t, ff;
	long        i;
	Vector3d         out;

	if(in.Y() <= 0.0) {       // < is bogus, just shuts up warnings
		out[0] = in[2];
		out[1] = in[2];
		out[2] = in[2];
		return out;
	}
	hh = in[0];
	if(hh >= 360.0) hh = 0.0;
	hh /= 60.0;
	i = (long)hh;
	ff = hh - i;
	p = in[2] * (1.0 - in[1]);
	q = in[2] * (1.0 - (in[1] * ff));
	t = in[2] * (1.0 - (in[1] * (1.0 - ff)));

	switch(i) {
	case 0:
		out[0] = in[2];
		out[1] = t;
		out[2] = p;
		break;
	case 1:
		out[0] = q;
		out[1] = in[2];
		out[2] = p;
		break;
	case 2:
		out[0] = p;
		out[1] = in[2];
		out[2] = t;
		break;

	case 3:
		out[0] = p;
		out[1] = q;
		out[2] = in[2];
		break;
	case 4:
		out[0] = t;
		out[1] = p;
		out[2] = in[2];
		break;
	case 5:
	default:
		out[0] = in[2];
		out[1] = p;
		out[2] = q;
		break;
	}
	return out;     
}

void GLWidget::gl_draw_facet(Facet_handle pFacet, double factor)
{
	Wm4::Vector3d cent(Wm4::Vector3d::ZERO);
	int degree(0);
	
	Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
	do
	{
		cent += to_wm4(pHalfedge->vertex()->point());
		degree ++;
	}
	while( ++pHalfedge != pFacet->facet_begin() );
	cent = cent / degree;
	
	if(m_smoothshading)
	{
		Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
		if(m_display_reverse_orientation == false)
		{
			//do
			//{
			//	glNormal3d(pHalfedge->vertex()->normal[0], pHalfedge->vertex()->normal[1], pHalfedge->vertex()->normal[2]);
			//	double per;
			//	if (m_display_slab_edge)
			//	{
			//		//per = pHalfedge->vertex()->slab_hausdorff_dist / (m_pThreeDimensionalShape->input.bb_diagonal_length * 0.020);
			//		//per = pHalfedge->vertex()->slab_hausdorff_dist / m_pThreeDimensionalShape->slab_mesh.maxhausdorff_distance;

			//		Vector3d rgb;


			//		//rgb = HSVToRGB(Vector3d(240 * (1 - per), 1.0, 1.0));
			//		//if (per <= 0.1)
			//		//{
			//		//	rgb = HSVToRGB(Vector3d(240 * (1 - per / 0.3) + 90 * (per / 0.3), 1.0, 1.0));
			//		//} 
			//		//else
			//		//{
			//		//	rgb = HSVToRGB(Vector3d(90 - 90 * per, 1.0, 1.0));
			//		//}

			//		//glColor3f(rgb[0], rgb[1], rgb[2]);


			//		float mc[3];
			//		m_ColorRamp.RedGreenBlue((pHalfedge->vertex()->slab_hausdorff_dist)*255. / 0.02, mc);
			//		//m_ColorRamp.RedGreenBlue((pHalfedge->vertex()->slab_hausdorff_dist)*255. / m_pThreeDimensionalShape->slab_mesh.maxhausdorff_distance, mc);
			//		glColor3fv(mc);


			//	}
			//	glVertex3d(
			//		(pHalfedge->vertex()->point()[0]* factor + cent[0] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length,
			//		(pHalfedge->vertex()->point()[1]* factor + cent[1] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length,
			//		(pHalfedge->vertex()->point()[2]* factor + cent[2] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length);
			//} while( ++pHalfedge != pFacet->facet_begin() );

			do
			{
				//glColor3f(1.0f, 0.95f, 0.95f);
				//glColor3f(0.9f, 0.9f, 0.9f);
				//glColor3f(.3f, .5f, .8f);
				//glColor3f(0.75f, 0.75f, 0.75f);
				//glColor3f(0.65f, 0.65f, 0.65f);
				glColor3f(136.0 / 255.0, 134.0 / 255.0, 119.0 / 255.0);
				glNormal3d(pHalfedge->vertex()->normal[0], pHalfedge->vertex()->normal[1], pHalfedge->vertex()->normal[2]);
				double per;
				if (m_display_slab_edge)
				{
					per = pHalfedge->vertex()->slab_hausdorff_dist / m_pThreeDimensionalShape->input.bb_diagonal_length;
					//per = pHalfedge->vertex()->slab_hausdorff_dist / m_pThreeDimensionalShape->slab_mesh.maxhausdorff_distance;
				}
				glVertex3d(
					(pHalfedge->vertex()->point()[0]* factor + cent[0] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length,
					(pHalfedge->vertex()->point()[1]* factor + cent[1] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length,
					(pHalfedge->vertex()->point()[2]* factor + cent[2] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length);
			} while( ++pHalfedge != pFacet->facet_begin() );
		}else
		{
			do
			{
				//glColor3f(1.0f, 0.95f, 0.95f);
				//glColor3f(0.9f, 0.9f, 0.9f);
				//glColor3f(.3f, .5f, .8f);
				glColor3f(0.75f, 0.75f, 0.75f);
				//glColor3f(0.65f, 0.65f, 0.65f);
				glNormal3d(-pHalfedge->vertex()->normal[0], -pHalfedge->vertex()->normal[1], -pHalfedge->vertex()->normal[2]);
				double per;
				if (m_display_slab_edge)
				{
					per = pHalfedge->vertex()->slab_hausdorff_dist / m_pThreeDimensionalShape->input.bb_diagonal_length;
					//per = pHalfedge->vertex()->slab_hausdorff_dist / m_pThreeDimensionalShape->slab_mesh.maxhausdorff_distance;
				}
				glVertex3d(
					(pHalfedge->vertex()->point()[0]* factor + cent[0] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length,
					(pHalfedge->vertex()->point()[1]* factor + cent[1] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length,
					(pHalfedge->vertex()->point()[2]* factor + cent[2] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length);
			} while( --pHalfedge != pFacet->facet_begin() );
		}
	}
	else
	{
		Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
		if(m_display_reverse_orientation == false)
		{
			glNormal3d(pFacet->normal[0], pFacet->normal[1], pFacet->normal[2]);
			do
			{
				glVertex3d(
					(pHalfedge->vertex()->point()[0]* factor + cent[0] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length,
					(pHalfedge->vertex()->point()[1]* factor + cent[1] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length,
					(pHalfedge->vertex()->point()[2]* factor + cent[2] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length);
			} while( ++pHalfedge != pFacet->facet_begin() );
		}
		else
		{
			glNormal3d(-pFacet->normal[0], -pFacet->normal[1], -pFacet->normal[2]);
			glColor3f(0.75f, 0.75f, 0.75f);
			//glColor3f(0.65f, 0.65f, 0.65f);
			//glColor3f(.3f, .5f, .8f);
			do
			{
				glVertex3d(
					(pHalfedge->vertex()->point()[0]* factor + cent[0] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length,
					(pHalfedge->vertex()->point()[1]* factor + cent[1] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length,
					(pHalfedge->vertex()->point()[2]* factor + cent[2] * (1. - factor)) / m_pThreeDimensionalShape->input.bb_diagonal_length);
			} while( --pHalfedge != pFacet->facet_begin() );
		}
	}
	
}

void GLWidget::gl_adjust_view(ThreeDimensionalShape * pThreeDimensionalShape)
{
	if(pThreeDimensionalShape == NULL)
		return;
	Mesh * pmesh = &(pThreeDimensionalShape->input);


	if(pmesh == NULL)
		return;
	
	double xmax = pmesh->m_max[0] / pmesh->bb_diagonal_length;
	double ymax = pmesh->m_max[1] / pmesh->bb_diagonal_length;
	double zmax = pmesh->m_max[2] / pmesh->bb_diagonal_length;
	
	double xmin = pmesh->m_min[0] / pmesh->bb_diagonal_length;
	double ymin = pmesh->m_min[1] / pmesh->bb_diagonal_length;
	double zmin = pmesh->m_min[2] / pmesh->bb_diagonal_length;
	
	double gsize[3];
	gsize[0] = xmax - xmin;
	gsize[1] = ymax - ymin;
	gsize[2] = zmax - zmin;
	
	double mscale = gsize[0]> gsize[1]?gsize[0]:gsize[1];
	mscale = mscale>gsize[2]?mscale:gsize[2];

	if(mscale != 0)
	{
		mscale = 10. / mscale;
		mscale *= scalefactor;
		glScaled(mscale, mscale, mscale);
	}
	
	glTranslated( -0.5*(xmax+xmin), -0.5*(ymax+ymin), -0.5*(zmax+zmin) );
	oldscale = mscale;
	oldtranslation[0] = -0.5*(xmax+xmin);
	oldtranslation[1] = -0.5*(ymax+ymin);
	oldtranslation[2] = -0.5*(zmax+zmin);
}

void GLWidget::gl_adjust_ma(ThreeDimensionalShape * pThreeDimensionalShape)
{

	// BY Zhiyang...
	// Already scaled...


	if(pThreeDimensionalShape == NULL)
		return;
	SlabMesh * pmesh = &(pThreeDimensionalShape->slab_mesh);


	if(pmesh == NULL)
		return;

	double xmax = pmesh->m_max[0];
	double ymax = pmesh->m_max[1];
	double zmax = pmesh->m_max[2];

	double xmin = pmesh->m_min[0];
	double ymin = pmesh->m_min[1];
	double zmin = pmesh->m_min[2];

	double gsize[3];
	gsize[0] = xmax - xmin;
	gsize[1] = ymax - ymin;
	gsize[2] = zmax - zmin;

	double mscale = gsize[0]> gsize[1]?gsize[0]:gsize[1];
	mscale = mscale>gsize[2]?mscale:gsize[2];

	if(mscale != 0)
	{
		mscale = 10. / mscale;
		mscale *= scalefactor;
		glScaled(mscale, mscale, mscale);
	}

	glTranslated( -0.5*(xmax+xmin), -0.5*(ymax+ymin), -0.5*(zmax+zmin) );
	oldscale = mscale;
	oldtranslation[0] = -0.5*(xmax+xmin);
	oldtranslation[1] = -0.5*(ymax+ymin);
	oldtranslation[2] = -0.5*(zmax+zmin);
}

void GLWidget::gl_draw_headbar()
{
	if(m_display_headbar == false)
		return;
	// draw menu bar
	glDisable(GL_LIGHTING);
	int nViewport[4];
	glGetIntegerv(GL_VIEWPORT, nViewport);
	int nWidth = nViewport[2];
	int nHeight = nViewport[3];

	glDisable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);

	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, nWidth, 0, nHeight,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushMatrix();
	glPopMatrix();
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_DEPTH_TEST);

	glColor3f(0.0f, 0.0f, 0.0f);
	QFont serifFont("Arial", 10, QFont::Normal);

	glEnable(GL_LIGHTING);

	if(m_pThreeDimensionalShape == NULL)
		return;

	double len[4];
	len[0] = m_pThreeDimensionalShape->input.m_max[0] - m_pThreeDimensionalShape->input.m_min[0];
	len[1] = m_pThreeDimensionalShape->input.m_max[1] - m_pThreeDimensionalShape->input.m_min[1];
	len[2] = m_pThreeDimensionalShape->input.m_max[2] - m_pThreeDimensionalShape->input.m_min[2];
	len[3] = sqrt(len[0]*len[0]+len[1]*len[1]+len[2]*len[2]);

	this->renderText(
		5,
		height() - 85,
		QString(tr("Mesh: Num_V: %1  Num_E: %2  Num_Tri: %3"))
		.arg(m_pThreeDimensionalShape->input.size_of_vertices())
		.arg(m_pThreeDimensionalShape->input.size_of_halfedges()/2)
		.arg(m_pThreeDimensionalShape->input.size_of_facets())
		.arg(len[0])
		.arg(len[1])
		.arg(len[2])
		,
		serifFont);

	this->renderText(
		5,
		height() - 55,
		QString(tr("Initial Medial Mesh : Num_V: %1  Num_E: %2  Num_Tri: %3"))
		.arg(m_pThreeDimensionalShape->slab_mesh.iniNumVertices)
		.arg(m_pThreeDimensionalShape->slab_mesh.iniNumEdges)
		.arg(m_pThreeDimensionalShape->slab_mesh.iniNumFaces)
		.arg(len[0])
		.arg(len[1])
		.arg(len[2])
		,
		serifFont);

	//this->renderText(
	//		5,
	//		height() - 125,
	//		QString(tr("Sphere Mesh: Num_V: %1  Num_E: %2  Num_Tri: %3"))
	//		.arg(m_pThreeDimensionalShape->input_sphere_mesh.numVertices)
	//		.arg(m_pThreeDimensionalShape->input_sphere_mesh.numEdges)
	//		.arg(m_pThreeDimensionalShape->input_sphere_mesh.numFaces)
	//		.arg(len[0])
	//		.arg(len[1])
	//		.arg(len[2])
	//		,
	//		serifFont);

	//	this->renderText(
	//		5,
	//		height() - 105,
	//		QString(tr("QEM Mesh: Num_V: %1  Num_E: %2  Num_Tri: %3"))
	//		.arg(m_pThreeDimensionalShape->qem_mesh.numVertices)
	//		.arg(m_pThreeDimensionalShape->qem_mesh.numEdges)
	//		.arg(m_pThreeDimensionalShape->qem_mesh.numFaces)
	//		.arg(len[0])
	//		.arg(len[1])
	//		.arg(len[2])
	//		,
	//		serifFont);

	//	this->renderText(
	//		5,
	//		height() - 85,
	//		QString(tr("QA QEM Mesh: Num_V: %1  Num_E: %2  Num_Tri: %3"))
	//		.arg(m_pThreeDimensionalShape->ma_qem_mesh.numVertices)
	//		.arg(m_pThreeDimensionalShape->ma_qem_mesh.numEdges)
	//		.arg(m_pThreeDimensionalShape->ma_qem_mesh.numFaces)
	//		.arg(len[0])
	//		.arg(len[1])
	//		.arg(len[2])
	//		,
	//		serifFont);.	
	
	//if (m_pThreeDimensionalShape->slab_initial)
	//{
	//	EdgeInfo topEdge = m_pThreeDimensionalShape->slab_mesh.edge_collapses_queue.top();
	//	this->renderText(
	//		5,
	//		height() - 125,
	//		QString(tr("Next Simplify Cost: %1  "))
	//		.arg(topEdge.collapse_cost)
	//		.arg(len[0])
	//		,
	//		serifFont);
	//
	//	if (m_pThreeDimensionalShape->slab_mesh.edges[topEdge.edge_num].first)
	//	{
	//		this->renderText(
	//			5,
	//			height() - 105,
	//			QString(tr("Next Ratio between Hyperbolic and Euclid: %1"))
	//			.arg(m_pThreeDimensionalShape->slab_mesh.edges[topEdge.edge_num].second->hyperbolic_weight)
	//			.arg(len[0])
	//			,
	//			serifFont);

	//		this->renderText(
	//			5,
	//			height() - 85,
	//			QString(tr("Next QEM error: %1"))
	//			.arg((topEdge.collapse_cost / m_pThreeDimensionalShape->slab_mesh.edges[topEdge.edge_num].second->hyperbolic_weight) - m_pThreeDimensionalShape->slab_mesh.k)
	//			.arg(len[0])
	//			,
	//			serifFont);
	//	}
	//}

		this->renderText(
			5,
			height() - 25,
			QString(tr("Simplified Medial Mesh: Num_V: %1  Num_E: %2  Num_Tri: %3"))
			.arg(m_pThreeDimensionalShape->slab_mesh.numVertices)
			.arg(m_pThreeDimensionalShape->slab_mesh.numEdges)
			.arg(m_pThreeDimensionalShape->slab_mesh.numFaces)
			.arg(len[0])
			.arg(len[1])
			.arg(len[2])
			,
			serifFont);

		//this->renderText(
		//	5,
		//	height() - 45,
		//	QString(tr("Hausdorff: VQEM: %1	 Slab: %2"))
		//	//.arg(m_pThreeDimensionalShape->slab_mesh.maxhausdorff_distance / m_pThreeDimensionalShape->input.bb_diagonal_length)
		//	.arg(m_pThreeDimensionalShape->slab_mesh.maxhausdorff_distance)
		//	.arg(len[0])
		//	,
		//	serifFont);

		//this->renderText(
		//	5,
		//	height() - 25,
		//	QString(tr("Mean_Hausdorff_Distance: VQEM: %1	 Slab: %2"))
		//	.arg(m_pThreeDimensionalShape->slab_mesh.meanhausdorff_distance / m_pThreeDimensionalShape->input.bb_diagonal_length)
		//	.arg(len[0])
		//	,
		//	serifFont);
		//
		//this->renderText(
		//	5,
		//	height() - 5,
		//	QString(tr("Max_Mean_Square_Error: VQEM: %1	 Slab: %2"))
		//	.arg(sqrt(m_pThreeDimensionalShape->slab_mesh.max_mean_squre_error) / m_pThreeDimensionalShape->input.bb_diagonal_length)
		//	.arg(len[0])
		//	,
		//	serifFont);
	return;
}

void GLWidget::gl_draw_colorbar()
{
	if(m_display_colorbar == false)
		return;

	int nViewport[4];
	glGetIntegerv(GL_VIEWPORT, nViewport);

	int nWidth = nViewport[2], nHeight = nViewport[3];
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, nWidth, 0, nHeight,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushMatrix();

	int color_height=nWidth/3;
	int color_width=nHeight/30;

	int color_base_x = 30;
	int color_base_y = nHeight/2-color_height/2;
	float mc[3];
	for(int i=0; i<=color_height;i++)
	{
		int l = int(i*255.0/color_height);
		m_ColorRamp.RedGreenBlue(l,mc);

		glColor3fv(mc);
		glBegin(GL_LINES);
		glVertex2i(color_base_x, color_base_y+i);
		glVertex2i(color_base_x+color_width, color_base_y+i);
		glEnd();
	}

	QFont serifFont("Arial", 25, QFont::Bold);

	glColor3f(0.0f, 0.0f, 0.0f);

	//renderText(color_base_x+color_width+10, color_base_y+color_height, QString(tr("%1")).arg(m_Colorbar_MinValue), serifFont);
	//renderText(color_base_x+color_width+10, color_base_y+color_height/2, QString(tr("%1")).arg(0.5*(m_Colorbar_MinValue+m_Colorbar_MaxValue)), serifFont);
	//renderText(color_base_x+color_width+10, color_base_y, QString(tr("%1")).arg(m_Colorbar_MaxValue), serifFont);
	renderText(color_base_x+color_width+10, color_base_y+color_height, QString(tr("%1")).arg(0), serifFont);
	renderText(color_base_x+color_width+10, color_base_y+color_height/2, QString(tr("%1")).arg(0.5), serifFont);
	renderText(color_base_x+color_width+10, color_base_y, QString(tr("%1")).arg(1), serifFont);

	glPopMatrix();
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_DEPTH_TEST);	
	glEnable(GL_LIGHTING);
}

//处理鼠标选择，左键触发
void GLWidget::processLeftButton(int xPos, int yPos)
{
	GLfloat fAspect;

	//点击计数器和视口存储
	GLint viewPort[4];

	//设置缓冲区
	glSelectBuffer(BUFFER_LENGTH, selectBuffer);

	//获得视口
	glGetIntegerv(GL_VIEWPORT,viewPort);

	//修改渲染模式
	glRenderMode(GL_SELECT);

	glInitNames();     //初始化名字栈
	glPushName(0);     //在名字栈中放入一个初始化名字，这里为‘0’

	//切换到投影模式，并保存矩阵
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();

	//围绕鼠标点(xPos，yPos)建立新的单位立方体裁剪区域
	//并在垂直和水平方向扩展M个像素
	glLoadIdentity();
	gluPickMatrix(xPos,viewPort[3]-yPos,5.0,5.0,viewPort);

	//应用透视投影矩阵
	fAspect = (float)viewPort[2]/(float)viewPort[3];
	gluPerspective(45.0f,fAspect,1.0,425.0);

	//绘制场景
	//display();
	//ply2Model->DrawModel(selected);
	//updateGL();
	//gl_draw_vertices(m_pThreeDimensionalShape);
	gl_draw_slab_vertex(m_pThreeDimensionalShape);

	//收集点击记录,hits为选中点的个数
	hits = glRenderMode(GL_RENDER);

	//当有点被选中时，
	if (hits != 0)
	{
		//hits为被选中点的个数
		std::cout << "hits=" << hits << std::endl;

		//选择深度最小的点为选中目标点，
		//selectBuffer[]中分别填充的是：
		//stelect[0]:被选中标识
		//select[1]:最小深度
		//select[2]:最大深度
		//select[3]:标号ID
		int idOfSelectedPoint=0; //选中点的ID
		//int i=0;
		//int tempi = 0;
		//for (i=1;i<hits;i++)
		//{
		//	float temp = (float)selectBuffer[1];//第一个点的最小深度

		//	if (temp > selectBuffer[i*4+1])
		//	{
		//		temp = selectBuffer[i*4+1];
		//		tempi = i*4+1;
		//	}
		//}
		idOfSelectedPoint = selectBuffer[3];

		std::cout << "idOfSelectedPoint=" << idOfSelectedPoint << std::endl;
		char temp[20];
		sprintf(temp, "%d", idOfSelectedPoint);
		string idstring(temp);
		//QMessageBox::StandardButton reply;
		//reply = QMessageBox::information(this, tr("Time"), tr(idstring.c_str()));


		int i = selectBuffer[3];
		double min_coll = DBL_MAX;
		int index = 0;
		set<unsigned> eset = m_pThreeDimensionalShape->slab_mesh.vertices[i].second->edges_;
		for (set<unsigned>::iterator si = eset.begin(); si != eset.end(); si++)
		{
			if (m_pThreeDimensionalShape->slab_mesh.edges[*si].second->collapse_cost < min_coll)
			{
				index = *si;
				min_coll = m_pThreeDimensionalShape->slab_mesh.edges[*si].second->collapse_cost;
			}
		}		
		
		sprintf(temp, "%f", min_coll * 1e5);
		string collstring(temp);

		double ratio = m_pThreeDimensionalShape->slab_mesh.GetRatioHyperbolicEuclid(index);
		sprintf(temp, "%f", ratio);
		string ratiostring(temp);

		double qem = m_pThreeDimensionalShape->slab_mesh.edges[index].second->qem_error;
		sprintf(temp, "%f", qem * 1e5);
		string qemstring(temp);

		unsigned rrank = 1;
		unsigned crank = 1;
		unsigned qrank = 1;
		for (unsigned ri = 0; ri < m_pThreeDimensionalShape->slab_mesh.edges.size(); ri++)
		{
			if (m_pThreeDimensionalShape->slab_mesh.edges[ri].first)
			{
				if (m_pThreeDimensionalShape->slab_mesh.edges[ri].second->hyperbolic_weight < ratio)
					rrank++;
				if (m_pThreeDimensionalShape->slab_mesh.edges[ri].second->collapse_cost < min_coll)
					crank++;
				if (m_pThreeDimensionalShape->slab_mesh.edges[ri].second->qem_error < qem)
					qrank++;
			}
		}

		sprintf(temp, "%d", rrank);
		string rrankstring(temp);

		sprintf(temp, "%d", crank);
		string crankstring(temp);

		sprintf(temp, "%d", qrank);
		string qrankstring(temp);

		unsigned v1 = m_pThreeDimensionalShape->slab_mesh.edges[index].second->vertices_.first;
		unsigned v2 = m_pThreeDimensionalShape->slab_mesh.edges[index].second->vertices_.second;
		int relate_face = m_pThreeDimensionalShape->slab_mesh.vertices[v1].second->related_face + m_pThreeDimensionalShape->slab_mesh.vertices[v2].second->related_face;

		sprintf(temp, "%d", relate_face);
		string facestring(temp);

		sprintf(temp, "%d", index);
		string edgestring(temp);

		QMessageBox::StandardButton reply;
		reply = QMessageBox::information(this, tr("Selected Info"), tr("VNum: ") + tr(idstring.c_str()) + 
			//tr("; Related Face: ") + tr(facestring.c_str()) + 			
			tr("; Edge number: ") + tr(facestring.c_str()) + 
			tr("; Min_Edge_Coll: ") + tr(collstring.c_str()) + 
			tr("; Min_Edge_Weight: ") + tr(ratiostring.c_str()) + 
			tr("; Min_Edge_QEM: ") + tr(qemstring.c_str()) + 
			tr("; ratio rank: ") + tr(rrankstring.c_str()) + 
			tr("; qem rank: ") + tr(qrankstring.c_str()) + 
			tr("; cost rank: ") + tr(crankstring.c_str()));

	}

	//回复投影矩阵
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	//恢复到模型视图，以便进行正常的渲染
	glMatrixMode(GL_MODELVIEW);
}

void GLWidget::SaveViewpoint(const char * fout)
{
	std::ofstream viewout(fout);
	viewout << oldscale << std::endl;
	viewout << oldtranslation[0] << std::endl;
	viewout << oldtranslation[1] << std::endl;
	viewout << oldtranslation[2] << std::endl;
	viewout << scalefactor << std::endl;
	viewout << xRot << std::endl;
	viewout << yRot << std::endl;
	viewout << zRot << std::endl;
	viewout << xTrans << std::endl;
	viewout << yTrans << std::endl;
	viewout.close();
}

void GLWidget::LoadViewpoint(const char * fin)
{
	std::ifstream viewin(fin);
	viewin >> oldscale;
	viewin >> oldtranslation[0];
	viewin >> oldtranslation[1];
	viewin >> oldtranslation[2];
	viewin >> scalefactor;
	viewin >> xRot;
	viewin >> yRot;
	viewin >> zRot;
	viewin >> xTrans;
	viewin >> yTrans;
	viewin.close();
}
