// *** Camera.CPP ***
// Definitions for Camera class

#include "Camera.h"


Camera::Camera():translate(0,0,0),focus(0,0,0)
{
	width=0;
	height=0;

	zmin=0.0;
	zmax=0.0;
	xmin=0.0;
	xmax=0.0;
	ymin=0.0;
	ymax=0.0;
	len = 0;
	scale = 1.0;
	
}


void Camera::setBox(Wm4::Vector3d pmin, Wm4::Vector3d pmax)
{
	double temp;
	if (pmin.X() > pmax.X())
	{
		temp = pmin.X();
		pmin.X() = pmax.X();
		pmax.X() = temp;
	}
	if (pmin.Y() > pmax.Y())
	{
		temp = pmin.Y();
		pmin.Y() = pmax.Y();
		pmax.Y() = temp;
	}
	if (pmin.Z() > pmax.Z())
	{
		temp = pmin.Z();
		pmin.Z() = pmax.Z();
		pmax.Z() = temp;
	}

	//double fluff = 0.0001;
	double fluff = 0.1*(pmax-pmin).Length();
	//double fluff = 0.5;
	Wm4::Vector3d fluffer(fluff, fluff, fluff);

	//add in some fudge factor
	pmin = pmin + fluffer/(-1.0);
	pmax = pmax + fluffer;

	Wm4::Vector3d size = pmax - pmin;

	focus = pmin + size/2.0;

	len =size.Length()/2;

	xmin = focus.X() - len;
	xmax = focus.X() + len;
	ymin = focus.Y() - len;
	ymax = focus.Y() + len;
	zmin = focus.Z() - len;
	zmax = focus.Z() + len;

}

void Camera::setBox(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
{
	Wm4::Vector3d c1(xmin,ymin,zmin);
	Wm4::Vector3d c2(xmax,ymax,zmax);
	setBox(c1,c2);
}

void Camera::SetViewPort(int x, int y, int width, int height)
{
	this->x = (GLint)x;
	this->y = (GLint)y;
	this->width = (GLsizei)width;
	this->height = (GLsizei)height;
	viewport[0] = (GLint)x;
	viewport[1] = (GLint)y;
	viewport[2] = (GLsizei)width;
	viewport[3] = (GLsizei)height;
	glViewport(x,y,width,height);
}

void Camera::Translate(double startx,double starty,double endx,double endy)
{
	GLint startrealy,endrealy;
	startrealy = viewport[3]-(GLint)starty - 1;
	endrealy = viewport[3]-(GLint)endy -1;

	double proj_matrix[16];
	//Setup();
	glGetDoublev(GL_PROJECTION_MATRIX, proj_matrix);
	Wm4::Matrix4d modelView;
	//modelView.Scale(this->scale,this->scale,this->scale);
	Wm4::Vector3d start,end;
	
	gluUnProject((GLdouble)startx,(GLdouble)startrealy,0.0,modelView,proj_matrix,viewport,&start.X(),&start.Y(),&start.Z());
	gluUnProject((GLdouble)endx,(GLdouble)endrealy,0.0,modelView.M,proj_matrix,viewport,&end.X(),&end.Y(),&end.Z());

	mymath::Mat44d mat;
	
	mat.Translate(end - start);
	modelview_matrix = mat*modelview_matrix;//平移
	this->translate += end - start;
}

Wm4::Vector3d Camera::unProject(int x, int y, int z)
{
	Wm4::Vector3d p;
	GLint realy = viewport[3]-(GLint)y -1;

	double proj_matrix[16];
	
	glGetDoublev(GL_PROJECTION_MATRIX, proj_matrix);
	mymath::Mat44d modelView;
	//modelView.Scale(this->scale,this->scale,this->scale);

	gluUnProject((GLdouble)x,(GLdouble)realy,0.0,modelView.M,proj_matrix,viewport,&p.x,&p.y,&p.z);

	this->modelview_matrix.Inverse(modelView);
	p.z = z;
	Wm4::Vector3d posP = modelView * p;
	return posP;
}
void Camera::Scale(double delta)
{
	scale *= delta;
	Setup();
	Vec3d temp(focus + translate);

	mymath::Mat44d mat;
	mat.invTranslate(temp);
	modelview_matrix = mat*modelview_matrix;//平移到原点，因为是左乘，不需要颠倒顺序
	mat.Scale(delta,delta,delta);
	modelview_matrix = mat*modelview_matrix;//缩放
	mat.Translate(temp);
	modelview_matrix = mat*modelview_matrix;//平移回来
}

void Camera::Scale(double deltax, double deltay, double deltaz)
{
	Vec3d temp(focus + translate);
	mymath::Mat44d mat;
	mat.invTranslate(temp);
	modelview_matrix = mat*modelview_matrix;//平移到原点，因为是左乘，不需要颠倒顺序
	mat.Scale(deltax,deltay,deltaz);
	modelview_matrix = mat*modelview_matrix;//缩放
	mat.Translate(temp);
	modelview_matrix = mat*modelview_matrix;//平移回来
}
void Camera::addMatrix(mymath::Mat44d mat)
{
	this->modelview_matrix = mat * this->modelview_matrix;
}
mymath::Mat44d Camera::Rotate(double startx,double starty,double endx,double endy)
{
	startx = ((startx/((width-1)/2))-1);
	starty = -((starty/((height-1)/2))-1);	// OGL需要, Y 轴反向
	endx = ((endx/((width-1)/2))-1);
	endy = -((endy/((height-1)/2))-1);	// OGL需要, Y 轴反向

	//构造虚球体上对应的世界坐标点.
	double startlen = sqrt(startx*startx + starty*starty);
	double endlen = sqrt(endx*endx + endy*endy );

	Vec3d start,end;

	if(startlen > 1.0)
	{
		start.x = startx/startlen;
		start.y = starty/startlen;
		start.z = 0;
	}
	else
	{
		start.x = startx;
		start.y = starty;
		start.z = sqrt(1.0-startlen*startlen);
	}
	if(endlen > 1.0)
	{
		end.x = endx/endlen;
		end.y = endy/endlen;
		end.z = 0;
	}
	else
	{
		end.x = endx;
		end.y = endy;
		end.z = sqrt(1.0-endlen*endlen);
	}
	//求转动轴,
	start = start.Normalize();
	end = end.Normalize();
	Vec3d axis = start/end;//叉乘得到转轴.
	//求转动角度
	double angle = acos(start*end);//弧度
	angle = 180*angle/3.1415926;
	//先将物体平移到原点,position和focus类型不同意,以后修改

	Vec3d temp(focus + translate);

	mymath::Mat44d mat,resultM;
	mat.invTranslate(temp);
	resultM = mat;
	//modelview_matrix = mat*modelview_matrix;//平移到原点，因为是左乘，不需要颠倒顺序
	mat.Rotate(angle,axis);
	//modelview_matrix = mat*modelview_matrix;//旋转
	resultM = mat * resultM;
	mat.Translate(temp);
	//modelview_matrix = mat*modelview_matrix;//平移回去
	resultM = mat * resultM;
	modelview_matrix = resultM * modelview_matrix;
	return resultM;
}

void Camera::Reset() {
	modelview_matrix.Identity();
}

void Camera::SetLookAt(const Wm4::Vector3d& Eye, const Wm4::Vector3d& ViewUp)
{
	modelview_matrix.Identity();
	mymath::Mat44d mat;
	mat.LookAt(Eye,Wm4::Vector3d(0,0,0),ViewUp);

	Vec3d temp(focus + translate);
	mat.invTranslate(temp);
	modelview_matrix = mat * modelview_matrix;
	mat.LookAt(Eye,Vec3d(0,0,0),ViewUp);
	modelview_matrix = mat*modelview_matrix;
	mat.Translate(temp);
	modelview_matrix = mat*modelview_matrix;
	
	///*modelview_matrix.M[12] = -focus.x;
	//modelview_matrix.M[13] = -focus.y;
	//modelview_matrix.M[14] = -focus.z;*/
	//modelview_matrix.M[12] = 0;
	//modelview_matrix.M[13] = 0;
	//modelview_matrix.M[14] = 0;
	//this->translate = Vertex(0,0,0);
}
void Camera::SetupModelView()
{
	//glLoadIdentity();
	glMultMatrixd(modelview_matrix.M);
}
void Camera::getBox(double& l,double& r,double& b,double& t)
{
	GLdouble n = -zmax;
	GLdouble f = -zmin;

	if((xmax-xmin)/(ymax-ymin) > (GLdouble)width/(GLdouble)height)
	{
		l = xmin;
		r = xmax;
		b = (ymax+ymin)/2-(xmax-xmin)*(GLdouble)height/(2*(GLdouble)width);
		t = (ymax+ymin)/2+(xmax-xmin)*(GLdouble)height/(2*(GLdouble)width);
	}
	else
	{
		b = ymin;
		t = ymax;
		l = (xmax+xmin)/2-(ymax-ymin)*(GLdouble)width/(2*(GLdouble)height);
		r = (xmax+xmin)/2+(ymax-ymin)*(GLdouble)width/(2*(GLdouble)height);
	}
}
void Camera::Setup()
{
	glViewport(x,y,width,height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	GLdouble n = -zmax;
	GLdouble f = -zmin;
	GLdouble l,r,b,t;
	if((xmax-xmin)/(ymax-ymin) > (GLdouble)width/(GLdouble)height)
	{
		l = xmin;
		r = xmax;
		b = (ymax+ymin)/2-(xmax-xmin)*(GLdouble)height/(2*(GLdouble)width);
		t = (ymax+ymin)/2+(xmax-xmin)*(GLdouble)height/(2*(GLdouble)width);
	}
	else
	{
		b = ymin;
		t = ymax;
		l = (xmax+xmin)/2-(ymax-ymin)*(GLdouble)width/(2*(GLdouble)height);
		r = (xmax+xmin)/2+(ymax-ymin)*(GLdouble)width/(2*(GLdouble)height);
	}
	
	if(scale < 1)
		glOrtho(l,r,b,t,n,f);
	else
		glOrtho(l,r,b,t,n*scale*20,f*scale*20);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixd(modelview_matrix.M);
}