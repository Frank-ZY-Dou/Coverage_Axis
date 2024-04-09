#ifndef CAMERA_H_
#define CAMERA_H_

#include <CGAL/glu.h>
#include <string>
#include <fstream>
#include <cmath>
#include <cassert>
#include "LinearAlgebra/Wm4Vector.h"
#include "LinearAlgebra/Wm4Matrix.h"



using namespace std;

class Camera
{
public:
	Wm4::Vector3d focus;//对于屏幕中心点的世界坐标
	Wm4::Vector3d translate;//平移
	GLint viewport[4];//窗口

	GLint x;
	GLint y;
	GLsizei width;
	GLsizei height;

	GLdouble zmin;
	GLdouble zmax;
	GLdouble xmin;
	GLdouble xmax;
	GLdouble ymin;
	GLdouble ymax;
	GLdouble len;//包围盒半径

	Wm4::Matrix4d modelview_matrix;

	GLdouble scale;

public:
	Camera();
	double getNear(){return -zmax*scale;};
	double getFar(){return -zmin*scale;};
	void setBox(Wm4::Vector3d pmin,Wm4::Vector3d pmax);
	void setBox(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax);
	void SetViewPort(int x,int y,int width,int height);
	void Translate(double startx,double starty,double endx,double endy);//鼠标点击的两个屏幕坐标.
	void Scale(double delta);//滚轮放大.
	void Scale(double deltax, double deltay, double deltaz);
	Wm4::Vector3d unProject(int x, int y, int z);
	Wm4::Matrix4d Rotate(double startx,double starty,double endx,double endy);
	void SetLookAt(const Wm4::Vector3d& Eye,const Wm4::Vector3d& ViewUp);
	void Setup();
	void SetupModelView();
	void getBox(double& l,double& r,double& b,double& t);
	void addMatrix(Wm4::Matrix4d mat);
	void Reset();
};

#endif