#ifndef _GEOMETRYOBJECTS_H
#define _GEOMETRYOBJECTS_H

#include "LinearAlgebra/Wm4Matrix.h"
#include "LinearAlgebra/Wm4Vector.h"

using namespace Wm4;

class SimpleTriangle
{
public:
	SimpleTriangle():normal(Vector3d(0, 0, 0)){};

public:
	Vector3d v[3];
	Vector3d normal;
	bool ProjectOntoSimpleTriangle(const Vector3d & p, Vector3d & fp, double & dist);
	void UpdateNormal(bool reverse = false);
};

class Sphere
{
public:
	Sphere(){}
	Sphere(Vector3d center_, double radius_) {center = center_; radius = radius_;}
public:
	Vector3d center;
	double radius;

	void ProjectOntoSphere(const Vector3d & p, Vector3d & fp, double & signeddist);
	double DistanceToSphere(const Vector3d & p);
public:
	Sphere operator*(const double d)
	{
		return Sphere(d*center,d*radius);
	}
	Sphere operator+(const Sphere & other)
	{
		return Sphere(center + other.center, radius + other.radius);
	}
};

class Cone
{
public:
	Vector3d apex;
	Vector3d axis;
	double base;
	double top;
	double height;

	// 1 -- height = 0.0;
	// 2 -- cylinder
	// 3 -- general regular cone
	int type;

	// the rotation axis to rotate the axis to z-axis
	Vector3d rot_axis;
	double rot_angle;
public:
	Cone(){}
	Cone(Wm4::Vector3d c0, double r0, Wm4::Vector3d c1, double r1);
	bool ProjectOntoCone(const Vector3d & p, Vector3d & fp, double & signeddist);
	Sphere BoundingSphere();

	std::vector< Sphere > SampleSpheres(unsigned num);
};

class SamplePoint : public Wm4::Vector3d
{
public:
	SamplePoint() {fp[0] = X(); fp[1] = Y(); fp[2] = Z(); onregion = 0;};
	SamplePoint(double x,double y,double z) : Vector3d(x,y,z) {fp = Wm4::Vector3d(x,y,z);};
	~SamplePoint() {};

public:
	Wm4::Vector3d fp;
	double fpdist;
	int tag;
	std::set<unsigned int> adjacentmedialpointlist;
	unsigned int connecting_medialpoint;
	int connecting_edge[2];
	int connecting_tri[3];
	// 0 -- sphere
	// 1 -- cone
	// 2 -- plane
	int onregion;

public:
	int tag_int;
	double tag_double;
	std::set<unsigned int> tag_ui_set;
	unsigned tag_ui;

};
class Triangle
{
public:
	Triangle();
	Triangle(Vector3d v0, Vector3d v1, Vector3d v2) {v[0] = v0; v[1] = v1; v[2] = v2;};
	~Triangle();
	double Area();
	Vector3d Centroid();
public:
	Vector3d v[3]; // three vertices
	Vector3d n;  // normal vector
	double d; // density

public:
	bool is_obtuse();
	double anglev(int i);
	double voronoi_area_Meyer(int i);
	double voronoi_area_new(int i);

public:
	double v_interpolationerror[3];
	double f_interpolationerror;
	double f_normalerror;
	int v_tag[3];
	int f_i_tag;
	int f_n_tag;
};

class AnisoTriangle : public Triangle
{
public:
	AnisoTriangle() : Triangle() {};
	AnisoTriangle(Vector3d v0, Vector3d v1, Vector3d v2, Matrix3d m0, Matrix3d m1, Matrix3d m2) : Triangle(v0,v1,v2) {m[0] = m0; m[1] = m1; m[2] = m2;};
public:
	Matrix3d m[3]; // metric defined on vertices 
	Matrix3d metric; // metric of the triangle 
};

class LineSegment
{
public:
	LineSegment();
	LineSegment(Vector3d v0, Vector3d v1) {v[0] = v0; v[1] = v1;};
	~LineSegment();
	double Length();
public:
	Vector3d v[2];
};

class AnisoLineSegment : public LineSegment
{
public:
	AnisoLineSegment() : LineSegment() {};
	AnisoLineSegment(Vector3d v0, Vector3d v1, Matrix3d m0, Matrix3d m1) : LineSegment(v0,v1) {m[0] = m0; m[1] = m1;};
public:
	Matrix3d m[2];
};

class PNPlane
{
public:
	Wm4::Vector3d p; // a point on the plane
	Wm4::Vector3d n; // the normal of the plane
	int id;
public:
	PNPlane(Wm4::Vector3d p_, Wm4::Vector3d n_){ p=p_; n=n_;}
	PNPlane(){}
};


struct DualEdge
{
	unsigned int v[2];
	bool operator < (const DualEdge & e) const
	{
		if(v[0] < e.v[0]) return true;
		else if(v[0] > e.v[0]) return false;
		
		if(v[1] < e.v[1]) return true;
		else if(v[1] > e.v[1]) return false;
		
		return false;
	}
	bool operator > (const DualEdge & e) const
	{
		if(v[0] > e.v[0]) return true;
		else if(v[0] < e.v[0]) return false;
	
		if(v[1] > e.v[1]) return true;
		else if(v[1] < e.v[1]) return false;
	
		return false;
	}
};

struct DualTriangle
{
	unsigned int v[3];
	double max_interpolationerror;
	double total_interpolationerror;
	double max_normalerror;
	double total_normalerror;
	bool operator < (const DualTriangle & t) const
	{
		if(v[0] < t.v[0]) return true;
		else if(v[0] > t.v[0]) return false;

		if(v[1] < t.v[1]) return true;
		else if(v[1] > t.v[1]) return false;
		
		if(v[2] < t.v[2]) return true;
		else if(v[2] > t.v[2]) return false;
		
		return false;
	}
	bool operator > (const DualTriangle & t) const
	{
		if(v[0] > t.v[0]) return true;
		else if(v[0] < t.v[0]) return false;
	
		if(v[1] > t.v[1]) return true;
		else if(v[1] < t.v[1]) return false;

		if(v[2] > t.v[2]) return true;
		else if(v[2] < t.v[2]) return false;
		
		return false;
	}
};

int clip(PNPlane & pnplane, Triangle t, std::vector<Triangle> & inner, std::vector<Triangle> & outer);
double TriangleArea(Vector3d &v0, Vector3d &v1, Vector3d &v2);
Vector3d ThreeAxesScaling(Vector3d & axis0, double & scale0, Vector3d & axis1, double & scale1, Vector3d & axis2, double & scale2, Vector3d & v);
//Vector3d TransformToEuclidean(Vector3d & axis0, double & scale0, Vector3d & axis1, double & scale1, Vector3d & axis2, double & scale2, Vector3d & p);

// project the point p onto the triangle v0-v1-v2, fp is the footpoint, dist is the distance between p and fp. 
// if fp lies in the triangle, return true; otherwise (on the edges or the vertices), return false
bool ProjectOntoLineSegment(const Vector3d & p, const Vector3d & v0, const Vector3d & v1, Vector3d & fp, double & dist);

// check whether p1, p2 lie on the same side of the line a, b, assuming that p1, p2, a and b lie in the same plane
bool SameSide(const Wm4::Vector3d & p1, const Wm4::Vector3d & p2, const Wm4::Vector3d & a, const Wm4::Vector3d & b);

// check whether p lies in the triangle v0-v1-v2
// p, v0, v1, v2 should be on the same plane 
bool InsideTriangle(const Wm4::Vector3d &p, const Wm4::Vector3d &v0, const Wm4::Vector3d &v1, const Wm4::Vector3d &v2);

// project the point onto the line segment v0-v1, fp is the footpoint, dist is the distance between p and fp.
// if p lies in the segment, return true; otherwise (v0 or v1), return false
bool ProjectOntoTriangle(const Vector3d & p, const Vector3d & v0, const Vector3d & v1, const Vector3d & v2, Vector3d & fp, double & dist);

bool obtuse_triangle(Vector3d & v0, Vector3d & v1, Vector3d & v2);
	
bool acute_triangle(Vector3d & v0, Vector3d & v1, Vector3d & v2);

double cotan(Vector3d v, Vector3d v1, Vector3d v2);

double angle_from_cotan(Vector3d v, Vector3d v0, Vector3d v1);

// project the point p onto the ellipsoid x^2 / a^2 + y^2 / b^2 + z^2 / c^2 = 1
Wm4::Vector3d ProjectPointOntoEllipsoid(const double a, const double b, const double c, const Wm4::Vector3d & p);

// project the point p onto the elliptic cone x^2 / a^2 + z^2 / c^2 = y^2 / b^2
Wm4::Vector3d ProjectPointOntoEllipticCone(const double a, const double b, const double c, const Wm4::Vector3d & p);

/*
// project the point p onto the ellipsoid x^2 / a^2 + y^2 / b^2 + z^2 / c^2 = 1
void ProjectPointOntoEllipsoid(double a, double b, double c, Wm4::Vector3d & p);
*/

bool DifferentialInfoOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p,
	Wm4::Vector3d & pu, Wm4::Vector3d & pv, Wm4::Vector3d & normdir, 
	double & maxc, double & minc, Wm4::Vector3d & maxdir, Wm4::Vector3d & mindir);

double GaussianCurvatureOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p);
double MeanCurvatureOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p);
double MaxCurvatureOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p);
double MinCurvatureOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p);
Wm4::Vector3d NormalDirectionOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p);
Wm4::Vector3d MaxPrincipalDirectionOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p);
Wm4::Vector3d MinPrincipalDirectionOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p);


Wm4::Vector3d LinearInterpolation(Wm4::Vector3d & originalpoint, Wm4::Vector3d & offsetpoint, double offsetfactor);


Wm4::Vector3d RandomPointonSphere();

bool FastNoIntersect( const Wm4::Vector3d & p, const Wm4::Vector3d & dir, const Wm4::Vector3d & q, const Wm4::Vector3d & norm);

bool RayTriangleIntersect(const Wm4::Vector3d & p, const Wm4::Vector3d & dir, const Wm4::Vector3d & v0, const Wm4::Vector3d & v1, const Wm4::Vector3d & v2);
bool RayTriangleIntersectv2(const Wm4::Vector3d & p, const Wm4::Vector3d & dir, const Wm4::Vector3d & v0, const Wm4::Vector3d & v1, const Wm4::Vector3d & v2);



// the distance from p to the line v0v1
bool DistanceToLine(const Vector3d & p, const Vector3d & v0, const Vector3d & v1, double & dist, Vector3d & fp);

bool TriangleFromThreeSpheres(const Vector3d & c0,
							  const double & r0,
							  const Vector3d & c1,
							  const double & r1,
							  const Vector3d & c2,
							  const double & r2,
							  SimpleTriangle & st0,
							  SimpleTriangle & st1);

Vector3d TriangleNormal(const Vector3d & v0, const Vector3d & v1, const Vector3d & v2);

double VectorAngle(const Vector3d & v0, const Vector3d & v1);

#endif // _GEOMETRYOBJECTS_H
