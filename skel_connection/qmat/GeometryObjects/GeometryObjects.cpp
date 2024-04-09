#include "GeometryObjects.h"

bool SimpleTriangle::ProjectOntoSimpleTriangle(const Vector3d & p, Vector3d & fp, double & dist)
{
	return ProjectOntoTriangle(p,v[0],v[1],v[2],fp,dist);
}

void SimpleTriangle::UpdateNormal(bool reverse)
{
	normal = TriangleNormal(v[0],v[1],v[2]);
	if (reverse)
	{
		normal = normal * -1;
	}
}

void Sphere::ProjectOntoSphere(const Vector3d & p, Vector3d & fp, double & signeddist)
{
	Vector3d ncp(p-center);
	double ncp_sl = ncp.SquaredLength();
	ncp.Normalize();
	fp = center + radius * ncp;
	signeddist = (ncp_sl < pow(radius,2))?-(p-fp).Length():(p-fp).Length();
}

double Sphere::DistanceToSphere(const Vector3d & p)
{
	Vector3d ncp(p-center);
	ncp.Normalize();
	return (center+radius*ncp-p).Length();
}

Cone::Cone(Wm4::Vector3d c0, double r0, Wm4::Vector3d c1, double r1)
{
	Vector3d c0c1 = c1-c0;
	// one sphere is included in another sphere
	if (c0c1.Length() - abs(r1 - r0) < 1e-8)
	{
		apex = r1 > r0 ? c1 : c0;
		axis = Vector3d(0,0,1);
		base = r1 > r0 ? r1 : r0;
		top = r1 > r0 ? r1 : r0;
		height = 0.;
		type = 1;
		rot_axis = Vector3d(0,0,1);
		rot_angle = 0.;
		return;
	}

	if(c0c1.Length() < 1e-8)
	{
		apex = c0;
		axis = Vector3d(0,0,1);
		base = r0;
		top = r0;
		height = 0.;
		type = 1;
		rot_axis = Vector3d(0,0,1);
		rot_angle = 0.;
	}
	else
	{
		double dr0r1 = fabs(r0-r1);
		if(dr0r1 < 1e-8)
		{
			apex = c0;
			axis = c1-c0;
			axis.Normalize();
			base = r0;
			top = r0;
			height = (c1-c0).Length();
			type = 2;
		}
		else
		{
			apex = (r1 * c0 - r0 * c1) / (r1 - r0);
			axis = (r0<r1)?(c1-c0):(c0-c1);
			axis.Normalize();

			double cangle;
			Vector3d apexc0 = c0 - apex;
			double vc0len = apexc0.Length();
			Vector3d apexc1 = c1 - apex;
			double vc1len = apexc1.Length();
			cangle = sqrt(1.-r0*r0/vc0len/vc0len);

			if(r0 < r1)
			{
				//cangle = sqrt(1.-r1*r1/vc1len/vc1len);
				apex = apex + apexc0 * cangle * cangle;
				base = r0 * cangle;
				top = r1 * cangle;
				height = (vc1len - vc0len) * cangle * cangle;
			}
			else
			{
				//cangle = sqrt(1.-r0*r0/vc0len/vc0len);
				apex = apex + apexc1 * cangle * cangle;
				base = r1 * cangle;
				top = r0 * cangle;
				height = (vc0len - vc1len) * cangle * cangle;
			}
			type = 3;
		}

		Vector3d za(0,0,1);
		rot_angle = acos(axis.Dot(za));
		if( (fabs(rot_angle) < 1e-12) || (fabs(rot_angle - Wm4::Mathd::PI) < 1e-12) )
			rot_axis = Vector3d(1,0,0);
		else
			rot_axis = axis.Cross(za);
		rot_axis.Normalize();
		rot_angle *= (180./Wm4::Mathd::PI);
	}
}

bool Cone::ProjectOntoCone(const Vector3d & p, Vector3d & fp, double & signeddist)
{
	Vector3d apexp;
	apexp = p - apex;
	Vector3d apexpCaxis = apexp.Cross(axis);

	Wm4::Vector3d v0,v1;
	// on the axis
	if(fabs(apexpCaxis.Length()) < 1e-12)
	{
		Vector3d bu,bv,bw;
		bw = axis;
		bw.GenerateComplementBasis(bu,bv,bw);

		v0 = apex + bv * base;
		v1 = apex + axis * height + bv * top;
		double dist;
		bool res = ProjectOntoLineSegment(p,v0,v1,fp,dist);
		signeddist = res?-dist:dist;
		return res;
	}
	else
	{
		bool interior = true;
		double dp = apexp.Dot(axis);
		if( (dp < 0) || (dp > height) ) 
			interior = false;
		else
		{
			double cone_dist_to_axis = dp/height * top + (height-dp)/height * base;
			double dist_to_axis;
			Vector3d tfp;
			ProjectOntoLineSegment(p,apex,apex+axis*height,tfp,dist_to_axis);
			if(dist_to_axis > cone_dist_to_axis)
				interior = false;
		}
		apexp = apexp - apexp.Dot(axis)*axis;
		apexp.Normalize();

		v0 = apex + apexp * base;
		v1 = apex + axis * height + apexp * top;
		double dist;
		bool res = ProjectOntoLineSegment(p,v0,v1,fp,dist);
		if(interior)
			signeddist = -dist;
		else
			signeddist = dist;
		return res;
	}
}

Sphere Cone::BoundingSphere()
{
	Vector3d ep[2];
	ep[0] = apex;
	ep[1] = apex + axis * height;
	double rad;
	rad = 0.5*(height + base + top);
	Sphere s;
	s.center = 0.5 * (ep[0] + ep[1]);
	s.radius = rad;
	return s;
}

std::vector< Sphere > Cone::SampleSpheres( unsigned num)
{
	std::vector< Sphere > vs;

	Sphere s[2];
	s[0].center = apex;
	s[1].center = apex + axis * height;
	s[0].radius = base;
	s[1].radius = top;

	for(unsigned i = 0; i <= num; i ++)
	{
		Sphere ss;
		double t = (double) i / num;
		ss = s[0] * (1.-t)+ s[1] * t;
		vs.push_back(ss);
	}
	return vs;
}

Triangle::Triangle()
{
}

Triangle::~Triangle()
{
}

double Triangle::Area()
{
	return TriangleArea(v[0], v[1], v[2]);
}

Vector3d Triangle::Centroid()
{
	return ( v[0] + v[1] + v[2] ) / 3.;
}

bool Triangle::is_obtuse()
{
	return obtuse_triangle(v[0], v[1], v[2]);
}

double Triangle::anglev(int i)
{
	Vector3d vec0, vec1;
	vec0 = v[(i+2)%3] - v[i];
	vec1 = v[(i+1)%3] - v[i];
	return acos( vec0.Dot(vec1) / vec0.Length() / vec1.Length() );
}

double Triangle::voronoi_area_Meyer(int i)
{
	if(is_obtuse())
	{
		if(anglev(i) > Wm4::Mathd::PI / 2.)
			return Area() / 2.;
		else
			return Area() / 4.;
	}
	else
	{
		double cotanip1, cotania1;
		cotanip1 = cotan(v[(i+2)%3], v[i], v[(i+1)%3]);
		cotania1 = cotan(v[(i+1)%3], v[i], v[(i+2)%3]);
		return ( cotanip1 * (v[i] - v[(i+1)%3]).SquaredLength() + cotania1 * (v[i] - v[(i+2)%3]).SquaredLength() )/ 8.; 
	}
}

double Triangle::voronoi_area_new(int i)
{
	return Area() / 3.;
}


LineSegment::LineSegment()
{
}

LineSegment::~LineSegment()
{
}

double LineSegment::Length()
{
	return (v[0] - v[1]).Length();
}

int clip(PNPlane & pnplane, Triangle t, std::vector<Triangle> & inner, std::vector<Triangle> & outer)
{
	if(t.Area() < 1e-16)
		return -1;
		
	Vector3d p(pnplane.p), n(pnplane.n);
	Vector3d v[3], pv[3];
	double pvdn[3];
	int sign[3], positive_count(0);
	int conf(0);
	
	int pow2i = 1;
	for(int i = 0; i < 3; i ++)
	{
		v[i] = t.v[i];
		pv[i] = v[i] - p;
		pvdn[i] = pv[i].Dot(n);
		if(pvdn[i] < 0.0)
			sign[i] = 0; // inside
		else 
		{
			sign[i] = 1; // outside
			positive_count ++;
		}
		conf += sign[i] * pow2i;
		pow2i *= 2;
	}
	if(positive_count == 3)
	{
		inner.push_back(t); // the whole triangle is in the positive half-space
		return 3;
	}
	else if(positive_count == 2)
	{
		// inner - a quad, i.e., two triangles
		// outer -- a triangle
		Vector3d vo, vi[2];
		if(conf == 3){
			vo = v[2]; vi[0] = v[0]; vi[1] = v[1];}
		else if(conf == 5){
			vo = v[1]; vi[0] = v[2]; vi[1] = v[0];}
		else { // conf == 6
			vo = v[0]; vi[0] = v[1]; vi[1] = v[2];}
		Vector3d intv[2];
		for(int i = 0; i < 2; i ++)
		{
			double deno = (vo-vi[i]).Dot(n);
			if( fabs(deno) == 0.0)
				intv[i] = 0.5 * (vo+vi[i]);
			else
			{
				double u = (vo.Dot(n) - p.Dot(n)) / deno;
				if(u < 0.0) u = 0.0;
				if(u > 1.0) u = 1.0;
				intv[i] = (1. - u) * vo + u * vi[i];
			}
		}
		Triangle to(vo,intv[0], intv[1]), ti0(intv[0], vi[0], vi[1]), ti1(intv[0], vi[1], intv[1]);
		inner.push_back(ti0);
		inner.push_back(ti1);
		outer.push_back(to);
		return 2;
	}
	else if(positive_count == 1)
	{
		// inner - a triangle
		// outer -- a quad, i.e., two triangles
		Vector3d vi, vo[2];
		if(conf == 1){
			vi = v[0]; vo[0] = v[1]; vo[1] = v[2];}
		else if(conf == 2){
			vi = v[1]; vo[0] = v[2]; vo[1] = v[0];}
		else{ // conf == 4
			vi = v[2]; vo[0] = v[0]; vo[1] = v[1];}
		Vector3d intv[2];
		for(int i = 0; i < 2; i ++)
		{
			double deno = (vi-vo[i]).Dot(n);
			if( fabs(deno) == 0.0)
				intv[i] = 0.5 * (vi+vo[i]);
			else
			{
				double u = (vi.Dot(n) - p.Dot(n)) / deno;
				if(u < 0.0) u = 0.0;
				if(u > 1.0) u = 1.0;
				intv[i] = (1. - u) * vi + u * vo[i];
			}
		}
		Triangle ti(vi,intv[0],intv[1]), to0(intv[0], vo[0], vo[1]), to1(intv[0], vo[1], intv[1]);
		inner.push_back(ti);
		outer.push_back(to0);
		outer.push_back(to1);
		return 1;
	}else // positive_count == 0
	{
		outer.push_back(t); // the whole triangle is in the negative half-space
		return 0;
	}
}

double TriangleArea(Vector3d & v0, Vector3d & v1, Vector3d & v2)
{
	return 0.5 * fabs( ( (v1-v0).Cross(v2-v0) ).Length() );
}

Vector3d ThreeAxesScaling(Vector3d & axis0, double & scale0, Vector3d & axis1, double & scale1, Vector3d & axis2, double & scale2, Vector3d & v)
{
	Vector3d result(Vector3d::ZERO);
	result += v.Dot(axis0)*scale0*axis0;
	result += v.Dot(axis1)*scale1*axis1;
	result += v.Dot(axis2)*scale2*axis2;
	return result;
}

// 
//Vector3d TransformToEuclidean(Vector3d & axis0, double & scale0, Vector3d & axis1, double & scale1, Vector3d & axis2, double & scale2, Vector3d & p)
//{
//	return Vector3d(p.Dot(axis0)*scale0, p.Dot(axis1)*scale1, p.Dot(axis2)*scale2);
//}


// project the point p onto the triangle v0-v1-v2, fp is the footpoint, dist is the distance between p and fp. 
// if fp lies in the triangle, return true; otherwise (on the edges or the vertices), return false
bool ProjectOntoLineSegment(const Vector3d & p, const Vector3d & v0, const Vector3d & v1, Vector3d & fp, double & dist)
{
	double t((p-v0).Dot(v1-v0) / (v1-v0).SquaredLength());
	if( (t >= 0.0) && (t <= 1.0) )
	{
		fp = (1.0-t)*v0 + t*v1;
		dist = (p-fp).Length();
		return true;
	}
	else if( t < 0.0)
	{
		fp = v0;
		dist = (p-v0).Length();
		return false;
	}
	else // t > 1.0
	{
		fp = v1;
		dist = (p-v1).Length();
		return false;
	}
}

bool SameSide(const Wm4::Vector3d & p1, const Wm4::Vector3d & p2, const Wm4::Vector3d & a, const Wm4::Vector3d & b)
{
	Vector3d cp1 = (b-a).Cross(p1-a);
	//std::cout << "cp1: " << cp1 << std::endl;
	Vector3d cp2 = (b-a).Cross(p2-a);
	//std::cout << "cp2: " << cp2 << std::endl;
	//std::cout << "cp1.dot(cp2): " << cp1.Dot(cp2) << std::endl;
	if( cp1.Dot(cp2) >= 1e-12)
		return true;
	else
		return false;
}

bool InsideTriangle(const Wm4::Vector3d &p, const Wm4::Vector3d &v0, const Wm4::Vector3d &v1, const Wm4::Vector3d &v2)
{
	if( SameSide(p, v0, v1, v2) && SameSide(p, v1, v0, v2) && SameSide(p, v2, v0, v1) )
		return true;
	else
		return false;
}

// project the point onto the line segment v0-v1, fp is the footpoint, dist is the distance between p and fp.
// if p lies in the segment, return true; otherwise (v0 or v1), return false
bool ProjectOntoTriangle(const Vector3d & p, const Vector3d & v0, const Vector3d & v1, const Vector3d & v2, Vector3d & fp, double & dist)
{
	Vector3d norm( (v0-v1).Cross(v0-v2) );
	norm.Normalize();
	fp = p - (p-v1).Dot(norm)*norm;
	if( InsideTriangle(fp, v0, v1, v2) )
	{
		dist = (fp-p).Length();
		return true;
	}
	else
	{
		Vector3d fp01, fp02, fp12;
		double dist01, dist02, dist12;
		ProjectOntoLineSegment(p,v0,v1,fp01,dist01);
		ProjectOntoLineSegment(p,v0,v2,fp02,dist02);
		ProjectOntoLineSegment(p,v1,v2,fp12,dist12);
		if( (dist01 <= dist02) && (dist01 <= dist12) )
		{
			dist = dist01;
			fp = fp01;
			return false;
		}
		else if( (dist02 <= dist01) && (dist02 <= dist12) )
		{
			dist = dist02;
			fp = fp02;
			return false;
		}
		else
		{
			dist = dist12;
			fp = fp12;
			return false;
		}
	}
}

bool obtuse_triangle(Vector3d & v0, Vector3d & v1, Vector3d & v2)
{
	Vector3d v0v1, v0v2, v1v2;
	v0v1 = v1 - v0;
	v0v2 = v2 - v0;
	v1v2 = v2 - v1;
	
	double slv0v1, slv0v2, slv1v2;
	slv0v1 = v0v1.SquaredLength();
	slv0v2 = v0v2.SquaredLength();
	slv1v2 = v1v2.SquaredLength();
	
	if(slv0v1 + slv0v2 < slv1v2)
		return true;
	
	if(slv0v1 + slv1v2 < slv0v2)
		return true;
	
	if(slv0v2 + slv1v2 < slv0v1)
		return true;	
	
	return false;
}

bool acute_triangle(Vector3d & v0, Vector3d & v1, Vector3d & v2)
{
	return !obtuse_triangle(v0,v1,v2);
}

double cotan(Vector3d v, Vector3d v1, Vector3d v2)
{
	Vector3d vv1, vv2;
	vv1 = v1 - v;
	vv2 = v2 - v;
	
	double numerator;
	double denominator;
	numerator = vv1.Dot(vv2);
	denominator = sqrt(vv1.SquaredLength() * vv2.SquaredLength() - numerator*numerator);
	if(denominator == 0.0) return (0.0);

	return (numerator / denominator);
}

double angle_from_cotan(Vector3d v, Vector3d v0, Vector3d v1)
{
	Vector3d vv0, vv1;
	vv0 = v0 - v;
	vv1 = v1 - v;
	
	double numerator(vv0.Dot(vv1));
	double denominator(sqrt(vv0.SquaredLength() * vv1.SquaredLength() - numerator * numerator));
	
	return ( fabs( atan2(denominator, numerator) ) );
}

// project the point p onto the ellipsoid x^2 / a^2 + y^2 / b^2 + z^2 / c^2 = 1
Wm4::Vector3d ProjectPointOntoEllipsoid(const double a, const double b, const double c, const Wm4::Vector3d & p)
{

	// find the intersection fp of the ellipsoid and the ray from the origin to p
	Wm4::Vector3d fp = p;

	double sum(0.0);
	sum = fp[0] * fp[0] / a / a + fp[1] * fp[1] / b / b + fp[2] * fp[2] / c / c;
	sum = sqrt(sum);

	if(sum > 1e-12)
		fp /= sum;
	else
	{
		// p is the origin
		if((a <= b) && (a <= c))
		{
			return Wm4::Vector3d(a,0,0);
		}
		else if((b <= a) && (b <= c))
		{
			return Wm4::Vector3d(0,b,0);
		}
		else
		{
			return Wm4::Vector3d(0,0,c);
		}
	}


	// find the parameters of the intersecting point fp
	double u,v;
	double su = fp[2]/c;
	u = asin(su);

	double cu = cos(u);
	double cv,sv;
	cv = fp[0]/a/cu;
	sv = fp[1]/b/cu;

	v = acos(cv);
	if(sv < 0)
		v = -v;

	// move fp to reduce the distance 
	double current_sqdist = (fp-p).SquaredLength();
	
	double step = 1e-2;

	bool goodstep;
	do
	{
		bool decreased;
		do
		{
			decreased = false;

			double nu[4], nv[4];
			nu[0] = u+step; nv[0] = v;
			nu[1] = u-step; nv[1] = v;
			nu[2] = u; nv[2] = v+step;
			nu[3] = u; nv[3] = v-step;

			double minsqdist(1e20);
			unsigned min_idx(0);
			for(unsigned k = 0; k < 4; k ++)
			{
				Wm4::Vector3d nup;
				nup[0] = a * cos(nu[k]) * cos(nv[k]);
				nup[1] = b * cos(nu[k]) * sin(nv[k]);
				nup[2] = c * sin(nu[k]);
				double tempsqdist;
				tempsqdist = (nup-p).SquaredLength();
				if(tempsqdist < minsqdist)
				{
					minsqdist = tempsqdist;
					min_idx = k;
				}
			}

			if(minsqdist < current_sqdist)
			{
				u = nu[min_idx];
				v = nv[min_idx];
				current_sqdist = minsqdist;
				decreased = true;
			}

		}while(decreased);

		step /= 2.;
		if(step > 1e-4)
			goodstep = true;
		else 
			goodstep = false;
	}while(goodstep);



	fp[0] = a * cos(u) * cos(v);
	fp[1] = b * cos(u) * sin(v);
	fp[2] = c * sin(u);

	return fp;
}

// project the point p onto the elliptic cone x^2 / a^2 + z^2 / c^2 = y^2 / b^2
Wm4::Vector3d ProjectPointOntoEllipticCone(const double a, const double b, const double c, const Wm4::Vector3d & p)
{
	// find the intersection fp of the ellipsoid and the ray from the origin to p
	Wm4::Vector3d fp = p;

	double sum(0.0);
	sum = fp[0] * fp[0] / a / a - fp[1] * fp[1] / b / b + fp[2] * fp[2] / c / c;
	sum = sqrt(sum);

	if(sum > 1e-12)
		fp /= sum;
	else
	{
		return fp;
	}


	// find the parameters of the intersecting point fp
	double u,v;

	u = fp[1] / b;


	double cv,sv;
	cv = fp[0]/a/u;
	sv = fp[2]/c/u;

	v = acos(cv);
	if(sv < 0)
		v = -v;

	// move fp to reduce the distance 
	double current_sqdist = (fp-p).SquaredLength();
	
	double step = 1e-2;

	bool goodstep;
	do
	{
		bool decreased;
		do
		{
			decreased = false;

			double nu[4], nv[4];
			nu[0] = u+step; nv[0] = v;
			nu[1] = u-step; nv[1] = v;
			nu[2] = u; nv[2] = v+step;
			nu[3] = u; nv[3] = v-step;

			double minsqdist(1e20);
			unsigned min_idx(0);
			for(unsigned k = 0; k < 4; k ++)
			{
				Wm4::Vector3d nup;
				/*
				nup[0] = a * cos(nu[k]) * cos(nv[k]);
				nup[1] = b * cos(nu[k]) * sin(nv[k]);
				nup[2] = c * sin(nu[k]);
				*/
				nup[0] = nu[k]*a*cos(nv[k]);
				nup[1] = nu[k]*b;
				nup[2] = nu[k]*c*sin(nv[k]);

				double tempsqdist;
				tempsqdist = (nup-p).SquaredLength();
				if(tempsqdist < minsqdist)
				{
					minsqdist = tempsqdist;
					min_idx = k;
				}
			}

			if(minsqdist < current_sqdist)
			{
				u = nu[min_idx];
				v = nv[min_idx];
				current_sqdist = minsqdist;
				decreased = true;
			}

		}while(decreased);

		step /= 2.;
		if(step > 1e-4)
			goodstep = true;
		else 
			goodstep = false;
	}while(goodstep);



	fp[0] = u*a*cos(v);
	fp[1] = u*b;
	fp[2] = u*c*sin(v);

	return fp;
}


/*
void ProjectPointOntoEllipsoid(double a, double b, double c, Wm4::Vector3d & p)
{
	double sum(0.0);
	sum = p[0] * p[0] / a / a + p[1] * p[1] / b / b + p[2] * p[2] / c / c;
	sum = sqrt(sum);

	if(sum != 0.0)
		p /= sum;
}
*/
bool DifferentialInfoOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p, 
	Wm4::Vector3d & pu, Wm4::Vector3d & pv, Wm4::Vector3d & normdir, 
	double & maxc, double & minc, Wm4::Vector3d & maxdir, Wm4::Vector3d & mindir)
{
	// the point does not lie on the ellipsoid
	if(fabs(p[0]*p[0]/a/a+p[1]*p[1]/b/b+p[2]*p[2]/c/c-1.) > 1e-12)
		return false;

	// x = a cosu sinv
	// y = b sinu sinv
	// z = c cosv

	//double v;
	//double u;

	double cosv = p[2] / c;
	double sinv = sqrt(1. - cosv*cosv);
	double cosu = p[0] / a / sinv;
	double sinu = p[1] / b / sinv;


	//v = acos( p[2] / c);
	//u = atan( p[1] * a / p[0] / b);
	
	pu[0] = -a * sinu * sinv;
	pu[1] = b * cosu * sinv;
	pu[2] = 0.0;

	pv[0] = a * cosu * cosv;
	pv[1] = b * sinu * cosv;
	pv[2] = -c * sinv;

	normdir = pu.Cross(pv);
	normdir.Normalize();

	Wm4::Vector3d puu, puv, pvv;

	puu[0] = -a * cosu * sinv;
	puu[1] = -b * sinu * sinv;
	puu[2] = 0.0;

	puv[0] = -a * sinu * cosv;
	puv[1] = b * cosu * cosv;
	puv[2] = 0.0;

	pvv[0] = -a * cosu * sinv;
	pvv[1] = -b * sinu * sinv;
	pvv[2] = -c * cosv;

	Wm4::Matrix2d fff;
	fff[0][0] = pu.Dot(pu);
	fff[0][1] = fff[1][0] = pu.Dot(pv);
	fff[1][1] = pv.Dot(pv);

	Wm4::Matrix2d sff;
	sff[0][0] = puu.Dot(normdir);
	sff[0][1] = sff[1][0] = puv.Dot(normdir);
	sff[1][1] = pvv.Dot(normdir);

	Wm4::Matrix2d so;

	so[0][0] = sff[0][0]*fff[1][1] - sff[0][1]*fff[0][1];
	so[0][1] = sff[0][1]*fff[1][1] - sff[1][1]*fff[0][1];
	so[1][0] = sff[0][1]*fff[0][0] - sff[0][0]*fff[0][1];
	so[1][1] = sff[1][1]*fff[0][0] - sff[0][1]*fff[0][1];

	so *= 1. / fff.Determinant();

	double gc, mc;
	gc = so.Determinant();
	mc = 0.5 * (so[0][0] + so[1][1]);

	Wm4::Matrix2d rotmat, diagmat;
	so.EigenDecomposition(rotmat, diagmat);

	if(diagmat[0][0] > diagmat[1][1])
	{
		maxc = diagmat[0][0];
		minc = diagmat[1][1];
		maxdir = rotmat[0][0] * pu+ rotmat[0][1] * pv;
		mindir = rotmat[1][0] * pu+ rotmat[1][1] * pv;
	}
	else
	{
		maxc = diagmat[1][1];
		minc = diagmat[0][0];		
		maxdir = rotmat[1][0] * pu+ rotmat[1][1] * pv;
		mindir = rotmat[0][0] * pu+ rotmat[0][1] * pv;
	}

	return true;
}

double GaussianCurvatureOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p)
{
	Wm4::Vector3d pu,pv,normdir,maxdir,mindir;
	double maxc, minc;
	DifferentialInfoOnEllipsoid(a, b, c, p, pu, pv, normdir, maxc, minc, maxdir, mindir);
	return maxc * minc;
}

double MeanCurvatureOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p)
{
	Wm4::Vector3d pu,pv,normdir,maxdir,mindir;
	double maxc, minc;
	DifferentialInfoOnEllipsoid(a, b, c, p, pu, pv, normdir, maxc, minc, maxdir, mindir);
	return 0.5*(maxc + minc);
}

double MaxCurvatureOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p)
{
	Wm4::Vector3d pu,pv,normdir,maxdir,mindir;
	double maxc, minc;
	DifferentialInfoOnEllipsoid(a, b, c, p, pu, pv, normdir, maxc, minc, maxdir, mindir);
	return maxc;
}

double MinCurvatureOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p)
{
	Wm4::Vector3d pu,pv,normdir,maxdir,mindir;
	double maxc, minc;
	DifferentialInfoOnEllipsoid(a, b, c, p, pu, pv, normdir, maxc, minc, maxdir, mindir);
	return minc;
}

Wm4::Vector3d NormalDirectionOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p)
{
	Wm4::Vector3d pu,pv,normdir,maxdir,mindir;
	double maxc, minc;
	DifferentialInfoOnEllipsoid(a, b, c, p, pu, pv, normdir, maxc, minc, maxdir, mindir);
	return normdir;
}

Wm4::Vector3d MaxPrincipalDirectionOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p)
{
	Wm4::Vector3d pu,pv,normdir,maxdir,mindir;
	double maxc, minc;
	DifferentialInfoOnEllipsoid(a, b, c, p, pu, pv, normdir, maxc, minc, maxdir, mindir);
	return maxdir;
}

Wm4::Vector3d MinPrincipalDirectionOnEllipsoid(double a, double b, double c, Wm4::Vector3d & p)
{
	Wm4::Vector3d pu,pv,normdir,maxdir,mindir;
	double maxc, minc;
	DifferentialInfoOnEllipsoid(a, b, c, p, pu, pv, normdir, maxc, minc, maxdir, mindir);
	return mindir;
}

Wm4::Vector3d LinearInterpolation(Wm4::Vector3d & originalpoint, Wm4::Vector3d & offsetpoint, double offsetfactor)
{
	return offsetfactor * originalpoint + (1. - offsetfactor) * offsetpoint;
}

Wm4::Vector3d RandomPointonSphere()
{
	Wm4::Vector3d p;

	srand((unsigned int)time(NULL));


	double angle = Wm4::Mathd::IntervalRandom(0., Wm4::Mathd::TWO_PI);
	double cs = cos(angle);
	double sn = sin(angle);

	double value = Wm4::Mathd::SymmetricRandom();
	double complement = Wm4::Mathd::Sqrt(Wm4::Mathd::FAbs(1. - value * value) );
	p[0] = value;
	p[1] = complement * cs;
	p[2] = complement * sn;

	return p;
}

bool FastNoIntersect(const Wm4::Vector3d & p, const Wm4::Vector3d & dir, const Wm4::Vector3d & q, const Wm4::Vector3d & norm)
{
	Wm4::Vector3d qp;
	qp = p - q;

	double qp_d_norm = qp.Dot(norm);
	double dir_d_norm = dir.Dot(norm);

	if( ( qp_d_norm > 0) && (dir_d_norm >= 0) )
		return true;
	if( (qp_d_norm < 0 ) && (dir_d_norm <= 0) )
		return true;

	return false;
}

bool RayTriangleIntersect(const Wm4::Vector3d & p, const Wm4::Vector3d & dir, const Wm4::Vector3d & v0, const Wm4::Vector3d & v1, const Wm4::Vector3d & v2)
{
	Wm4::Vector3d norm;
	norm = (v1 - v0).Cross(v2 - v0);
	norm.Normalize();

	if(abs(dir.Dot(norm)) < 1e-18)
		return false;

	double intt;
	intt = (v0-p).Dot(norm) / dir.Dot(norm);

	Wm4::Vector3d intp;
	intp = p + intt * dir;

	if(intt <= 0)
		return false;

	if(InsideTriangle(intp, v0, v1, v2))
		return true;
	else
		return false;	
}

bool RayTriangleIntersectv2(const Wm4::Vector3d & p, const Wm4::Vector3d & dir, const Wm4::Vector3d & v0, const Wm4::Vector3d & v1, const Wm4::Vector3d & v2)
{
	Wm4::Vector3d diff = p - v0;
	Wm4::Vector3d edge1 = v1 - v0;
	Wm4::Vector3d edge2 = v2 - v0;
	Wm4::Vector3d norm = edge1.Cross(edge2);

	double DdN = dir.Dot(norm);
	double sign;

	if(DdN > 0)
		sign = 1;
	else if(DdN < -0)
			sign = -1;
	else
		return false;

	double DdQxE2 = sign * dir.Dot(diff.Cross(edge2));

	if(DdQxE2 > 0)
	{
		double DdE1xQ = sign * dir.Dot(edge1.Cross(diff));
		if(DdE1xQ > 0)
		{
			if( DdQxE2 + DdE1xQ <= DdN)
			{
				double QdN = -sign * diff.Dot(norm);
				if(QdN > 0)
				{
					return true;
				}
			}
		}
	}

	return false;
}

bool DistanceToLine(const Vector3d & p, const Vector3d & v0, const Vector3d & v1, double & dist, Vector3d & fp)
{
	Vector3d v0v1(v1-v0), pv0(v0-p), pv1(v1-p);
	double area = fabs(v0v1.Cross(pv0).Length());
	if(v0v1.Length() > 1e-12)
	{
		dist = area / v0v1.Length();
		double t = (pv0.Dot(pv0)-pv0.Dot(pv1)) / (pv0.Dot(pv0)+pv1.Dot(pv1)-2*pv0.Dot(pv1));
		fp = (1-t)*v0+t*v1;
		return true;
	}
	else
		return false;
}

bool TriangleFromThreeSpheres(const Vector3d & c0,
							  const double & r0,
							  const Vector3d & c1,
							  const double & r1,
							  const Vector3d & c2,
							  const double & r2,
							  SimpleTriangle & st0,
							  SimpleTriangle & st1)
{
	//Vector3d c0c1(c1-c0),c0c2(c2-c0),c1c2(c2-c1);
	//double c0c1len(c0c1.Length()),c0c2len(c0c2.Length()),c1c2len(c1c2.Length());
	//double dr0r1(fabs(r0-r1)),dr0r2(fabs(r0-r2)),dr1r2(fabs(r1-r2));

	//// some spheres are concentric and there are no triangles.
	//if( (c0c1len < 1e-8) || (c0c2len < 1e-8) || (c1c2len < 1e-8) )
	//	return false;

	//Vector3d norm;
	//norm = c0c1.Cross(c0c2);
	//norm.Normalize();

	//// equal-radius spheres
	//if( (dr0r1 < 1e-8) && (dr0r2 < 1e-8) && (dr1r2 < 1e-8) )
	//{
	//	st0.v[0] = c0 + norm * r0;
	//	st0.v[1] = c1 + norm * r1;
	//	st0.v[2] = c2 + norm * r2;
	//	st0.UpdateNormal();

	//	st1.v[0] = c0 - norm * r0;
	//	st1.v[1] = c1 - norm * r1;
	//	st1.v[2] = c2 - norm * r2;
	//	st1.UpdateNormal();
	//	return true;
	//}
	//else 
	//{
	//	// two points on the tangent plane
	//	Vector3d apex0,apex1;

	//	if ((r0 < 1e-8 && r1 < 1e-8) || (r0 < 1e-8 && r2 < 1e-8) || (r1 < 1e-8 && r2 < 1e-8))
	//	{
	//		double distc0;
	//		Vector3d fp;
	//		double sangle;
	//		double cangle;
	//		Vector3d norfpc0;
	//		// two spheres are equal-radius
	//		if (r0 > 1e-8)
	//		{
	//			apex0 = (r0 * c2 - r2 * c0) / (r0 - r2);
	//			apex1 = (r0 * c1 - r1 * c0) / (r0 - r1);
	//			apex0 = (r1 * c0 - r0 * c1) / (r1 - r0);
	//			apex1 = (r1 * c2 - r2 * c1) / (r1 - r2);
	//			DistanceToLine(c0,apex0,apex1,distc0,fp);
	//			sangle = r0/distc0;
	//			if(fabs(sangle) > 1.)
	//				return false;
	//			cangle = sqrt(1.-r0*r0/distc0/distc0);
	//			norfpc0 = (c0-fp);
	//			norfpc0.Normalize();
	//		}
	//		else if (r1 > 1e-8)
	//		{
	//			apex0 = (r1 * c0 - r0 * c1) / (r1 - r0);
	//			apex1 = (r1 * c2 - r2 * c1) / (r1 - r2);
	//			DistanceToLine(c1,apex0,apex1,distc0,fp);
	//			sangle = r1/distc0;
	//			if(fabs(sangle) > 1.)
	//				return false;
	//			cangle = sqrt(1.-r1*r1/distc0/distc0);
	//			norfpc0 = (c1-fp);
	//			norfpc0.Normalize();
	//		}
	//		else if (r2 > 1e-8)
	//		{
	//			apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
	//			apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
	//			DistanceToLine(c2,apex0,apex1,distc0,fp);
	//			sangle = r2/distc0;
	//			if(fabs(sangle) > 1.)
	//				return false;
	//			cangle = sqrt(1.-r2*r2/distc0/distc0);
	//			norfpc0 = (c2-fp);
	//			norfpc0.Normalize();
	//		}

	//		Vector3d newnorm[2];
	//		newnorm[0] = norm*cangle - norfpc0*sangle;
	//		newnorm[1] = -norm*cangle - norfpc0*sangle;
	//		st0.v[0] = c0 + r0*newnorm[0];
	//		st0.v[1] = c1 + r1*newnorm[0];
	//		st0.v[2] = c2 + r2*newnorm[0];
	//		st0.UpdateNormal();

	//		st1.v[0] = c0 + r0*newnorm[1];
	//		st1.v[1] = c1 + r1*newnorm[1];
	//		st1.v[2] = c2 + r2*newnorm[1];
	//		st1.UpdateNormal();
	//	}
	//	else
	//	{
	//		// two spheres are equal-radius
	//		if(dr0r1 < 1e-8)
	//		{
	//			apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
	//			apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
	//		}
	//		else if (dr0r2 < 1e-8)
	//		{
	//			apex0 = (r1 * c0 - r0 * c1) / (r1 - r0);
	//			apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
	//		}
	//		else if (dr1r2 < 1e-8)
	//		{
	//			apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
	//			apex1 = (r0 * c1 - r1 * c0) / (r0 - r1);
	//		}
	//		else
	//		{
	//			apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
	//			apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
	//		}

	//		double distc0;
	//		Vector3d fp;
	//		DistanceToLine(c0,apex0,apex1,distc0,fp);
	//		double sangle = r0/distc0;
	//		if(fabs(sangle) > 1.)
	//			return false;
	//		double cangle = sqrt(1.-r0*r0/distc0/distc0);
	//		Vector3d norfpc0(c0-fp);
	//		norfpc0.Normalize();
	//		Vector3d newnorm[2];
	//		newnorm[0] = norm*cangle - norfpc0*sangle;
	//		newnorm[1] = -norm*cangle - norfpc0*sangle;
	//		st0.v[0] = c0 + r0*newnorm[0];
	//		st0.v[1] = c1 + r1*newnorm[0];
	//		st0.v[2] = c2 + r2*newnorm[0];
	//		st0.UpdateNormal();

	//		st1.v[0] = c0 + r0*newnorm[1];
	//		st1.v[1] = c1 + r1*newnorm[1];
	//		st1.v[2] = c2 + r2*newnorm[1];
	//		st1.UpdateNormal();
	//	}
	//}

	//return true;

	Vector3d c0c1(c1-c0),c0c2(c2-c0),c1c2(c2-c1);
	double c0c1len(c0c1.Length()),c0c2len(c0c2.Length()),c1c2len(c1c2.Length());
	double dr0r1(fabs(r0-r1)),dr0r2(fabs(r0-r2)),dr1r2(fabs(r1-r2));

	// some spheres are concentric and there are no triangles.
	if( (c0c1len < 1e-8) || (c0c2len < 1e-8) || (c1c2len < 1e-8) )
		return false;

	//// some spheres are included in some other spheres 
	//if ((c0c1len - abs(r0 - r1) < 1e-8) || (c0c2len - abs(r0 - r2) < 1e-8) || (c1c2len - abs(r1 - r2) < 1e-8))
	//	return false;

	Vector3d norm;
	norm = c0c1.Cross(c0c2);
	norm.Normalize();

	// equal-radius spheres
	if( (dr0r1 < 1e-8) && (dr0r2 < 1e-8) && (dr1r2 < 1e-8) )
	{
		st0.v[0] = c0 + norm * r0;
		st0.v[1] = c1 + norm * r1;
		st0.v[2] = c2 + norm * r2;
		st0.UpdateNormal();

		st1.v[0] = c0 - norm * r0;
		st1.v[1] = c1 - norm * r1;
		st1.v[2] = c2 - norm * r2;
		st1.UpdateNormal(true);
		return true;
	}
	else 
	{
		// two points on the tangent plane
		Vector3d apex0,apex1;

		// two spheres are equal-radius
		if(dr0r1 < 1e-8)
		{
			apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
			apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
		}
		else if (dr0r2 < 1e-8)
		{
			apex0 = (r1 * c0 - r0 * c1) / (r1 - r0);
			apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
		}
		else if (dr1r2 < 1e-8)
		{
			apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
			apex1 = (r0 * c1 - r1 * c0) / (r0 - r1);
		}
		else
		{
			apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
			apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
		}

		double distc0;
		Vector3d fp;
		DistanceToLine(c0,apex0,apex1,distc0,fp);
		double sangle = r0/distc0;
		if(fabs(sangle) > 1.)
			return false;
		double cangle = sqrt(1.-r0*r0/distc0/distc0);
		Vector3d norfpc0(c0-fp);
		norfpc0.Normalize();
		Vector3d newnorm[2];
		newnorm[0] = norm*cangle - norfpc0*sangle;
		newnorm[1] = -norm*cangle - norfpc0*sangle;
		st0.v[0] = c0 + r0*newnorm[0];
		st0.v[1] = c1 + r1*newnorm[0];
		st0.v[2] = c2 + r2*newnorm[0];
		st0.UpdateNormal();

		st1.v[0] = c0 + r0*newnorm[1];
		st1.v[1] = c1 + r1*newnorm[1];
		st1.v[2] = c2 + r2*newnorm[1];
		st1.UpdateNormal(true);
	}

	return true;
}

Vector3d TriangleNormal(const Vector3d & v0, const Vector3d & v1, const Vector3d & v2)
{
	Vector3d v01, v02;
	v01 = v1 - v0;
	v02 = v2 - v0;
	Vector3d norm;
	norm = v01.Cross(v02);
	norm.Normalize();
	return norm;
}


double VectorAngle(const Vector3d & v0, const Vector3d & v1)
{
	double dot_multi_result = v0.Dot(v1) / v0.Length() / v1.Length();
	
	return acos(dot_multi_result);
}