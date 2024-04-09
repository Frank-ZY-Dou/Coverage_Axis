#ifndef _MESH_H
#define _MESH_H

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/convex_hull_2.h>


#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include <CGAL/IO/Polyhedron_iostream.h>



#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

#include <CGAL/Cartesian_d.h>
#include <CGAL/Min_sphere_d.h>
#include <CGAL/Min_sphere_annulus_d_traits_d.h>

#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <fstream>
#include <string>
#include <cmath>


#include "LinearAlgebra/Wm4Vector.h"
#include "LinearAlgebra/Wm4Matrix.h"
#include "GeometryObjects/GeometryObjects.h"


typedef double simple_numbertype;
typedef CGAL::Cartesian<simple_numbertype> simple_kernel;
typedef CGAL::Polyhedron_3<simple_kernel, CGAL::Polyhedron_items_3> SimpleMesh;

typedef CGAL::Point_3< CGAL::Cartesian<simple_numbertype> > Point;

using namespace std;
using namespace Wm4;

# define MY_PI 3.14159265358979323846


class VertexInfo
{
public:
	int id;
	int tag;
	std::vector<int> tagvec;
	std::set<unsigned int> pole_bplist;
	double var_double;
};

class CellInfo
{
public:
	bool inside;
	int id;
	int tag;
	bool is_pole;
	std::set<unsigned int> pole_bplist;
	double dist_center_to_boundary; // approximate 
};


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K> Vb;
typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo, K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
//typedef CGAL::Delaunay_triangulation_3<K, Tds> Triangulation;


class Triangulation : public CGAL::Delaunay_triangulation_3<K, Tds>
{
public:
	double TetCircumRadius(const Tetrahedron & tet);
	double TetLargestEdgeLength(const Tetrahedron & tet);
	double TetBoundingSphereRadius(const Tetrahedron & tet);
};



typedef Triangulation::Vertex_handle Vertex_handle_t;
typedef Triangulation::Cell_handle Cell_handle_t;

typedef Triangulation::Vertex_iterator Vertex_iterator_t;
typedef Triangulation::Edge_iterator Edge_iterator_t;
typedef Triangulation::Facet_iterator Facet_iterator_t;
typedef Triangulation::Cell_iterator Cell_iterator_t;

typedef Triangulation::Finite_cells_iterator Finite_cells_iterator_t;
typedef Triangulation::Finite_facets_iterator Finite_facets_iterator_t;
typedef Triangulation::Finite_edges_iterator Finite_edges_iterator_t;
typedef Triangulation::Finite_vertices_iterator Finite_vertices_iterator_t;

typedef Triangulation::Facet_circulator Facet_circulator_t;
typedef Triangulation::Cell_circulator Cell_circulator_t;

typedef Triangulation::Point Point_t;


typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;



typedef CGAL::Cartesian_d<double> CarK;
typedef CGAL::Min_sphere_annulus_d_traits_d<CarK> MSTraits;
typedef CGAL::Min_sphere_d<MSTraits> Min_sphere;

typedef CarK::Point_d MSPoint;



// conversion between CGAL::Point and Wm4::Vector3d

inline Wm4::Vector3d to_wm4(const Point & p)
{
	return Wm4::Vector3d( p.x(), p.y(), p.z() );
}

inline Wm4::Vector3d to_wm4(const CGAL::Vector_3< CGAL::Cartesian<simple_numbertype> > & p)
{
	return Wm4::Vector3d( p.x(), p.y(), p.z() );
}

inline Wm4::Vector3d to_wm4(const Point_t & p)
{
	return Wm4::Vector3d( p.x(), p.y(), p.z() );
}

inline Point to_cgal(const Wm4::Vector3d & p) 
{
	return Point( p[0], p[1], p[2] );
}



enum DENSITYPOLICY
{
	ONE,
	SQRTK1K2,
	K1K2	
};


enum METRICPOLICY
{
	K1_K2,
	SQUAREK1_SQUAREK2
};

enum VertexDiffType
{
	ELLIPTIC,
	HYPERBOLIC,
	PARABOLIC,
	PLANAR,
	BORDER,
	UNDEFINED,
	BADDDT
};

template <class Refs, class T, class P, class Norm>
class MPFacet: public CGAL::HalfedgeDS_face_base<Refs, T>
{
public:
	Vector3d color;
	int id;
	Norm normal;
	typedef Norm Normal;
	int tag;

	Matrix3d metric;
	double density;

	// curvature
	Vector3d maxdir;
	Vector3d mindir;
	Vector3d normdir;
	double maxcurvature;
	double mincurvature;
	// curvature 

public:
	MPFacet(){ tag = 0; }
};

template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class MPHalfedge : public CGAL::HalfedgeDS_halfedge_base<Refs, Tprev, Tvertex, Tface>
{
public:
	int id;
	int tag;
public:
	MPHalfedge(){ tag = 0; }
};

template <class Refs, class T, class P, class N, class Kernel>
class MPVertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
public:
	Vector3d color;
	int id;
	int tag;
	typedef N Normal;
	Normal normal;
	Matrix3d metric;
	double density;

	double vqem_hausdorff_dist;
	unsigned vqem_hansdorff_index;

	double slab_hausdorff_dist;
	unsigned slab_hansdorff_index;

	// the matrix of A and b
	Wm4::Matrix4d A;
	Wm4::Vector4d b;
	double c;

	// curvature
	Vector3d maxdir;
	Vector3d mindir;
	Vector3d normdir;
	double maxcurvature;
	double mincurvature;
	// curvature 
	
	// a scalar value defined on a vertex
	double functionvalue;
	// a scalar value defined on a vertex

	// estimated differential information from [Meyer et al 2002]
	Vector3d normal_Meyer;
	Vector3d maxdir_Meyer, mindir_Meyer;
	double maxcurvature_Meyer, mincurvature_Meyer;
	double meancurvature_Meyer;
	double gaussiancurvature_Meyer;
	double gaussiancurvature_new;
	// estimated differential information from [Meyer et al 2002]


	// differential type of the vertex
	VertexDiffType cdt;
	VertexDiffType ddt;
	// differential type of the vertex
	
	
public:
	MPVertex(){ tag = 0; functionvalue = 0.0; cdt = UNDEFINED; ddt = UNDEFINED; c = 0;}
	MPVertex(const P & pt) : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt) {c = 0;}	
};

struct MPItems : public CGAL::Polyhedron_items_3
{
	template<class Refs, class Traits>
	struct Vertex_wrapper
	{
		typedef typename Traits::Point_3 Point;
		typedef typename Traits::Vector_3 Normal;
		typedef MPVertex<Refs, CGAL::Tag_true, Point, Normal, Traits> Vertex;
	};

	template <class Refs, class Traits>
	struct Face_wrapper
	{
		typedef typename Traits::Point_3 Point;
		typedef typename Traits::Vector_3 Normal;
		typedef MPFacet<Refs, CGAL::Tag_true, Point, Normal> Face;
	};

	template <class Refs, class Traits>
	struct Halfedge_wrapper
	{
		typedef typename Traits::Vector_3 Normal;
		typedef MPHalfedge<Refs, CGAL::Tag_true, CGAL::Tag_true, CGAL::Tag_true, Normal> Halfedge;
	};
};

class MPMesh : public CGAL::Polyhedron_3<simple_kernel, MPItems>
{
public:
	typedef simple_kernel::FT FT;
	typedef simple_kernel::Point_3 Point;
	typedef simple_kernel::Vector_3 Vector;

public:
	simple_kernel::FT m_min[3];
	simple_kernel::FT m_max[3];
	unsigned int m_nb_components;
	unsigned int m_nb_boundaries;
	int m_genus;

	std::vector<Vertex_iterator> pVertexList;
	std::vector<Facet_iterator> pFaceList;

public:
	DENSITYPOLICY m_density_policy;
	METRICPOLICY m_metric_policy;
public:
	MPMesh();
	~MPMesh()
	{
		if(domain != NULL)
			delete domain;
	}

	void computebb(); // implemented
	void GenerateVertexList(); // implemented
	void GenerateFaceList(); // implemented
	void GenerateList(); // implemented
	void GenerateRandomColor(); // implemented

	

	void compute_normals_per_facet(); // implemented
	void compute_normals_per_vertex(); // implemented
	void compute_normals(); // implemented
	void copybb(MPMesh * pmesh); // implemented

	// compute the matrix of A and b for sphere mesh
	void compute_sphere_matrix();

	static unsigned int degree(Facet_handle pFace)
	{
		return (unsigned int)(CGAL::circulator_size(pFace->facet_begin()));
	}
	static unsigned int valence(Vertex_handle pVertex)
	{
		return (unsigned int)(CGAL::circulator_size(pVertex->vertex_begin()));
	}
	static bool is_border(Vertex_handle pVertex)
	{
		Halfedge_around_vertex_circulator pHalfedge = pVertex->vertex_begin();
		if(pHalfedge == NULL) // isolated vertex
			return true;
		Halfedge_around_vertex_circulator begin = pHalfedge;
		CGAL_For_all(pHalfedge, begin)
			if(pHalfedge->is_border())
				return true;
		return false;
	}
	void tag_facets(const int tag); // implemented
	void tag_halfedges(const int tag); // implemented
	void tag_vertices(const int tag); // implemented
	void compute_facet_center(Facet_handle pFace, Point & center); // implemented
	Halfedge_handle get_border_halfedge_tag(int tag); // implemented
	Facet_handle get_facet_tag(const int tag); // implemented
	void tag_component(Facet_handle pSeedFacet, const int tag_free, const int tag_done); // implemented
	unsigned int nb_boundaries(); // implemented
	unsigned int nb_components(); // implemented
	bool is_simple_watertight();
	int genus(); // implemented
	int genus(int c, int v, int f, int e, int b); // implemented
	void compute_components_boundaries_genus(); // implemented

	unsigned int get_nb_components(); // implemented
	unsigned int get_nb_boundaries(); // implemented
	int get_genus(); // implemented
	Wm4::Vector3d GetCentroid(Face_iterator pFace); // implemented
	double GetFaceLargestAngle(Facet_handle pFace); // implemented
	double GetFaceSmallestAngle(Facet_handle pFace); // implemented
	double GetArea(Facet_handle pFace);


	// check point inside the polyhedron
	bool inside(const Wm4::Vector3d & p);


	// check whether the point is inside the bounding box
	bool inside_boundingbox(const Wm4::Vector3d & p);

	// curvature
	double m_max_v_gaussiancurvature;
	double m_min_v_gaussiancurvature;
	double m_max_v_meancurvature;
	double m_min_v_meancurvature;
	double m_max_v_abs_gaussiancurvature;
	double m_min_v_abs_gaussiancurvature;
	double m_max_v_ratio_principalcurvatures;
	double m_min_v_ratio_principalcurvatures;
	double m_max_v_maxcurvature;
	double m_min_v_maxcurvature;
	double m_max_v_mincurvature;
	double m_min_v_mincurvature;

	double m_max_f_gaussiancurvature;
	double m_min_f_gaussiancurvature;
	double m_max_f_meancurvature;
	double m_min_f_meancurvature;
	double m_max_f_abs_gaussiancurvature;
	double m_min_f_abs_gaussiancurvature;
	double m_max_f_ratio_principalcurvatures;
	double m_min_f_ratio_principalcurvatures;
	double m_max_f_maxcurvature;
	double m_min_f_maxcurvature;
	double m_max_f_mincurvature;
	double m_min_f_mincurvature;

	void EvaluateCurvatureInterval();
	void ComputeFaceCurvatures(double MinThreshold, double RatioThreshold);
	void ComputeFaceDensity();
	void ComputeVertexDensity();
	void BuildVertexMetricTensors(double normal_coefficient);
	void BuildFaceMetricTensors(double normal_coefficient);
	// curvature

	// IO
	void LoadPrincipalCurvatures(string filename, double MinThreshold, double RatioThreshold, double IsotropicCoefficient);
	// IO

	// curvature estimation [Meyer et al 2002]
	double vertex_voronoi_area_Meyer(Vertex_handle vh);
	double vertex_voronoi_area_new(Vertex_handle vh);
	double vertex_angle(Vertex_handle vh);

	bool vertex_gaussian_curvature_Meyer(Vertex_handle vh);
	void ComputeVertexGaussianCurvature_Meyer();

	bool vertex_gaussian_curvature_new(Vertex_handle vh);
	void ComputeVertexGaussianCurvature_new();
		
	Wm4::Vector3d vertex_cotan(Vertex_handle vh);
	bool vertex_mean_curvature_normal_Meyer(Vertex_handle vh);
	void ComputeVertexMeanCurvatureNormal_Meyer();

	bool vertex_principal_curvature_Meyer(Vertex_handle vh);
	void ComputeVertexPrincipalCurvature_Meyer();

	void EstimateNormalCurvature_Meyer();

	void MeyerEstimationReverseOrientation();

	double m_max_gaussiancurvature_Meyer;
	double m_min_gaussiancurvature_Meyer;

	double m_max_gaussiancurvature_new;
	double m_min_gaussiancurvature_new;

	double m_max_meancurvature_Meyer;
	double m_min_meancurvature_Meyer;

	double m_max_maxcurvature_Meyer;
	double m_min_maxcurvature_Meyer;

	double m_max_mincurvature_Meyer;
	double m_min_mincurvature_Meyer;
	
	// curvature estimation [Meyer et al 2002]


	// discrete convex-concave
	int number_of_bad_vertices;
	int number_of_real_bad_vertices;
	// discrete convex-concave


public:
	
	Mesh_domain * domain;
	Triangulation dt;
	Triangulation dt_dt;



	unsigned NearestVertexId(Vector3d p);
	Vector3d NearestVertex(Vector3d p);
	Vector3d NearestPoint(Vector3d p);

	double maxhausdorff_distance;
	double bb_diagonal_length;
	double mean_distance_to_mesh;
	double mean_distance_to_ma;
	double mean_hausdorff_distance;

	void computesimpledt();
	void computedt();
	void markpoles();

public:
	int LocalFlipCount(Vertex_handle vh);
	bool insideout[50][50][50];
};

struct Facet_normal // (functor)
{
	template <class Facet>
	void operator()(Facet & f)
	{
		typename Facet::Normal sum = CGAL::NULL_VECTOR;
		typename Facet::Halfedge_around_facet_circulator h = f.facet_begin();
		do
		{
			typename Facet::Normal normal = CGAL::cross_product(
				h->next()->vertex()->point() - h->vertex()->point(),
				h->next()->next()->vertex()->point() - h->next()->vertex()->point());
			double sqnorm = normal * normal;
			if(sqnorm != 0)
				normal = normal / (float)std::sqrt(sqnorm);
			sum = sum + normal;
		}while( ++h != f.facet_begin());
		double sqnorm = sum * sum;
		if(sqnorm != 0)
			f.normal = sum / std::sqrt(sqnorm);
		else
			f.normal = CGAL::NULL_VECTOR;
	}
};

struct Vertex_normal // (functor)
{
	template <class Vertex>
	void operator()(Vertex & v)
	{
		typename Vertex::Normal normal = CGAL::NULL_VECTOR;
		typename Vertex::Halfedge_around_vertex_const_circulator pHalfedge = v.vertex_begin();
		typename Vertex::Halfedge_around_vertex_const_circulator begin = pHalfedge;
		CGAL_For_all(pHalfedge, begin)
			if(!pHalfedge->is_border())
				normal = normal + pHalfedge->facet()->normal;
		double sqnorm  = normal * normal;
		if(sqnorm != 0)
			v.normal = normal / (float)std::sqrt(sqnorm);
		else
			v.normal = CGAL::NULL_VECTOR;
	}
};

typedef MPMesh Mesh;

typedef simple_kernel::FT FT;
typedef simple_kernel::RT RT;
typedef simple_kernel::Point_3 Point;
typedef simple_kernel::Vector_3 Vector;
typedef simple_kernel::Direction_3 Direction;

typedef Mesh::Facet_handle Facet_handle;
typedef Mesh::Vertex_handle Vertex_handle;
typedef Mesh::Halfedge_handle Halfedge_handle;

typedef Mesh::Facet_iterator Facet_iterator;
typedef Mesh::Vertex_iterator Vertex_iterator;
typedef Mesh::Halfedge_iterator Halfedge_iterator;

typedef Mesh::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
typedef Mesh::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;

typedef Mesh::Edge_iterator Edge_iterator;






#endif // _MESH_H
