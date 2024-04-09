#ifndef _PRIMMESH_H
#define _PRIMMESH_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <queue>
#include <string>
#include "LinearAlgebra/Wm4Vector.h"
#include "Mesh.h"
using namespace std;
class PrimVertex{
public:
	std::set<unsigned> edges_;	//edge list
	std::set<unsigned> faces_; //triangle list
	bool HasEdge(unsigned eid){return (edges_.find(eid) != edges_.end());}
	bool HasFace(unsigned fid){return (faces_.find(fid) != faces_.end());}
	PrimVertex() : fake_boundary_vertex(false), boundary_vertex(false), saved_vertex(false), 
		non_manifold_vertex(false), collaspe_weight(0.0), boundVec(Vector3d(0, 0, 0)), mean_square_error(0.0), related_face(0){};
	virtual ~PrimVertex(){};

public:
	// record the information of related square error 
	double mean_square_error;
	unsigned related_face;

public:
	// vertex property
	bool fake_boundary_vertex;
	bool boundary_vertex;
	bool saved_vertex;
	bool non_manifold_vertex;
	set<unsigned> boundary_edge_vec;
	Wm4::Vector3d boundVec;

public:
	Sphere sphere;
	double collaspe_weight;
	unsigned index;
	std::set<unsigned> bplist; // boundary point list of the original mesh

	int tag;
	Wm4::Vector3d normal;
};


class PrimEdge{
public:
	std::pair<unsigned,unsigned> vertices_; // vertex list
	std::set<unsigned> faces_; // triangle list
	bool HasVertex(unsigned vid){return ( (vertices_.first == vid) || (vertices_.second == vid));}
	bool HasFace(unsigned fid){return (faces_.find(fid) != faces_.end());}
	PrimEdge(): fake_boundary_edge(false), boundary_edge(false), non_manifold_edge(false), topo_contractable(true){};
	virtual ~PrimEdge(){};

public:
	Sphere sphere;
	Cone cone;
	bool valid_cone;

public:
	double collapse_cost;
	double qem_error;
	unsigned index;

public:
	bool fake_boundary_edge;
	bool boundary_edge;
	bool non_manifold_edge;
	bool topo_contractable;
};


class PrimFace{
public:
	std::set<unsigned> vertices_; // vertex list
	std::set<unsigned> edges_; // edge list
	bool HasVertex(unsigned vid){return (vertices_.find(vid) != vertices_.end());}
	bool HasEdge(unsigned eid){return (edges_.find(eid) != edges_.end());}
	virtual ~PrimFace(){};

public:
	unsigned index;
	Sphere centroid;
	Wm4::Vector3d normal;

public:
	SimpleTriangle st[2];
	bool valid_st;
};


class EdgeInfo
{
	friend bool operator < (EdgeInfo e1, EdgeInfo e2)
	{
		return e1.collapse_cost > e2.collapse_cost;
	}
public:
	unsigned edge_num;
	double collapse_cost;

	EdgeInfo(){};
	EdgeInfo(unsigned num, double cost){edge_num = num; collapse_cost = cost;};
};

template<class Real>
class cmpByValue{
public:
	bool operator() (const pair<unsigned, Real> lhs, const pair<unsigned, Real> rhs){
		return lhs.second < rhs.second;
	}
};

template<class Real>
class cmpByBiggerValue{
public:
	bool operator() (const pair<unsigned, Real> lhs, const pair<unsigned, Real> rhs){
		return lhs.second > rhs.second;
	}
};

//typedef std::pair<bool, PrimVertex*> Bool_PrimVertexPointer;
//typedef std::pair<bool, PrimEdge*> Bool_PrimEdgePointer;
//typedef std::pair<bool, PrimFace*> Bool_PrimFacePointer;

class PrimMesh{
public:
	Mesh * pmesh;
	Mesh_domain * domain;
	std::string meshname;
	//public:
	//	std::vector<Bool_PrimVertexPointer> vertices;
	//	std::vector<Bool_PrimEdgePointer> edges;
	//	std::vector<Bool_PrimFacePointer> faces;

public:
	unsigned num_f_manifolds;
	unsigned num_e_manifolds;
	unsigned num_v_manifolds;

	unsigned iniNumVertices;
	unsigned iniNumEdges;
	unsigned iniNumFaces;

	unsigned numVertices;
	unsigned numEdges;
	unsigned numFaces;

	unsigned simp_record;
	unsigned temp_value;

	std::vector<pair<unsigned, unsigned>> bound_vector;
	std::vector<pair<unsigned, double>> boundVec_vector;
	std::set<unsigned> boundary_vertexes;

	double meanhausdorff_distance;
	double maxhausdorff_distance;
	double initialhausdorff_distance;
	Triangulation dt;

	// 0 for medial mesh
	// 1 for slab mesh
	unsigned type;

	double start_multi;
	double end_multi;

	double max_mean_squre_error;
	double min_mean_squre_error;

	long merge_time;
	long cal_time;

public:
	PrimMesh() : iniNumVertices(0), iniNumEdges(0), iniNumFaces(0), merge_time(0.0), cal_time(0.0), start_multi(5.0), end_multi(10.0), max_mean_squre_error(0.0), min_mean_squre_error(0.0){};
	virtual ~PrimMesh(){};

public:
	std::priority_queue<EdgeInfo> edge_collapses_queue;
	std::priority_queue<EdgeInfo> boundary_edge_collapses_queue;

};

#endif