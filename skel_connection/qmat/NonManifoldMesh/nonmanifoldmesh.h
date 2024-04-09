#ifndef _NONMANIFOLDMESH_H
#define _NONMANIFOLDMESH_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

#include "LinearAlgebra/Wm4Vector.h"
#include "LinearAlgebra/Wm4Matrix.h"
#include "Mesh.h"

class NonManifoldMesh_Vertex
{
public:
	std::set<unsigned> edges_; // edge list
	std::set<unsigned> faces_; // triangle list
	bool HasEdge(unsigned eid){return (edges_.find(eid) != edges_.end());}
	bool HasFace(unsigned fid){return (faces_.find(fid) != faces_.end());}

public:
	//Wm4::Vector3d pos;
	int tag;
	std::set<unsigned> bplist; // boundary point list
	
	// distance criterion
	//double radius;

	// vetex-related sphere
	Sphere sphere;

	// angle criterion
	double angle;

	// lambda criterion
	double lambda;

	// scale axis 
	double scale_axis_factor;

	// union of balls criterion
	double unionofballs_error;

	// matrix of A
	//Wm4::Matrix4d A;

	//// matrix of b
	//Wm4::Vector4d b;
	
	// the Q matrix for each vertex
	Wm4::Matrix4d Q;

	bool is_pole;
	bool is_non_manifold;
	bool is_disk;
	bool is_boundary;

	std::set<unsigned> pole_bplist;

	int vmanifoldid;

	// for edge contraction error evaluation 
	double v_evaluated_distance_error_envelope;

	// 
	vector<Sphere> mergedspheres;
};

class NonManifoldMesh_Edge
{
public:
	std::pair<unsigned,unsigned> vertices_; // vertex list
	std::set<unsigned> faces_; // triangle list
	bool HasVertex(unsigned vid){return ( (vertices_.first == vid) || (vertices_.second == vid));}
	bool HasFace(unsigned fid){return (faces_.find(fid) != faces_.end());}

public:
	NonManifoldMesh_Edge(){validenvelope = true;}
public:
	
	int tag;
	double length;

	bool validenvelope;

	// merging
	//Wm4::Vector2d contract_pos;
	//double contract_rad;

	//double contraction_error;
	//unsigned contraction_target_id;

	double contraction_error_f_s; // first -> second;
	double contraction_error_s_f; // second -> first;

	double de_contraction_error_f_s; // first -> second;
	double de_contraction_error_s_f; // second -> first;


	Cone cone;
	bool valid_cone;

	//int num_adjacent_triangles;
	bool is_dangling;
	bool is_boundary;
	bool is_disk;
	bool is_non_manifold;

	int emanifoldid;

	bool has_footpoints;

	double collapse_cost;
	Sphere sphere;
	Wm4::Matrix4d Q;
	//double de_error;
	//unsigned de_target;
};

class NonManifoldMesh_Face
{
public:
	bool validenvelope;
public:
	std::set<unsigned> vertices_; // vertex list
	std::set<unsigned> edges_; // edge list
	bool HasVertex(unsigned vid){return (vertices_.find(vid) != vertices_.end());}
	bool HasEdge(unsigned eid){return (edges_.find(eid) != edges_.end());}

public:
	Wm4::Vector3d centroid;
	Wm4::Vector3d normal;
	int tag;

	int fmanifoldid;

	bool is_boundary;

	// Mp Matrix for each plane
	Wm4::Matrix4d Mp;

public:
	SimpleTriangle st[2];
	bool valid_st;

	bool has_footpoints;
};

typedef std::pair<bool, NonManifoldMesh_Vertex*> Bool_VertexPointer;
typedef std::pair<bool, NonManifoldMesh_Edge*> Bool_EdgePointer;
typedef std::pair<bool, NonManifoldMesh_Face*> Bool_FacePointer;

class cmp_qem{
public:
	bool operator() (NonManifoldMesh_Edge &a, NonManifoldMesh_Edge &b){
		return a.collapse_cost > b.collapse_cost;
	}
};

class NonManifoldMesh
{
public:
	Mesh * pmesh;
	Mesh_domain * domain;
public:
	double diameter;
public:
	double seconds;
public:
	std::string meshname;
public:
	NonManifoldMesh(){error_threshold = 1e-4; seconds = 0.;numVertices = numEdges = numFaces = 0;}
public:
	std::vector<Bool_VertexPointer> vertices;
	std::vector<Bool_EdgePointer> edges;
	std::vector<Bool_FacePointer> faces;


public:
	std::vector< Sphere > SampleSpheres(unsigned fid, unsigned num)
	{
		Sphere s[3];
		unsigned count = 0;
		for(set<unsigned>::iterator si = faces[fid].second->vertices_.begin(); si != faces[fid].second->vertices_.end(); si ++)
			s[count++] = vertices[*si].second->sphere;
		
		std::vector<Sphere> vs;
		for(unsigned i = 0; i <= num; i ++)
			for(unsigned j = 0; j <= num-i; j ++)
			{
				double u = (double)i/num;
				double v = (double)j/num;
				double w = 1. - u - v;
				Sphere samp_s;
				samp_s = s[0] * u + s[1] * v + s[2] * w;
				vs.push_back(samp_s);
			}
		return vs;
	}

public:
	std::vector<SamplePoint> BoundaryPoints;
	unsigned num_f_manifolds;
	unsigned num_e_manifolds;
	unsigned num_v_manifolds;

	unsigned numVertices;
	unsigned numEdges;
	unsigned numFaces;

	double min_edge_contraction_error;
	double max_edge_contraction_error;

	double max_boundarypoint_footpoint_distance;
	double max_boundarypoint_footpoint_idx;

	double mean_boundarypoint_footpoint_distance;

public:
	void AdjustStorage();

public:
	bool ValidVertex(unsigned vid);
	bool Edge(unsigned vid0, unsigned vid1, unsigned & eid);
	bool Face(const std::set<unsigned> & vset, unsigned & fid);
	void UpdateCentroid(unsigned fid);
	void ComputeFacesCentroid();
	void UpdateNormal(unsigned fid);
	void ComputeFacesNormal();
	void GetNeighborVertices(unsigned vid, std::set<unsigned> & neighborvertices);
	void GetLinkedEdges(unsigned eid, std::set<unsigned> & neighboredges);
	void GetAdjacentFaces(unsigned fid, std::set<unsigned> & neighborfaces);
	bool Contractible(unsigned vid_src, unsigned vid_tgt);
	bool MergeVertices(unsigned vid_src1, unsigned vid_src2, unsigned &vid_tgt);

	bool non_pole_only;

	double error_threshold;

public:
	void Export(std::string fname);

	void DeleteFace(unsigned fid);
	void DeleteEdge(unsigned eid);
	void DeleteVertex(unsigned vid);

	void InsertVertex(NonManifoldMesh_Vertex *vertex, unsigned &vid);
	void InsertEdge(unsigned vid0, unsigned vid1, unsigned & eid);
	void InsertFace(std::set<unsigned> vset);

	unsigned VertexIncidentEdgeCount(unsigned vid);
	unsigned VertexIncidentFaceCount(unsigned vid);
	unsigned EdgeIncidentFaceCount(unsigned eid);

public:
	void ComputeEdgeCone(unsigned eid);
	void ComputeEdgesCone();

	void ComputeFaceSimpleTriangles(unsigned fid);
	void ComputeFacesSimpleTriangles();

public:
	std::priority_queue<NonManifoldMesh_Edge, vector<NonManifoldMesh_Edge>, cmp_qem> edge_collapses_queue;
	void initCollapseQueue();
	void EvaluateEdgeCollapseCost(unsigned eid);
	void MinCostEdgeCollapse(unsigned & eid);

	void Simplify(int threshold);
};


#endif // _NONMANIFOLDMESH_H
