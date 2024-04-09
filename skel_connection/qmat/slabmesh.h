#ifndef _SLABMESH_H
#define _SLABMESH_H

#include "PrimMesh.h"


class SlabPrim
{
public:
	Wm4::Matrix4d slab_A;
	Wm4::Vector4d slab_b;
	double slab_c;

	Wm4::Matrix4d add_A;
	Wm4::Vector4d add_b;
	double add_c;

	// hyperbolic weight
	double hyperbolic_weight;

	SlabPrim() : slab_c(0.0), add_c(0.0), hyperbolic_weight(0.0){}
};

class SlabVertex : public PrimVertex, public SlabPrim
{
public:
	bool is_pole;
	bool is_non_manifold;
	bool is_disk;
	bool is_boundary;
	bool is_selected;
};

class SlabEdge : public PrimEdge, public SlabPrim
{
};

class SlabFace : public PrimFace, public SlabPrim
{
};

typedef std::pair<bool, SlabVertex*> Bool_SlabVertexPointer;
typedef std::pair<bool, SlabEdge*> Bool_SlabEdgePointer;
typedef std::pair<bool, SlabFace*> Bool_SlabFacePointer;

class SlabMesh : public PrimMesh
{
public:
	std::vector<Bool_SlabVertexPointer> vertices;
	std::vector<Bool_SlabEdgePointer> edges;
	std::vector<Bool_SlabFacePointer> faces;

public:
	// 1. preserve method one
	// 2. preserve method two
	// 3. preserve method three
	int preserve_boundary_method;

	// 0. no hyperbolic weight
	// 1. add hyperbolic distance to the related edges
	// 2. add hyperbolic area to the related face
	// 3. add ratio of hyperbolic and Euclid to the related edges
	int hyperbolic_weight_type;

	// 1.compute the boundary vertices only
	// 2.compute the boundary vertices and its related vertices
	int boundary_compute_scale;

	// the influence factor to control the ratio between hyperbolic and Euclid distance
	double k;

	// whether clear the error before initialization
	bool clear_error;

	bool preserve_saved_vertices;

	bool compute_hausdorff;

	bool prevent_inversion;

	bool initial_boundary_preserve;

	double m_min[3];
	double m_max[3];

	unsigned simplified_inside_edges;
	unsigned simplified_boundary_edges;

	double bound_weight;

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
	void UpdateVertexNormal(unsigned vid);
	void ComputeVerticesNormal();
	void GetNeighborVertices(unsigned vid, std::set<unsigned> & neighborvertices);
	void GetLinkedEdges(unsigned eid, std::set<unsigned> & neighboredges);
	void GetAdjacentFaces(unsigned fid, std::set<unsigned> & neighborfaces);
	bool Contractible(unsigned vid_src, unsigned vid_tgt);
	bool Contractible(unsigned vid_src1, unsigned vid_src2, Vector3d &v_tgt);
	bool MergeVertices(unsigned vid_src1, unsigned vid_src2, unsigned &vid_tgt);
	vector<string> my_split(const string &str, const string &pattern); // Zhiyang added.
	double p_distance(double x1, double y1, double z1, double x2, double y2, double z2); // Zhiyang added.

	unsigned VertexIncidentEdgeCount(unsigned vid);
	unsigned VertexIncidentFaceCount(unsigned vid);
	unsigned EdgeIncidentFaceCount(unsigned eid);


public:
	void DeleteFace(unsigned fid);
	void DeleteEdge(unsigned eid);
	void DeleteVertex(unsigned vid);

	void InsertVertex(SlabVertex *vertex, unsigned &vid);
	void InsertEdge(unsigned vid0, unsigned vid1, unsigned & eid);
	void InsertFace(std::set<unsigned> vset);

	void ComputeEdgeCone(unsigned eid);
	void ComputeEdgesCone();

	void ComputeVertexProperty(unsigned vid);
	void ComputeVerticesProperty();
	void ComputeFaceSimpleTriangles(unsigned fid);
	void ComputeFacesSimpleTriangles();

public:
	void initBoundaryCollapseQueue();
	void initCollapseQueue();
	void Simplify_with_Selected_Pole(int threshold, vector<vector<double>> &selected_pole);
	void Simplify(int threshold);
	void SimplifyBoudary(int threshold);
	bool MinCostBoundaryEdgeCollapse(unsigned & eid);
	bool MinCostEdgeCollapse(unsigned & eid);
	void EvaluateEdgeCollapseCost(unsigned eid);
	void EvaluateEdgeHausdorffCost(unsigned eid);
	void ReEvaluateEdgeHausdorffCost(unsigned eid);
	void readMA_ball_diff_radius(string objpath, vector<vector<double> >& vset, vector<vector<int> >& eset, vector<vector<int> >& fset);
public: 
	void DistinguishVertexType();
	unsigned GetSavedPointNumber();
	unsigned GetConnectPointNumber();
	void InsertSavedPoint(unsigned vid);
	double NearestPoint(Vector3d point, unsigned vid);

public:
	void PreservBoundaryMethodOne();
	void PreservBoundaryMethodTwo();
	void PreservBoundaryMethodThree();
	void PreservBoundaryMethodFour();
	void PreservBoundaryMethodFive();

	void GetEnvelopeSet(const Vector4d & lamder, const set<unsigned> & neighbor_v, const set< std::set<unsigned> > & adj_faces, vector<Sphere> & sph_vec, vector<Cone> & con_vec, vector<SimpleTriangle> & st_vec);
	double EvaluateVertexDistanceErrorEnvelope(Vector4d & lamder, set<unsigned> & neighbor_vertices, set< set<unsigned> > & neighbor_faces, set<unsigned> & bplist);

	double GetHyperbolicLength(unsigned eid);
	double GetRatioHyperbolicEuclid(unsigned eid);

	void ExportSimplifyResult();
	void Export(string fname);
	void Export_OBJ(string fname);


public:
	void clear();
	void RecomputerVertexType();
	void computebb();

	void CleanIsolatedVertices();
	void InitialTopologyProperty(unsigned vid);
	void InitialTopologyProperty();
};

#endif