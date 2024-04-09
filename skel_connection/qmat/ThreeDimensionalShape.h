#ifndef _THREE_DIMENSIONAL_SHAPE_H_
#define _THREE_DIMENSIONAL_SHAPE_H_

#include "Mesh.h"
#include "NonManifoldMesh/nonmanifoldmesh.h"
#include "SlabMesh.h"

class ThreeDimensionalShape
{
public:
	ThreeDimensionalShape() : slab_initial(false) {}

	void ComputeInputNMM();
	
	// load the user simplified ma
	void LoadInputNMM(std::string fname);

	long LoadSlabMesh();

	// initial the matrix of each face and vertex for slab mesh
	void InitialSlabMesh();	

	double NearestPoint(Vector3d point, unsigned vid);

	// compute the Hausdorff distance
	void ComputeHausdorffDistance();	

	void PruningSlabMesh();

public:
	Mesh input;		// the mesh of the input shape

	unsigned num_vor_v, num_vor_e, num_vor_f;

	NonManifoldMesh input_nmm;

	SlabMesh slab_mesh;

	bool slab_initial;

};
#endif // _THREE_DIMENSIONAL_SHAPE_H_