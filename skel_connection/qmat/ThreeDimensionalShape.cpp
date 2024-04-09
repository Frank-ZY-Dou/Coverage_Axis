#include "ThreeDimensionalShape.h"

#include <QString>

void ThreeDimensionalShape::ComputeInputNMM()
{

	cout << "here!" << endl;
	input_nmm.numVertices = 0;
	input_nmm.numEdges = 0;
	input_nmm.numFaces = 0;

	//
	input_nmm.vertices.clear();
	input_nmm.edges.clear();
	input_nmm.faces.clear();

	Triangulation * pt = &(input.dt);

	//loading
	//slab_mesh.numVertices = 0;
	//slab_mesh.numEdges = 0;
	//slab_mesh.numFaces = 0;
	
	num_vor_v = 0;
	num_vor_e = 0;
	num_vor_f = 0;

	double len[4];
	len[0] = input.m_max[0] - input.m_min[0];
	len[1] = input.m_max[1] - input.m_min[1];
	len[2] = input.m_max[2] - input.m_min[2];
	len[3] = sqrt(len[0]*len[0]+len[1]*len[1]+len[2]*len[2]);
	input_nmm.diameter = len[3];
	int mysample_counter = 0;
	for (Finite_vertices_iterator_t fvi = pt->finite_vertices_begin(); fvi != pt->finite_vertices_end(); fvi++) {
		input_nmm.BoundaryPoints.push_back(SamplePoint(fvi->point()[0], fvi->point()[1], fvi->point()[2]));
	}
	cout << "boundary size: "<< input_nmm.BoundaryPoints.size() << endl;
	
	// Here is the key...


	int mas_vertex_count(0);
	//
	//slab_mesh.maxhausdorff_distance = 0;
	for(Finite_cells_iterator_t fci = pt->finite_cells_begin(); fci != pt->finite_cells_end(); fci ++)
	{
		if(fci->info().inside == false)
		{
			fci->info().tag = -1;
			continue;
		}
		fci->info().tag = mas_vertex_count ++;
		

		Bool_VertexPointer bvp;
		bvp.first = true;
		bvp.second = new NonManifoldMesh_Vertex;
		(*bvp.second).sphere.center = to_wm4(CGAL::circumcenter(pt->tetrahedron(fci)));
		(*bvp.second).is_pole = fci->info().is_pole;
		for(unsigned k = 0; k < 4; k ++)
			(*bvp.second).bplist.insert(fci->vertex(k)->info().id);
		(*bvp.second).sphere.radius = pt->TetCircumRadius(pt->tetrahedron(fci));
		//(*bvp.second).sphere.radius = fci->info().dist_center_to_boundary; // make sure that all the spheres are inside the domain
		//bvp.second->sphere.center = bvp.second->sphere.center;
		//bvp.second->sphere.radius = bvp.second->sphere.radius;
		input_nmm.vertices.push_back(bvp);
		input_nmm.numVertices ++;
		num_vor_v ++;

		// 计算slab vertices 
		//Bool_SlabVertexPointer bsvp2;
		//bsvp2.first = true;
		//bsvp2.second = new SlabVertex;
		//(*bsvp2.second).sphere.center[0] = (*bvp.second).sphere.center.X();
		//(*bsvp2.second).sphere.center[1] = (*bvp.second).sphere.center.Y();
		//(*bsvp2.second).sphere.center[2] = (*bvp.second).sphere.center.Z();
		//(*bsvp2.second).sphere.radius = (*bvp.second).sphere.radius;
		////(*bsvp2.second).sphere.radius = r * 1.2;
		//(*bsvp2.second).index = slab_mesh.vertices.size();
		//(*bsvp2.second).bplist = (*bvp.second).bplist;
		//slab_mesh.vertices.push_back(bsvp2);
		//slab_mesh.numVertices ++;

		//double min_dis = DBL_MAX;
		//for (set<unsigned>::iterator si = (*bvp.second).bplist.begin(); si != (*bvp.second).bplist.end(); si++)
		//{
		//	Vector3d bou_ver(input.pVertexList[*si]->point()[0], input.pVertexList[*si]->point()[1], input.pVertexList[*si]->point()[2]);
		//	Sphere ma_ver = bsvp2.second->sphere;
		//	double temp_length = abs((bou_ver - ma_ver.center).Length() - ma_ver.radius);
		//	min_dis = min(min_dis, temp_length);
		//}
		//slab_mesh.maxhausdorff_distance = max(slab_mesh.maxhausdorff_distance, min_dis);
	}
	//slab_mesh.initialhausdorff_distance = slab_mesh.maxhausdorff_distance;
	
	for(Finite_facets_iterator_t ffi = pt->finite_facets_begin(); ffi != pt->finite_facets_end(); ffi ++)
	{
		Triangulation::Object o = pt->dual(*ffi);
		if(const Triangulation::Segment *s = CGAL::object_cast<Triangulation::Segment>(&o))
		{
			if( (ffi->first->info().inside == false) || (pt->mirror_facet(*ffi).first->info().inside == false) )
				continue;
			Bool_EdgePointer bep;
			bep.first = true;
			bep.second = new NonManifoldMesh_Edge;
			(*bep.second).vertices_.first = ffi->first->info().tag;
			(*bep.second).vertices_.second = pt->mirror_facet(*ffi).first->info().tag;
			(*input_nmm.vertices[ffi->first->info().tag].second).edges_.insert(input_nmm.edges.size());
			(*input_nmm.vertices[pt->mirror_facet(*ffi).first->info().tag].second).edges_.insert(input_nmm.edges.size());
			input_nmm.edges.push_back(bep);
			input_nmm.numEdges ++;
			num_vor_e ++;

			// 计算slab edges
			//Bool_SlabEdgePointer bsep2;
			//bsep2.first = true;
			//bsep2.second = new SlabEdge;
			//(*bsep2.second).vertices_.first = (*bep.second).vertices_.first;
			//(*bsep2.second).vertices_.second = (*bep.second).vertices_.second;
			//(*slab_mesh.vertices[(*bsep2.second).vertices_.first].second).edges_.insert(slab_mesh.edges.size());
			//(*slab_mesh.vertices[(*bsep2.second).vertices_.second].second).edges_.insert(slab_mesh.edges.size());
			//(*bsep2.second).index = slab_mesh.edges.size();
			//slab_mesh.edges.push_back(bsep2);
			//slab_mesh.numEdges ++;
		}
	}
	
	for(Finite_edges_iterator_t fei = pt->finite_edges_begin(); fei != pt->finite_edges_end(); fei ++)
	{
		bool all_finite_inside = true;
		std::vector<Cell_handle_t> vec_ch;
		Cell_circulator_t cc = pt->incident_cells(*fei);
		do
		{
			if(pt->is_infinite(cc))
				all_finite_inside = false;
			else if(cc->info().inside == false)
				all_finite_inside = false;
			vec_ch.push_back(cc++);
		}while(cc != pt->incident_cells(*fei));
		if(!all_finite_inside)
			continue;

		for(unsigned k = 2; k < vec_ch.size() - 1; k ++)
		{
			Bool_EdgePointer bep;
			bep.first = true;
			bep.second = new NonManifoldMesh_Edge;
			(*bep.second).vertices_.first = vec_ch[0]->info().tag;
			(*bep.second).vertices_.second = vec_ch[k]->info().tag;
			(*input_nmm.vertices[vec_ch[0]->info().tag].second).edges_.insert(input_nmm.edges.size());
			(*input_nmm.vertices[vec_ch[k]->info().tag].second).edges_.insert(input_nmm.edges.size());
			input_nmm.edges.push_back(bep);
			input_nmm.numEdges ++;

			// 计算slab edges
			//Bool_SlabEdgePointer bsep2;
			//bsep2.first = true;
			//bsep2.second = new SlabEdge;
			//(*bsep2.second).vertices_.first = (*bep.second).vertices_.first;
			//(*bsep2.second).vertices_.second = (*bep.second).vertices_.second;
			//(*slab_mesh.vertices[(*bsep2.second).vertices_.first].second).edges_.insert(slab_mesh.edges.size());
			//(*slab_mesh.vertices[(*bsep2.second).vertices_.second].second).edges_.insert(slab_mesh.edges.size());
			//(*bsep2.second).index = slab_mesh.edges.size();
			//slab_mesh.edges.push_back(bsep2);
			//slab_mesh.numEdges ++;
		}

		for(unsigned k = 1; k < vec_ch.size() - 1; k ++)
		{
			Bool_FacePointer bfp;
			bfp.first = true;
			bfp.second = new NonManifoldMesh_Face;
			unsigned vid[3];
			vid[0] = vec_ch[0]->info().tag;
			vid[1] = vec_ch[k]->info().tag;
			vid[2] = vec_ch[k+1]->info().tag;
			(*bfp.second).vertices_.insert(vec_ch[0]->info().tag);
			(*bfp.second).vertices_.insert(vec_ch[k]->info().tag);
			(*bfp.second).vertices_.insert(vec_ch[k+1]->info().tag);
			unsigned eid[3];
			if(input_nmm.Edge(vid[0],vid[1],eid[0]))
				(*bfp.second).edges_.insert(eid[0]);
			if(input_nmm.Edge(vid[0],vid[2],eid[1]))
				(*bfp.second).edges_.insert(eid[1]);
			if(input_nmm.Edge(vid[1],vid[2],eid[2]))
				(*bfp.second).edges_.insert(eid[2]);
			input_nmm.vertices[vid[0]].second->faces_.insert(input_nmm.faces.size());
			input_nmm.vertices[vid[1]].second->faces_.insert(input_nmm.faces.size());
			input_nmm.vertices[vid[2]].second->faces_.insert(input_nmm.faces.size());
			input_nmm.edges[eid[0]].second->faces_.insert(input_nmm.faces.size());
			input_nmm.edges[eid[1]].second->faces_.insert(input_nmm.faces.size());
			input_nmm.edges[eid[2]].second->faces_.insert(input_nmm.faces.size());
			input_nmm.faces.push_back(bfp);
			input_nmm.numFaces ++;
			num_vor_f ++;

			// 计算slab face
			//Bool_SlabFacePointer bsfp2;
			//bsfp2.first = true;
			//bsfp2.second = new SlabFace;
			//(*bsfp2.second).vertices_.insert(vid[0]);
			//(*bsfp2.second).vertices_.insert(vid[1]);
			//(*bsfp2.second).vertices_.insert(vid[2]);
			//if(slab_mesh.Edge(vid[0],vid[1],eid[0]))
			//	(*bsfp2.second).edges_.insert(eid[0]);
			//if(slab_mesh.Edge(vid[0],vid[2],eid[1]))
			//	(*bsfp2.second).edges_.insert(eid[1]);
			//if(slab_mesh.Edge(vid[1],vid[2],eid[2]))
			//	(*bsfp2.second).edges_.insert(eid[2]);
			//(*bsfp2.second).index = slab_mesh.faces.size();
			//slab_mesh.vertices[vid[0]].second->faces_.insert(slab_mesh.faces.size());
			//slab_mesh.vertices[vid[1]].second->faces_.insert(slab_mesh.faces.size());
			//slab_mesh.vertices[vid[2]].second->faces_.insert(slab_mesh.faces.size());
			//slab_mesh.edges[eid[0]].second->faces_.insert(slab_mesh.faces.size());
			//slab_mesh.edges[eid[1]].second->faces_.insert(slab_mesh.faces.size());
			//slab_mesh.edges[eid[2]].second->faces_.insert(slab_mesh.faces.size());
			//slab_mesh.faces.push_back(bsfp2);
			//slab_mesh.numFaces++;
		}
	}
	input_nmm.Export(input_nmm.meshname);
	
	input_nmm.numVertices = 0;
	input_nmm.numEdges = 0;
	input_nmm.numFaces = 0;
	//input_nmm.ComputeFacesNormal();
	//input_nmm.ComputeFacesCentroid();
	//input_nmm.ComputeFacesSimpleTriangles();
	//input_nmm.ComputeEdgesCone();

	//slab_mesh.ComputeFacesCentroid();
	//slab_mesh.ComputeFacesNormal();
	//slab_mesh.ComputeVerticesNormal();

}

void ThreeDimensionalShape::LoadInputNMM(std::string fname){
	std::ifstream mastream(fname.c_str());
	NonManifoldMesh newinputnmm;
	newinputnmm.numVertices = 0;
	newinputnmm.numEdges = 0;
	newinputnmm.numFaces = 0;
	int nv, ne, nf;
	mastream >> nv >> ne >> nf;

	// slab mesh
	slab_mesh.numVertices = 0;
	slab_mesh.numEdges = 0;
	slab_mesh.numFaces = 0;

	double len[4];
	len[0] = input.m_max[0] - input.m_min[0];
	len[1] = input.m_max[1] - input.m_min[1];
	len[2] = input.m_max[2] - input.m_min[2];
	len[3] = sqrt(len[0]*len[0]+len[1]*len[1]+len[2]*len[2]);
	newinputnmm.diameter = len[3];
	slab_mesh.bound_weight = 0.1; 

	for(unsigned i = 0; i < input.pVertexList.size(); i ++)
		newinputnmm.BoundaryPoints.push_back(SamplePoint(
		input.pVertexList[i]->point()[0],
		input.pVertexList[i]->point()[1],
		input.pVertexList[i]->point()[2]
	));

	for(unsigned i = 0; i < nv; i ++)
	{
		char ch;
		double x,y,z,r;
		mastream >> ch >> x >> y >> z >> r;

		//Bool_VertexPointer bvp;
		//bvp.first = true;
		//bvp.second = new NonManifoldMesh_Vertex;
		//(*bvp.second).sphere.center[0] = x / input.bb_diagonal_length;
		//(*bvp.second).sphere.center[1] = y / input.bb_diagonal_length;
		//(*bvp.second).sphere.center[2] = z / input.bb_diagonal_length;
		//(*bvp.second).sphere.radius = r / input.bb_diagonal_length;
		//newinputnmm.vertices.push_back(bvp);
		//newinputnmm.numVertices ++;

		// handle the slab mesh
		Bool_SlabVertexPointer bsvp2;
		bsvp2.first = true;
		bsvp2.second = new SlabVertex;
		(*bsvp2.second).sphere.center[0] = x / input.bb_diagonal_length;
		(*bsvp2.second).sphere.center[1] = y / input.bb_diagonal_length;
		(*bsvp2.second).sphere.center[2] = z / input.bb_diagonal_length;
		(*bsvp2.second).sphere.radius = r / input.bb_diagonal_length;
		(*bsvp2.second).index = slab_mesh.vertices.size();
		slab_mesh.vertices.push_back(bsvp2);
		slab_mesh.numVertices ++;
	}

	for(unsigned i = 0; i < ne; i ++)
	{
		char ch;
		unsigned ver[2];
		mastream >> ch;
		mastream >> ver[0];
		mastream >> ver[1];

		//Bool_EdgePointer bep;
		//bep.first = true;
		//bep.second = new NonManifoldMesh_Edge;
		//(*bep.second).vertices_.first = ver[0];
		//(*bep.second).vertices_.second = ver[1];
		//(*newinputnmm.vertices[(*bep.second).vertices_.first].second).edges_.insert(newinputnmm.edges.size());
		//(*newinputnmm.vertices[(*bep.second).vertices_.second].second).edges_.insert(newinputnmm.edges.size());
		//newinputnmm.edges.push_back(bep);
		//newinputnmm.numEdges ++;

		// handle the slab mesh
		Bool_SlabEdgePointer bsep2;
		bsep2.first = true;
		bsep2.second = new SlabEdge;
		(*bsep2.second).vertices_.first = ver[0];
		(*bsep2.second).vertices_.second = ver[1];
		(*slab_mesh.vertices[(*bsep2.second).vertices_.first].second).edges_.insert(slab_mesh.edges.size());
		(*slab_mesh.vertices[(*bsep2.second).vertices_.second].second).edges_.insert(slab_mesh.edges.size());
		(*bsep2.second).index = slab_mesh.edges.size();
		slab_mesh.edges.push_back(bsep2);
		slab_mesh.numEdges ++;
	}

	for(unsigned i = 0; i < nf; i ++)
	{
		char ch;
		unsigned vid[3];
		unsigned eid[3];
		mastream >> ch >> vid[0] >> vid[1] >> vid[2];

		//Bool_FacePointer bfp;
		//bfp.first = true;
		//bfp.second = new NonManifoldMesh_Face;
		//(*bfp.second).vertices_.insert(vid[0]);
		//(*bfp.second).vertices_.insert(vid[1]);
		//(*bfp.second).vertices_.insert(vid[2]);
		//if(newinputnmm.Edge(vid[0],vid[1],eid[0]))
		//	(*bfp.second).edges_.insert(eid[0]);
		//if(newinputnmm.Edge(vid[0],vid[2],eid[1]))
		//	(*bfp.second).edges_.insert(eid[1]);
		//if(newinputnmm.Edge(vid[1],vid[2],eid[2]))
		//	(*bfp.second).edges_.insert(eid[2]);
		//newinputnmm.vertices[vid[0]].second->faces_.insert(newinputnmm.faces.size());
		//newinputnmm.vertices[vid[1]].second->faces_.insert(newinputnmm.faces.size());
		//newinputnmm.vertices[vid[2]].second->faces_.insert(newinputnmm.faces.size());
		//newinputnmm.edges[eid[0]].second->faces_.insert(newinputnmm.faces.size());
		//newinputnmm.edges[eid[1]].second->faces_.insert(newinputnmm.faces.size());
		//newinputnmm.edges[eid[2]].second->faces_.insert(newinputnmm.faces.size());
		//newinputnmm.faces.push_back(bfp);
		//newinputnmm.numFaces ++;

		// handle the slab mesh	
		Bool_SlabFacePointer bsfp2;
		bsfp2.first = true;
		bsfp2.second = new SlabFace;
		(*bsfp2.second).vertices_.insert(vid[0]);
		(*bsfp2.second).vertices_.insert(vid[1]);
		(*bsfp2.second).vertices_.insert(vid[2]);
		if(slab_mesh.Edge(vid[0],vid[1],eid[0]))
			(*bsfp2.second).edges_.insert(eid[0]);
		if(slab_mesh.Edge(vid[0],vid[2],eid[1]))
			(*bsfp2.second).edges_.insert(eid[1]);
		if(slab_mesh.Edge(vid[1],vid[2],eid[2]))
			(*bsfp2.second).edges_.insert(eid[2]);
		(*bsfp2.second).index = slab_mesh.faces.size();
		slab_mesh.vertices[vid[0]].second->faces_.insert(slab_mesh.faces.size());
		//slab_mesh.vertices[vid[0]].second->related_face += 2;
		slab_mesh.vertices[vid[1]].second->faces_.insert(slab_mesh.faces.size());
		//slab_mesh.vertices[vid[1]].second->related_face += 2;
		slab_mesh.vertices[vid[2]].second->faces_.insert(slab_mesh.faces.size());
		//slab_mesh.vertices[vid[2]].second->related_face += 2;
		slab_mesh.edges[eid[0]].second->faces_.insert(slab_mesh.faces.size());
		slab_mesh.edges[eid[1]].second->faces_.insert(slab_mesh.faces.size());
		slab_mesh.edges[eid[2]].second->faces_.insert(slab_mesh.faces.size());
		slab_mesh.faces.push_back(bsfp2);
		slab_mesh.numFaces++;
	}

	//newinputnmm.ComputeFacesNormal();
	//newinputnmm.ComputeFacesCentroid();
	//newinputnmm.ComputeFacesSimpleTriangles();
	//newinputnmm.ComputeEdgesCone();
	//input_nmm = newinputnmm;

	slab_mesh.iniNumVertices = slab_mesh.numVertices;
	slab_mesh.iniNumEdges = slab_mesh.numEdges;
	slab_mesh.iniNumFaces = slab_mesh.numFaces;

	slab_mesh.CleanIsolatedVertices();
	slab_mesh.computebb();
	slab_mesh.ComputeFacesCentroid();
	slab_mesh.ComputeFacesNormal();
	slab_mesh.ComputeVerticesNormal();
	slab_mesh.ComputeEdgesCone();
	slab_mesh.ComputeFacesSimpleTriangles();
	slab_mesh.DistinguishVertexType();
}

long ThreeDimensionalShape::LoadSlabMesh()
{
	slab_mesh.clear();
	long startt = clock();
	InitialSlabMesh();
	slab_mesh.initCollapseQueue();
	long endt = clock();
	return endt - startt;

	//if (slab_mesh.clear_error)
	//{
	//	slab_mesh.clear();
	//	long startt = clock();
	//	while(!slab_mesh.boundary_edge_collapses_queue.empty())
	//		slab_mesh.boundary_edge_collapses_queue.pop();

	//	InitialSlabMesh();
	//	slab_mesh.initCollapseQueue();
	//	long endt = clock();
	//	return endt - startt;
	//}
	//else
	//{

	//	slab_mesh.RecomputerVertexType();	
	//	long startt = clock();
	//	while(!slab_mesh.boundary_edge_collapses_queue.empty())
	//		slab_mesh.boundary_edge_collapses_queue.pop();
	//	if (slab_initial == false)
	//	{
	//		switch(slab_mesh.hyperbolic_weight_type)
	//		{
	//		case 1:
	//			InitialSlabMesh();
	//			break;
	//		case 2:
	//			InitialWeightedSlabMesh();
	//			break;
	//		case 3:
	//			InitialSlabMesh();
	//			break;
	//		default:
	//			InitialSlabMesh();
	//			break;
	//		}
	//		slab_initial = true;
	//	}
	//	slab_mesh.initCollapseQueue();
	//	long endt = clock();
	//	return endt - startt;
	//}
}

void ThreeDimensionalShape::InitialSlabMesh()
{
	// handle each face
	for(unsigned i = 0; i < slab_mesh.vertices.size(); i++)
	{ 
		if(!slab_mesh.vertices[i].first)
			continue;

		SlabVertex sv = *slab_mesh.vertices[i].second;
		std::set<unsigned> fset = sv.faces_;
		Vector4d C1(sv.sphere.center.X(), sv.sphere.center.Y(), sv.sphere.center.Z(), sv.sphere.radius);

		for (set<unsigned>::iterator si = fset.begin(); si != fset.end(); si++)
		{
			SlabFace sf = *slab_mesh.faces[*si].second;

			if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) || 
				sf.st[1].normal == Vector3d(0., 0., 0.))
				continue;

			Vector4d normal1(sf.st[0].normal.X(), sf.st[0].normal.Y(), sf.st[0].normal.Z(), 1.0);
			Vector4d normal2(sf.st[1].normal.X(), sf.st[1].normal.Y(), sf.st[1].normal.Z(), 1.0);

			// compute the matrix of A
			Matrix4d temp_A1, temp_A2;
			temp_A1.MakeTensorProduct(normal1, normal1);
			temp_A2.MakeTensorProduct(normal2, normal2);
			temp_A1 *= 2.0;
			temp_A2 *= 2.0;

			// compute the matrix of b
			double normal_mul_point1 = normal1.Dot(C1);
			double normal_mul_point2 = normal2.Dot(C1);
			Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
			Wm4::Vector4d temp_b2 = normal2 * 2 * normal_mul_point2;

			//compute c
			double temp_c1 = normal_mul_point1 * normal_mul_point1;
			double temp_c2 = normal_mul_point2 * normal_mul_point2;

			slab_mesh.vertices[i].second->slab_A += temp_A1;
			slab_mesh.vertices[i].second->slab_A += temp_A2;
			slab_mesh.vertices[i].second->slab_b += temp_b1;
			slab_mesh.vertices[i].second->slab_b += temp_b2;
			slab_mesh.vertices[i].second->slab_c += temp_c1;
			slab_mesh.vertices[i].second->slab_c += temp_c2;

			slab_mesh.vertices[i].second->related_face += 2;
		}
	}

	switch(slab_mesh.preserve_boundary_method)
	{
	case 1 :
		slab_mesh.PreservBoundaryMethodOne();
		break;
	case 2 :
		//slab_mesh.PreservBoundaryMethodTwo();
		break;
	case 3 :
		slab_mesh.PreservBoundaryMethodThree();
		break;
	default:
		slab_mesh.PreservBoundaryMethodFour();
		break;
	}

}

double ThreeDimensionalShape::NearestPoint(Vector3d point, unsigned vid)
{
	set<unsigned> faces = slab_mesh.vertices[vid].second->faces_;
	set<unsigned> edges = slab_mesh.vertices[vid].second->edges_;
	double mind = DBL_MAX;

	// calculation of related faces
	for (set<unsigned>::iterator si = faces.begin(); si != faces.end(); si++)
	{
		if (!slab_mesh.faces[*si].first)
			continue;
		SlabFace sf = *slab_mesh.faces[*si].second;
		if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) || 
			sf.st[1].normal == Vector3d(0., 0., 0.))
			continue;

		Vector3d v[3];
		Vector3d tfp;
		double td;
		for (int i = 0; i < 2; i++)
		{
			v[0] = sf.st[i].v[0];
			v[1] = sf.st[i].v[1];
			v[2] = sf.st[i].v[2];
			ProjectOntoTriangle(point, v[0], v[1], v[2], tfp, td);

			if(td < mind)	mind = td;
		}
	}

	// calculation of related edges
	for (set<unsigned>::iterator si = edges.begin(); si != edges.end(); si++)
	{
		if (!slab_mesh.edges[*si].first)
			continue;
		SlabEdge se = *slab_mesh.edges[*si].second;
		if (se.valid_cone == false)
			continue;

		Vector3d tfp;
		double td;
		se.cone.ProjectOntoCone(point, tfp, td);
		if (td < 0)
			td = -td;

		if(td < mind)	mind = td;
	}

	return mind;
}

void ThreeDimensionalShape::ComputeHausdorffDistance()
{
#if 0
	long start_time = clock();

	ma_qem_mesh.maxhausdorff_distance = 0.;
	for(unsigned i = 0; i < ma_qem_mesh.faces.size(); i ++)
	{
		Vector3d ve[8];
		if (ma_qem_mesh.faces[i].second->valid_st == false ||
			ma_qem_mesh.faces[i].second->st[0].normal == Vector3d(0., 0., 0.) || 
			ma_qem_mesh.faces[i].second->st[1].normal == Vector3d(0., 0., 0.))
			continue;

		ve[0] = ma_qem_mesh.faces[i].second->st[0].v[0];
		ve[1] = ma_qem_mesh.faces[i].second->st[0].v[1];
		ve[2] = ma_qem_mesh.faces[i].second->st[0].v[2];
		ve[3] = ma_qem_mesh.faces[i].second->st[1].v[0];
		ve[4] = ma_qem_mesh.faces[i].second->st[1].v[1];
		ve[5] = ma_qem_mesh.faces[i].second->st[1].v[2];
		ve[6] = (ve[0] + ve[1] +ve[2]) / 3.0;
		ve[7] = (ve[3] + ve[4] +ve[5]) / 3.0;
		double face_haus = 0.;
		for (int j = 0; j < 8; j++)
		{
			Vector3d fp = input.NearestVertex(ve[j]);
			double len = (ve[j] - fp).Length();
			face_haus = max(len, face_haus);
		}
		ma_qem_mesh.faces[i].second->hausdorff_dist = face_haus;

		ma_qem_mesh.maxhausdorff_distance = max(ma_qem_mesh.maxhausdorff_distance,face_haus);
	}

	long end_time = clock();
	long result = end_time - start_time;
#endif

	//ma_qem_mesh.maxhausdorff_distance = 0.;
	slab_mesh.maxhausdorff_distance = 0;
	double sumhausdorff_distance = 0;
	for (unsigned i = 0; i < input.pVertexList.size(); i++)
	{
		double min_dis = DBL_MAX;
		unsigned min_index = -1;
		Vector3d bou_ver(input.pVertexList[i]->point()[0], input.pVertexList[i]->point()[1], input.pVertexList[i]->point()[2]);
		bou_ver /=  input.bb_diagonal_length; 

		for (unsigned j = 0; j < slab_mesh.numVertices; j++)
		{
			Sphere ma_ver = slab_mesh.vertices[j].second->sphere;
			double temp_length = abs((bou_ver - ma_ver.center).Length() - ma_ver.radius);
			//if (temp_length >= 0 && temp_length < min_dis)
			if (temp_length < min_dis)
			{
				min_dis = temp_length;
				min_index = j;
			}

			//double temp_near_dis = slab_mesh.NearestPoint(bou_ver, min_index);
			//if (temp_near_dis < min_dis)
			//{
			//	min_dis = temp_near_dis;
			//	min_index = j;
			//}

		}

		//// 为何这里得出的结果比前面得出的结果还要小？
		//double nearest_dis = NearestPoint(bou_ver, min_index);
		//min_dis = min(nearest_dis, min_dis);

		//sumhausdorff_distance += min_dis;


		if (min_index != -1)
		{	
			double temp_near_dis = slab_mesh.NearestPoint(bou_ver, min_index);
			min_dis = min(temp_near_dis, min_dis);

			sumhausdorff_distance += min_dis;

			//ma_qem_mesh.vertices[min_index].second->bplist.push_back(i);
			//ma_qem_mesh.maxhausdorff_distance = max(ma_qem_mesh.maxhausdorff_distance, min_dis);

			slab_mesh.vertices[min_index].second->bplist.insert(i);
			slab_mesh.maxhausdorff_distance = max(slab_mesh.maxhausdorff_distance, min_dis);

			//input.pVertexList[i]->vqem_hausdorff_dist = min_dis / input.bb_diagonal_length;
			input.pVertexList[i]->vqem_hausdorff_dist = min_dis;
			input.pVertexList[i]->vqem_hansdorff_index = min_index;

			//input.pVertexList[i]->slab_hausdorff_dist = min_dis / input.bb_diagonal_length;
			input.pVertexList[i]->slab_hausdorff_dist = min_dis;
			input.pVertexList[i]->slab_hansdorff_index = min_index;
		}
	}
	
	//ma_qem_mesh.meanhausdorff_distance = sumhausdorff_distance / input.pVertexList.size();
	slab_mesh.meanhausdorff_distance = sumhausdorff_distance / input.pVertexList.size();
	//ma_qem_mesh.initialhausdorff_distance = ma_qem_mesh.maxhausdorff_distance;
	slab_mesh.initialhausdorff_distance = slab_mesh.maxhausdorff_distance;
}

void ThreeDimensionalShape::PruningSlabMesh()
{
	slab_mesh.ComputeVerticesProperty();

	bool has_boundary_non_pole;
	do
	{
		has_boundary_non_pole = false;
		unsigned vid;
		for(unsigned i = 0; i < slab_mesh.vertices.size(); i ++)
			if(slab_mesh.vertices[i].first)
				if((slab_mesh.vertices[i].second->is_boundary) &&
					(!slab_mesh.vertices[i].second->is_non_manifold) &&
					(!slab_mesh.vertices[i].second->is_pole) &&
					(slab_mesh.vertices[i].second->edges_.size() == 2))
				{
					vid = i;
					has_boundary_non_pole = true;
					break;
				}
				if(has_boundary_non_pole)
					slab_mesh.DeleteVertex(vid);
	}while(has_boundary_non_pole);


	slab_mesh.CleanIsolatedVertices();
	slab_mesh.ComputeEdgesCone();
	slab_mesh.ComputeFacesSimpleTriangles();
	slab_mesh.DistinguishVertexType();
	slab_mesh.computebb();
}