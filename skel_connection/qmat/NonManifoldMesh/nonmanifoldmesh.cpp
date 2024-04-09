#include "nonmanifoldmesh.h"
#include <ctime>
#include <cstdio>
#include <boost/lexical_cast.hpp>
#include <string>


void NonManifoldMesh::AdjustStorage()
{
	std::vector<unsigned> newv;
	std::vector<unsigned> newe;
	std::vector<unsigned> newf;

	newv.resize(vertices.size());
	newe.resize(edges.size());
	newf.resize(faces.size());

	unsigned count = 0;
	for(unsigned i = 0; i < vertices.size(); i ++)
		if(vertices[i].first)
			newv[i] = count ++;

	count = 0;
	for(unsigned i = 0; i < edges.size(); i ++)
		if(edges[i].first)
			newe[i] = count ++;

	count = 0;
	for(unsigned i = 0; i < faces.size(); i ++)
		if(faces[i].first)
			newf[i] = count ++;

	std::vector<Bool_VertexPointer> new_vertices;
	std::vector<Bool_EdgePointer> new_edges;
	std::vector<Bool_FacePointer> new_faces;

	for(unsigned i = 0; i < vertices.size(); i ++)
		if(vertices[i].first)
		{
			Bool_VertexPointer bvp;
			bvp = vertices[i];
			std::set<unsigned> neweset;
			std::set<unsigned> newfset;

			for(std::set<unsigned>::iterator si = bvp.second->edges_.begin();
				si != bvp.second->edges_.end(); si ++)
				neweset.insert(newe[*si]);
			for(std::set<unsigned>::iterator si = bvp.second->faces_.begin();
				si != bvp.second->faces_.end(); si ++)
				newfset.insert(newf[*si]);

			bvp.second->edges_ = neweset;
			bvp.second->faces_ = newfset;
			new_vertices.push_back(bvp);
		}

		for(unsigned i = 0; i < edges.size(); i ++)
			if(edges[i].first)
			{
				Bool_EdgePointer bep;
				bep = edges[i];
				std::pair<unsigned,unsigned> newvpair;
				std::set<unsigned> newfset;

				newvpair.first = newv[bep.second->vertices_.first];
				newvpair.second = newv[bep.second->vertices_.second];

				for(std::set<unsigned>::iterator si = bep.second->faces_.begin();
					si != bep.second->faces_.end(); si ++)
					newfset.insert(newf[*si]);

				bep.second->vertices_ = newvpair;
				bep.second->faces_ = newfset;
				new_edges.push_back(bep);
			}

			for(unsigned i = 0; i < faces.size(); i ++)
				if(faces[i].first)
				{
					Bool_FacePointer bfp;
					bfp = faces[i];
					std::set<unsigned> newvset;
					std::set<unsigned> neweset;

					for(std::set<unsigned>::iterator si = bfp.second->vertices_.begin();
						si != bfp.second->vertices_.end(); si ++)
						newvset.insert(newv[*si]);

					for(std::set<unsigned>::iterator si = bfp.second->edges_.begin();
						si != bfp.second->edges_.end(); si ++)
						neweset.insert(newe[*si]);

					bfp.second->vertices_ = newvset;
					bfp.second->edges_ = neweset;
					new_faces.push_back(bfp);
				}

				vertices = new_vertices;
				edges = new_edges;
				faces = new_faces;
}

bool NonManifoldMesh::ValidVertex(unsigned vid)
{
	if(vid > vertices.size())
		return false;

	return vertices[vid].first;
}

bool NonManifoldMesh::Edge(unsigned vid0, unsigned vid1, unsigned & eid)
{
	if(!ValidVertex(vid0) || !ValidVertex(vid1))
		return false;

	for(std::set<unsigned>::iterator si = (*vertices[vid0].second).edges_.begin(); si != (*vertices[vid0].second).edges_.end(); si ++)
	{
		if(edges[*si].first)
		{
			if(edges[*si].second->HasVertex(vid1))
			{
				eid = *si;
				return true;
			}
		}
	}

	return false;
}

bool NonManifoldMesh::Face(const std::set<unsigned> & vset, unsigned & fid)
{
	if(vset.size() <= 0)
		return false;

	for(std::set<unsigned>::iterator si = vset.begin(); si != vset.end(); si ++)
		if(!ValidVertex(*si))
			return false;

	unsigned vid0 = *(vset.begin());
	for(std::set<unsigned>::iterator si = vertices[vid0].second->faces_.begin();
		si != vertices[vid0].second->faces_.end(); si ++)
		if(faces[*si].first)
		{
			if(faces[*si].second->vertices_ == vset)
			{
				fid = *si;
				return true;
			}
		}

	return false;
}

void NonManifoldMesh::UpdateCentroid(unsigned fid)
{
	if(!faces[fid].first)
		return;

	faces[fid].second->centroid = Wm4::Vector3d::ZERO;
	unsigned count(0);
	for(std::set<unsigned>::iterator si = faces[fid].second->vertices_.begin(); si != faces[fid].second->vertices_.end(); si ++, count ++)
	{
		faces[fid].second->centroid += vertices[*si].second->sphere.center;
	}
	faces[fid].second->centroid /= count;
}

void NonManifoldMesh::ComputeFacesCentroid()
{
	for(unsigned i = 0; i < faces.size(); i ++)
		if(faces[i].first)
			UpdateCentroid(i);
}

void NonManifoldMesh::UpdateNormal(unsigned fid)
{
	if(!faces[fid].first)
		return;

	Vector3d v[3];
	std::set<unsigned>::iterator si = faces[fid].second->vertices_.begin();
	v[0] = vertices[*si].second->sphere.center;
	si ++;
	v[1] = vertices[*si].second->sphere.center;
	si ++;
	v[2] = vertices[*si].second->sphere.center;
	faces[fid].second->normal = (v[1]-v[0]).Cross(v[2]-v[0]);
	faces[fid].second->normal.Normalize();
}

void NonManifoldMesh::ComputeFacesNormal()
{
	for(unsigned i = 0; i < faces.size(); i ++)
		if(faces[i].first)
			UpdateNormal(i);
}

void NonManifoldMesh::GetNeighborVertices(unsigned vid, std::set<unsigned> & neighborvertices)
{
	if(!vertices[vid].first)
		return;
	for(std::set<unsigned>::iterator si = vertices[vid].second->edges_.begin();
		si != vertices[vid].second->edges_.end(); si ++)
	{
		if(edges[*si].first)
		{
			if(edges[*si].second->vertices_.first == vid)
				neighborvertices.insert(edges[*si].second->vertices_.second);
			else
				neighborvertices.insert(edges[*si].second->vertices_.first);
		}
	}
}

void NonManifoldMesh::GetLinkedEdges(unsigned eid, std::set<unsigned> & neighboredges)
{
	if(!edges[eid].first)
		return;
	neighboredges.clear();
	unsigned vid[2];
	vid[0] = edges[eid].second->vertices_.first;
	vid[1] = edges[eid].second->vertices_.second;
	for(unsigned k = 0; k < 2; k ++)
	{
		if(vertices[vid[k]].first)
		{
			for(std::set<unsigned>::iterator si = vertices[vid[k]].second->edges_.begin();
				si != vertices[vid[k]].second->edges_.end(); si ++)
				if(edges[*si].first)
					neighboredges.insert(*si);
		}
	}
	neighboredges.erase(eid);
}

void NonManifoldMesh::GetAdjacentFaces(unsigned fid, std::set<unsigned> & neighborfaces)
{
	if(!faces[fid].first)
		return;

	neighborfaces.clear();
	for(std::set<unsigned>::iterator si = faces[fid].second->edges_.begin();
		si != faces[fid].second->edges_.end(); si ++)
	{
		if(edges[*si].first)
		{
			for(std::set<unsigned>::iterator si2 = edges[*si].second->faces_.begin();
				si2 != edges[*si].second->faces_.end(); si2 ++)
				neighborfaces.insert(*si2);
		}
	}
	neighborfaces.erase(fid);
}

bool NonManifoldMesh::Contractible(unsigned vid_src, unsigned vid_tgt)
{
	if( !vertices[vid_src].first || !vertices[vid_tgt].first )
		return false;

	set<unsigned> ns;
	set<unsigned> nt;
	set<unsigned> ni;
	GetNeighborVertices(vid_src, ns);
	GetNeighborVertices(vid_tgt, nt);
	for(set<unsigned>::iterator si = ns.begin(); si != ns.end(); si ++)
		if(nt.find(*si) != nt.end())
			ni.insert(*si);
	unsigned eid;
	bool foundedge = Edge(vid_src, vid_tgt, eid);
	if(!foundedge)
		return false;
	if(edges[eid].second->faces_.size() != ni.size())
		return false;

	for(std::set<unsigned>::iterator si = vertices[vid_src].second->faces_.begin();
		si != vertices[vid_src].second->faces_.end(); si ++)
	{
		if(faces[*si].first)
		{
			if(!faces[*si].second->HasVertex(vid_tgt))
			{
				Vector3d pp[3], pa[3];
				unsigned count = 0;
				for(std::set<unsigned>::iterator si2 = faces[*si].second->vertices_.begin();
					si2 != faces[*si].second->vertices_.end(); si2 ++)
				{
					pp[count] = vertices[*si2].second->sphere.center;
					pa[count++] = (*si2 != vid_src)?vertices[*si2].second->sphere.center:vertices[vid_tgt].second->sphere.center;
				}
				Vector3d pnorm = TriangleNormal(pp[0],pp[1],pp[2]);
				Vector3d anorm = TriangleNormal(pa[0],pa[1],pa[2]);
				if(pnorm.Dot(anorm) < 0)
					return false;
			}
		}
	}
	return true;
}

bool NonManifoldMesh::MergeVertices(unsigned vid_src1, unsigned vid_src2, unsigned &vid_tgt)
{
	if(vid_src1 == vid_src2)
		return false;
	
	unsigned eid;
	if(!Edge(vid_src1, vid_src2, eid))
		return false;

	InsertVertex(new NonManifoldMesh_Vertex, vid_tgt);

	std::vector< std::set<unsigned> > tri_vec;
	for(std::set<unsigned>::iterator si = vertices[vid_src1].second->faces_.begin();
		si != vertices[vid_src1].second->faces_.end(); si ++)
		if(!faces[*si].second->HasVertex(vid_tgt))
		{
			std::set<unsigned> vset = faces[*si].second->vertices_;
			vset.erase(vid_src1);
			vset.insert(vid_tgt);
			tri_vec.push_back(vset);
		}

	for(std::set<unsigned>::iterator si = vertices[vid_src2].second->faces_.begin();
		si != vertices[vid_src2].second->faces_.end(); si ++)
		if(!faces[*si].second->HasVertex(vid_tgt))
		{
			std::set<unsigned> vset = faces[*si].second->vertices_;
			vset.erase(vid_src2);
			vset.insert(vid_tgt);
			tri_vec.push_back(vset);
		}

	std::vector< std::pair<unsigned,unsigned> > edge_vec;
	for(std::set<unsigned>::iterator si = vertices[vid_src1].second->edges_.begin();
		si != vertices[vid_src1].second->edges_.end(); si ++)
		if(!edges[*si].second->HasVertex(vid_tgt))
		{
			std::pair<unsigned, unsigned> vp = edges[*si].second->vertices_;
			if(vp.first == vid_src1)
				vp.first = vid_tgt;
			if(vp.second == vid_src1)
				vp.second = vid_tgt;
			edge_vec.push_back(vp);
		}

	for(std::set<unsigned>::iterator si = vertices[vid_src2].second->edges_.begin();
		si != vertices[vid_src2].second->edges_.end(); si ++)
		if(!edges[*si].second->HasVertex(vid_tgt))
		{
			std::pair<unsigned, unsigned> vp = edges[*si].second->vertices_;
			if(vp.first == vid_src2)
				vp.first = vid_tgt;
			if(vp.second == vid_src2)
				vp.second = vid_tgt;
			edge_vec.push_back(vp);
		}

	DeleteVertex(vid_src1);
	DeleteVertex(vid_src2);

	for(unsigned i = 0; i < tri_vec.size(); i ++)
		InsertFace(tri_vec[i]);

	for(unsigned i = 0; i < edge_vec.size(); i ++)
	{
		unsigned neweid;
		InsertEdge(edge_vec[i].first, edge_vec[i].second, neweid);
	}

	return true;
		
}

unsigned NonManifoldMesh::VertexIncidentEdgeCount(unsigned vid)
{
	if(!vertices[vid].first)
		return 0;
	return (unsigned)vertices[vid].second->edges_.size();
}

unsigned NonManifoldMesh::VertexIncidentFaceCount(unsigned vid)
{
	if(!vertices[vid].first)
		return 0;
	return (unsigned)vertices[vid].second->faces_.size();
}

unsigned NonManifoldMesh::EdgeIncidentFaceCount(unsigned eid)
{
	if(!edges[eid].first)
		return 0;
	return (unsigned)edges[eid].second->faces_.size();
}

void NonManifoldMesh::DeleteFace(unsigned fid)
{
	if(!faces[fid].first)
		return;

	for(std::set<unsigned>::iterator si = faces[fid].second->vertices_.begin();
		si != faces[fid].second->vertices_.end(); si ++)
			vertices[*si].second->faces_.erase(fid);

	for(std::set<unsigned>::iterator si = faces[fid].second->edges_.begin();
		si != faces[fid].second->edges_.end(); si ++)
			edges[*si].second->faces_.erase(fid);

	delete faces[fid].second;
	faces[fid].first = false;
	numFaces --;
}

void NonManifoldMesh::DeleteEdge(unsigned eid)
{
	if(!edges[eid].first)
		return;

	if(vertices[edges[eid].second->vertices_.first].first)
		vertices[edges[eid].second->vertices_.first].second->edges_.erase(eid);
	if(vertices[edges[eid].second->vertices_.second].first)
		vertices[edges[eid].second->vertices_.second].second->edges_.erase(eid);
	/*
	for(std::set<unsigned>::iterator si = edges[eid].second->faces_.begin();
		si != edges[eid].second->faces_.end(); si ++)
		if(faces[*si].first)
			faces[*si].second->edges_.erase(eid);
			*/
	std::set<unsigned> faces_del;
	for(std::set<unsigned>::iterator si = edges[eid].second->faces_.begin();
		si != edges[eid].second->faces_.end(); si ++)
			faces_del.insert(*si);
	for(std::set<unsigned>::iterator si = faces_del.begin(); si != faces_del.end(); si ++)
		DeleteFace(*si);
			//faces[*si].second->edges_.erase(eid);

	delete edges[eid].second;
	edges[eid].first = false;
	numEdges --;
}

void NonManifoldMesh::DeleteVertex(unsigned vid)
{
	if(!vertices[vid].first)
		return;

	std::set<unsigned> edges_del;
	for(std::set<unsigned>::iterator si = vertices[vid].second->edges_.begin();
		si != vertices[vid].second->edges_.end(); si ++)
		edges_del.insert(*si);

	std::set<unsigned> faces_del;
	for(std::set<unsigned>::iterator si = vertices[vid].second->faces_.begin();
		si != vertices[vid].second->faces_.end(); si ++)
		faces_del.insert(*si);

	for(std::set<unsigned>::iterator si = edges_del.begin(); si != edges_del.end(); si ++)
		DeleteEdge(*si);

	for(std::set<unsigned>::iterator si = faces_del.begin(); si != faces_del.end(); si ++)
		DeleteFace(*si);

	delete vertices[vid].second;
	vertices[vid].first = false;
	numVertices --;
}

void NonManifoldMesh::InsertVertex(NonManifoldMesh_Vertex *vertex, unsigned &vid){
	Bool_VertexPointer bvp;
	bvp.first = true;
	bvp.second = vertex;
	vid = (unsigned)vertices.size();
	vertices.push_back(bvp);
	numVertices ++;
}

void NonManifoldMesh::InsertEdge(unsigned vid0, unsigned vid1, unsigned & eid)
{
	if(Edge(vid0,vid1,eid))
		return;
	if (!ValidVertex(vid0) || !ValidVertex(vid1))
	{
		return;
	}
	Bool_EdgePointer bep;
	bep.first = true;
	bep.second = new NonManifoldMesh_Edge;
	bep.second->vertices_.first = vid0;
	bep.second->vertices_.second = vid1;
	vertices[vid0].second->edges_.insert((unsigned)edges.size());
	vertices[vid1].second->edges_.insert((unsigned)edges.size());
	eid = (unsigned)edges.size();
	edges.push_back(bep);
	//ComputeEdgeCone(eid);
	numEdges ++;
}

void NonManifoldMesh::InsertFace(std::set<unsigned> vset)
{
	unsigned fid;
	if(Face(vset,fid))
		return;

	unsigned vid[3];
	std::set<unsigned>::iterator si = vset.begin();
	vid[0] = *si;
	si ++;
	vid[1] = *si;
	si ++;
	vid[2] = *si;

	if(!vertices[vid[0]].first || !vertices[vid[1]].first || !vertices[vid[2]].first)
		return;

	for(std::set<unsigned>::iterator si = vertices[vid[0]].second->faces_.begin(); 
		si != vertices[vid[0]].second->faces_.end(); si ++)
		if(faces[*si].second->vertices_ == vset) // duplicate
			return;

	unsigned eid[3];
	InsertEdge(vid[0],vid[1],eid[0]);
	InsertEdge(vid[0],vid[2],eid[1]);
	InsertEdge(vid[1],vid[2],eid[2]);

	Bool_FacePointer bfp;
	bfp.first = true;
	bfp.second = new NonManifoldMesh_Face;
	bfp.second->vertices_.insert(vid[0]);
	bfp.second->vertices_.insert(vid[1]);
	bfp.second->vertices_.insert(vid[2]);
	bfp.second->edges_.insert(eid[0]);
	bfp.second->edges_.insert(eid[1]);
	bfp.second->edges_.insert(eid[2]);
	vertices[vid[0]].second->faces_.insert(faces.size());
	vertices[vid[1]].second->faces_.insert(faces.size());
	vertices[vid[2]].second->faces_.insert(faces.size());
	edges[eid[0]].second->faces_.insert(faces.size());
	edges[eid[1]].second->faces_.insert(faces.size());
	edges[eid[2]].second->faces_.insert(faces.size());
	faces.push_back(bfp);
	//UpdateCentroid((unsigned)faces.size()-1);
	//UpdateNormal((unsigned)faces.size()-1);
	ComputeFaceSimpleTriangles((unsigned)faces.size() - 1);
	numFaces ++;
}

void NonManifoldMesh::ComputeEdgeCone(unsigned eid)
{
	if(!edges[eid].first)
		return;
	Cone newc(vertices[edges[eid].second->vertices_.first].second->sphere.center, vertices[edges[eid].second->vertices_.first].second->sphere.radius,
		vertices[edges[eid].second->vertices_.second].second->sphere.center, vertices[edges[eid].second->vertices_.second].second->sphere.radius);
	edges[eid].second->cone = newc;
	if(newc.type == 1)
		edges[eid].second->valid_cone = false;
	else
		edges[eid].second->valid_cone = true;
}

void NonManifoldMesh::ComputeEdgesCone()
{
	for(unsigned i = 0; i < edges.size(); i ++)
		if(edges[i].first)
			ComputeEdgeCone(i);
}

void NonManifoldMesh::ComputeFaceSimpleTriangles(unsigned fid)
{
	if(!faces[fid].first)
		return;
	SimpleTriangle st0,st1;
	Wm4::Vector3d pos[3];
	double radius[3];
	unsigned count = 0;
	for(std::set<unsigned>::iterator si = faces[fid].second->vertices_.begin();
		si != faces[fid].second->vertices_.end(); si ++, count ++)
	{
		pos[count] = vertices[*si].second->sphere.center;
		radius[count] = vertices[*si].second->sphere.radius;
	}
	if(TriangleFromThreeSpheres(pos[0],radius[0],pos[1],radius[1],pos[2],radius[2],st0,st1))
	{
		faces[fid].second->st[0] = st0;
		faces[fid].second->st[1] = st1;
		faces[fid].second->valid_st = true;
	}
	else
		faces[fid].second->valid_st = false;

}

void NonManifoldMesh::ComputeFacesSimpleTriangles()
{
	for(unsigned i = 0; i < faces.size(); i ++)
		if(faces[i].first)
			ComputeFaceSimpleTriangles(i);
}

void NonManifoldMesh::MinCostEdgeCollapse(unsigned & eid){
	//merge 2 vertices of the edge first, then move the combined vertex to the preferred point and resize it.
	unsigned v1, v2;
	v1 = edges[eid].second->vertices_.first;
	v2 = edges[eid].second->vertices_.second;
	Wm4::Matrix4d Q = edges[eid].second->Q;
	Sphere sphere = edges[eid].second->sphere;
	if (!Contractible(v1, v2))
		return;

	unsigned vid_tgt;
	if(MergeVertices(v1, v2, vid_tgt)){
		vertices[vid_tgt].second->Q = Q;
		vertices[vid_tgt].second->sphere = sphere;

		for (std::set<unsigned>::iterator si = vertices[vid_tgt].second->edges_.begin(); si != vertices[vid_tgt].second->edges_.end(); si ++)
		{
			EvaluateEdgeCollapseCost(*si);
			edge_collapses_queue.push(*(edges[*si].second));
			ComputeEdgeCone(*si);
		}
	}

}

void NonManifoldMesh::Simplify(int threshold ){
	int deleteNum = 0;
	while (deleteNum < threshold && vertices.size() > 0)
	{
		NonManifoldMesh_Edge topEdge = edge_collapses_queue.top();
		edge_collapses_queue.pop();
		unsigned eid;
		if(Edge(topEdge.vertices_.first, topEdge.vertices_.second, eid)){
			MinCostEdgeCollapse(eid);
			deleteNum ++;
		}
	}
}

void NonManifoldMesh::EvaluateEdgeCollapseCost(unsigned eid){
	if (!edges[eid].first)
		return ;

	unsigned v1, v2;
	v1 = edges[eid].second->vertices_.first;
	v2 = edges[eid].second->vertices_.second;

	Wm4::Matrix4d Q1 = vertices[v1].second->Q;
	Wm4::Matrix4d Q2 = vertices[v2].second->Q;

	// get the location for min
	Wm4::Matrix4d temp_matrix = Q1 + Q2;
	temp_matrix.SetRow(3, Vector4d(0, 0, 0, 1));
	Wm4::Vector4d min_vertex = temp_matrix.Inverse() * Vector4d(0, 0, 0, 1);

	double coll_cost = temp_matrix.QForm(min_vertex, min_vertex);

	edges[eid].second->Q = Q1 + Q2;
	edges[eid].second->collapse_cost = coll_cost;
	edges[eid].second->sphere.center = Wm4::Vector3d(min_vertex.X(), min_vertex.Y(), min_vertex.Z());
	edges[eid].second->sphere.radius = min_vertex.W();

}

void NonManifoldMesh::initCollapseQueue(){
	for (int i = 0; i < numEdges; i ++)
	{
		if (edges[i].first)
		{
			EvaluateEdgeCollapseCost(i);
			edge_collapses_queue.push(*edges[i].second);
		}
	}
}

void NonManifoldMesh::Export(std::string fname){

	if (fname.find(".off") == std::string::npos)
	{
		fname += "___v_";
		fname += std::to_string(static_cast<long long>(numVertices));
		fname += "___e_";
		fname += std::to_string(static_cast<long long>(numEdges));
		fname += "___f_";
		fname += std::to_string(static_cast<long long>(numFaces));
	}
	else
	{
		fname = fname.substr(0, fname.find(".off"));
	}

	AdjustStorage();

	//std::string sphname = fname;
	//sphname += ".sph";

	//std::ofstream fsphout(sphname);

	//std::vector<Sphere> vecs;
	//for(unsigned i = 0; i < vertices.size(); i ++)
	//	vecs.push_back(vertices[i].second->sphere);
	//
	//for(unsigned i = 0; i < edges.size(); i ++)
	//{
	//	if(edges[i].second->validenvelope == false)
	//		continue;
	//	std::vector<Sphere> evs;
	//	evs = edges[i].second->cone.SampleSpheres(10);
	//	for(unsigned k = 0; k < evs.size(); k ++)
	//		vecs.push_back(evs[k]);
	//}
	//
	//for(unsigned i = 0; i < faces.size(); i ++)
	//{
	//	if(faces[i].second->validenvelope == false)
	//		continue;
	//	std::vector<Sphere> fvs = SampleSpheres(i,4);
	//	for(unsigned k = 0; k < fvs.size(); k ++)
	//		vecs.push_back(fvs[k]);
	//}
	///**/
	//fsphout << vecs.size() << std::endl;
	//for(unsigned i = 0; i < vecs.size(); i ++)
	//	fsphout << vecs[i].center << " " << vecs[i].radius << std::endl;
	//fsphout.close();



	std::string maname = fname;
	maname += ".ma";

	//std::ofstream fout("3dma.ma");
	std::ofstream fout(maname);

//	GraphVertexIterator gvi,gvi_end;

	fout << numVertices << " " << numEdges << " " << numFaces << std::endl;
	
	//fout << num_vertices(*g) << " " << num_edges(*g) << " " << g->tris.size() << std::endl;

	for(unsigned i = 0; i < vertices.size(); i ++)
		//fout << "v " << vertices[i].second->sphere.center << " " << vertices[i].second->sphere.radius << std::endl;
		fout << "v " << setiosflags(ios::fixed) << setprecision(15) << vertices[i].second->sphere.center << " " << vertices[i].second->sphere.radius << std::endl;

	for(unsigned i = 0; i < edges.size(); i ++)
		fout << "e " << edges[i].second->vertices_.first << " " << edges[i].second->vertices_.second << std::endl;
	for(unsigned i = 0; i < faces.size(); i ++)
	{
		fout << "f";
		for(std::set<unsigned>::iterator si = faces[i].second->vertices_.begin();
			si != faces[i].second->vertices_.end(); si ++)
			fout << " " << *si;
		fout << std::endl;
	}
	/*
	for(boost::tie(gvi,gvi_end) = vertices(*g); gvi != gvi_end; gvi ++)
		fout << "v "<< (*g)[*gvi].pos[0] << ' ' << (*g)[*gvi].pos[1] << ' ' << (*g)[*gvi].pos[2] << ' ' << (*g)[*gvi].radius << std::endl;
	for(std::pair<GraphEdgeIterator, GraphEdgeIterator> geip = edges(*g); geip.first != geip.second; geip.first ++)
		fout << "e " << boost::source(*geip.first, *g) << " " << boost::target(*geip.first, *g) << std::endl;// << " " << boost::target(*geip.first, *g) + 1 << std::endl;
	for(unsigned int i = 0; i < g->tris.size(); i ++)
		fout << "f " << MappingIdtoGVD(g,g->tris[i].vid[0]) << ' ' << MappingIdtoGVD(g,g->tris[i].vid[1]) << ' ' << MappingIdtoGVD(g,g->tris[i].vid[2]) << std::endl;
	*/
	fout.close();

	//std::string mappingname = fname;
	//mappingname += ".mapping";

	////std::ofstream fmapping("3dma.mapping");
	//std::ofstream fmapping(mappingname);
	//std::vector< std::set<unsigned int> > vs;
	//vs.resize(BoundaryPoints.size());
	////vs.resize(g->boundarysamplepoints.size());
	//fmapping << vs.size() << std::endl;
	//for(unsigned i = 0; i < vertices.size(); i ++)
	//	for(std::set<unsigned>::iterator si = vertices[i].second->bplist.begin();
	//		si != vertices[i].second->bplist.end(); si ++)
	//		vs[*si].insert(i);
	///*
	//for(boost::tie(gvi,gvi_end) = vertices(*g); gvi != gvi_end; gvi ++)
	//	for(unsigned int i = 0 ; i < (*g)[*gvi].v_samples.size(); i ++)
	//		vs[(*g)[*gvi].v_samples[i]].insert(*gvi);
	//		*/

	//for(unsigned int i = 0; i < vs.size(); i ++)
	//{
	//	fmapping << vs[i].size() << " ";
	//	for(std::set<unsigned int>::iterator si = vs[i].begin(); si != vs[i].end(); si ++)
	//		fmapping << *si << " " << 1. / vs[i].size() << " ";
	//	fmapping << std::endl;
	//}
	//	

	//fmapping.close();
}