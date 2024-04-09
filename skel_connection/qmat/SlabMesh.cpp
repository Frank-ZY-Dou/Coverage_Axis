#include "SlabMesh.h"
#include <omp.h>

void SlabMesh::AdjustStorage()
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

	std::vector<Bool_SlabVertexPointer> new_vertices;
	std::vector<Bool_SlabEdgePointer> new_edges;
	std::vector<Bool_SlabFacePointer> new_faces;

	for(unsigned i = 0; i < vertices.size(); i ++)
		if(vertices[i].first)
		{
			Bool_SlabVertexPointer bvp;
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
				Bool_SlabEdgePointer bep;
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
					Bool_SlabFacePointer bfp;
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

bool SlabMesh::ValidVertex(unsigned vid){
	if(vid > vertices.size())
		return false;

	return vertices[vid].first;
}

bool SlabMesh::Edge(unsigned vid0, unsigned vid1, unsigned & eid)
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

bool SlabMesh::Face(const std::set<unsigned> & vset, unsigned & fid)
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

void SlabMesh::UpdateCentroid(unsigned fid)
{
	if(!faces[fid].first)
		return;

	faces[fid].second->centroid.center = Wm4::Vector3d::ZERO;
	faces[fid].second->centroid.radius = 0.0;

	unsigned count(0);
	for(std::set<unsigned>::iterator si = faces[fid].second->vertices_.begin(); si != faces[fid].second->vertices_.end(); si ++, count ++)
	{
		faces[fid].second->centroid.center += vertices[*si].second->sphere.center;
		faces[fid].second->centroid.radius += vertices[*si].second->sphere.radius;
	}
	faces[fid].second->centroid.center /= count;
	faces[fid].second->centroid.radius /= count;
}

void SlabMesh::ComputeFacesCentroid()
{
	for(unsigned i = 0; i < faces.size(); i ++)
		if(faces[i].first)
			UpdateCentroid(i);
}

void SlabMesh::UpdateNormal(unsigned fid)
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

void SlabMesh::ComputeFacesNormal()
{
	for(unsigned i = 0; i < faces.size(); i ++)
		if(faces[i].first)
			UpdateNormal(i);
}

void SlabMesh::UpdateVertexNormal(unsigned vid)
{
	if(!vertices[vid].first)
		return;

	Vector3d v[3];
	std::set<unsigned> fs = vertices[vid].second->faces_;

	Vector3d vnormal;

	for (set<unsigned>::iterator si = fs.begin(); si != fs.end(); si++)
		vnormal += faces[*si].second->normal;

	vertices[vid].second->normal = vnormal;
	vertices[vid].second->normal.Normalize();
}

void SlabMesh::ComputeVerticesNormal()
{
	for(unsigned i = 0; i < vertices.size(); i ++)
		if(vertices[i].first)
			UpdateVertexNormal(i);
}

void SlabMesh::GetNeighborVertices(unsigned vid, std::set<unsigned> & neighborvertices)
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

void SlabMesh::GetLinkedEdges(unsigned eid, std::set<unsigned> & neighboredges)
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

void SlabMesh::GetAdjacentFaces(unsigned fid, std::set<unsigned> & neighborfaces)
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

bool SlabMesh::Contractible(unsigned vid_src, unsigned vid_tgt)
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

bool SlabMesh::MergeVertices(unsigned vid_src1, unsigned vid_src2, unsigned &vid_tgt)
{
	if(vid_src1 == vid_src2)
		return false;

	unsigned eid;
	InsertVertex(new SlabVertex, vid_tgt);
	// this is during edge collaspe.

	if (vertices[vid_src1].second->saved_vertex || vertices[vid_src2].second->saved_vertex)
		vertices[vid_tgt].second->saved_vertex = true;

	//if (vertices[vid_src1].second->fake_boundary_vertex || vertices[vid_src2].second->fake_boundary_vertex)
	//	vertices[vid_tgt].second->fake_boundary_vertex = true;

	//if (vertices[vid_src1].second->boundary_vertex || vertices[vid_src2].second->boundary_vertex)
	//	vertices[vid_tgt].second->boundary_vertex = true;

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

unsigned SlabMesh::VertexIncidentEdgeCount(unsigned vid)
{
	if(!vertices[vid].first)
		return 0;
	return (unsigned)vertices[vid].second->edges_.size();
}

unsigned SlabMesh::VertexIncidentFaceCount(unsigned vid)
{
	if(!vertices[vid].first)
		return 0;
	return (unsigned)vertices[vid].second->faces_.size();
}

unsigned SlabMesh::EdgeIncidentFaceCount(unsigned eid)
{
	if(!edges[eid].first)
		return 0;
	return (unsigned)edges[eid].second->faces_.size();
}

void SlabMesh::DeleteFace(unsigned fid)
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

void SlabMesh::DeleteEdge(unsigned eid)
{
	if(!edges[eid].first)
		return;

	// �����boundary_edge����������Եĸ���
	if (edges[eid].second->fake_boundary_edge)
	{
		if(vertices[edges[eid].second->vertices_.first].first)
		{
			vertices[edges[eid].second->vertices_.first].second->boundary_edge_vec.erase(eid);
			vertices[edges[eid].second->vertices_.first].second->fake_boundary_vertex = 
				vertices[edges[eid].second->vertices_.first].second->boundary_edge_vec.size() > 0 ? true : false;
		}
		if(vertices[edges[eid].second->vertices_.second].first)
		{	
			vertices[edges[eid].second->vertices_.second].second->boundary_edge_vec.erase(eid);
			vertices[edges[eid].second->vertices_.second].second->fake_boundary_vertex = 
				vertices[edges[eid].second->vertices_.second].second->boundary_edge_vec.size() > 0 ? true : false;
		}
	}

	if(vertices[edges[eid].second->vertices_.first].first)
		vertices[edges[eid].second->vertices_.first].second->edges_.erase(eid);
	if(vertices[edges[eid].second->vertices_.second].first)
		vertices[edges[eid].second->vertices_.second].second->edges_.erase(eid);
	std::set<unsigned> faces_del;
	for(std::set<unsigned>::iterator si = edges[eid].second->faces_.begin();
		si != edges[eid].second->faces_.end(); si ++)
		faces_del.insert(*si);
	for(std::set<unsigned>::iterator si = faces_del.begin(); si != faces_del.end(); si ++)
		DeleteFace(*si);

	delete edges[eid].second;
	edges[eid].first = false;
	numEdges --;
}

void SlabMesh::DeleteVertex(unsigned vid)
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

void SlabMesh::InsertVertex(SlabVertex *vertex, unsigned &vid){
	Bool_SlabVertexPointer bvp;
	bvp.first = true;
	bvp.second = vertex;
	vid = (unsigned)vertices.size();
	bvp.second->index = vid;
	vertices.push_back(bvp);
	numVertices ++;

}

void SlabMesh::InsertEdge(unsigned vid0, unsigned vid1, unsigned & eid)
{
	if(Edge(vid0,vid1,eid))
		return;
	if (!ValidVertex(vid0) || !ValidVertex(vid1))
	{
		return;
	}
	Bool_SlabEdgePointer bep;
	bep.first = true;
	bep.second = new SlabEdge;
	bep.second->vertices_.first = vid0;
	bep.second->vertices_.second = vid1;
	vertices[vid0].second->edges_.insert((unsigned)edges.size());
	vertices[vid1].second->edges_.insert((unsigned)edges.size());
	eid = (unsigned)edges.size();
	bep.second->index = eid;
	edges.push_back(bep);
	ComputeEdgeCone(eid);
	numEdges ++;
}

void SlabMesh::InsertFace(std::set<unsigned> vset)
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

	Bool_SlabFacePointer bfp;
	bfp.first = true;
	bfp.second = new SlabFace;
	bfp.second->vertices_.insert(vid[0]);
	bfp.second->vertices_.insert(vid[1]);
	bfp.second->vertices_.insert(vid[2]);
	bfp.second->edges_.insert(eid[0]);
	bfp.second->edges_.insert(eid[1]);
	bfp.second->edges_.insert(eid[2]);
	bfp.second->index = faces.size();
	vertices[vid[0]].second->faces_.insert(faces.size());
	vertices[vid[1]].second->faces_.insert(faces.size());
	vertices[vid[2]].second->faces_.insert(faces.size());
	edges[eid[0]].second->faces_.insert(faces.size());
	edges[eid[1]].second->faces_.insert(faces.size());
	edges[eid[2]].second->faces_.insert(faces.size());
	faces.push_back(bfp);
	UpdateCentroid((unsigned)faces.size()-1);
	UpdateNormal((unsigned)faces.size()-1);
	ComputeFaceSimpleTriangles((unsigned)faces.size() - 1);
	numFaces ++;
}

void SlabMesh::CleanIsolatedVertices()
{
	for(unsigned i = 0; i < vertices.size(); i ++)
		if(vertices[i].first)
		{
			if( (vertices[i].second->edges_.size() == 0) && (vertices[i].second->faces_.size() == 0) )
			{
				delete vertices[i].second;
				vertices[i].first = false;
				numVertices --;
			}
		}
}

void SlabMesh::ComputeVertexProperty(unsigned vid)
{
	if(!vertices[vid].first)
		return;
	std::set<unsigned> neighbor_vertices;
	GetNeighborVertices(vid,neighbor_vertices);
	for(std::set<unsigned>::iterator si = neighbor_vertices.begin();
		si != neighbor_vertices.end(); si ++)
		vertices[*si].second->tag = 0;

	for(std::set<unsigned>::iterator si = vertices[vid].second->faces_.begin();
		si != vertices[vid].second->faces_.end(); si ++)
		for(std::set<unsigned>::iterator ssi = faces[*si].second->vertices_.begin();
			ssi != faces[*si].second->vertices_.end(); ssi ++)
			vertices[*ssi].second->tag ++;

	bool has_one_tag(false);
	bool has_two_plus_tag(false);

	for(std::set<unsigned>::iterator si = neighbor_vertices.begin(); si != neighbor_vertices.end(); si ++)
		if( vertices[*si].second->tag > 2)
			has_two_plus_tag = true; // non-manifold
		else if( vertices[*si].second->tag < 2)
			has_one_tag = true;

	vertices[vid].second->is_boundary = has_one_tag?true:false;
	vertices[vid].second->is_disk = (has_one_tag||has_two_plus_tag)?false:true;
	vertices[vid].second->is_non_manifold = has_two_plus_tag?true:false;
}

void SlabMesh::ComputeVerticesProperty()
{
	for(unsigned i = 0; i < vertices.size(); i ++)
		if(vertices[i].first)
			ComputeVertexProperty(i);
}

void SlabMesh::ComputeEdgeCone(unsigned eid)
{
	if(!edges[eid].first)
		return;

	// test validation
	Vector3d c0 = vertices[edges[eid].second->vertices_.first].second->sphere.center;
	Vector3d c1 = vertices[edges[eid].second->vertices_.second].second->sphere.center;
	double r0 = vertices[edges[eid].second->vertices_.first].second->sphere.radius;
	double r1 = vertices[edges[eid].second->vertices_.second].second->sphere.radius;
	Vector3d c0c1 = c1-c0;
	double templeng = c0c1.Length() - abs(r1 - r0);


	Cone newc(vertices[edges[eid].second->vertices_.first].second->sphere.center, vertices[edges[eid].second->vertices_.first].second->sphere.radius,
		vertices[edges[eid].second->vertices_.second].second->sphere.center, vertices[edges[eid].second->vertices_.second].second->sphere.radius);
	edges[eid].second->cone = newc;
	if(newc.type == 1)
		edges[eid].second->valid_cone = false;
	else
		edges[eid].second->valid_cone = true;
}

void SlabMesh::ComputeEdgesCone()
{
	for(unsigned i = 0; i < edges.size(); i ++)
		if(edges[i].first)
			ComputeEdgeCone(i);
}

void SlabMesh::ComputeFaceSimpleTriangles(unsigned fid)
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

void SlabMesh::ComputeFacesSimpleTriangles()
{
	for(unsigned i = 0; i < faces.size(); i ++)
		if(faces[i].first)
			ComputeFaceSimpleTriangles(i);
}

void SlabMesh::DistinguishVertexType()
{
	for (unsigned i = 0; i != edges.size(); i++)
	{
		if (edges[i].first && (edges[i].second->faces_.size() == 1 || edges[i].second->faces_.size() == 0))
			//if (edges[i].first && edges[i].second->faces_.size() == 1)
		{
			// fake boundary edge and fake boundary vertex
			edges[i].second->fake_boundary_edge = true;
			vertices[edges[i].second->vertices_.first].second->fake_boundary_vertex = true;
			vertices[edges[i].second->vertices_.second].second->fake_boundary_vertex = true;
			vertices[edges[i].second->vertices_.first].second->boundary_edge_vec.insert(i);
			vertices[edges[i].second->vertices_.second].second->boundary_edge_vec.insert(i);
			//boundary_vertexes.insert(edges[i].second->vertices_.first);
			//boundary_vertexes.insert(edges[i].second->vertices_.second);
		}
		else if (edges[i].first && edges[i].second->faces_.size() >= 3)
		{
			// non manifold edge
			edges[i].second->non_manifold_edge = true;
			vertices[edges[i].second->vertices_.first].second->non_manifold_vertex = true;
			vertices[edges[i].second->vertices_.second].second->non_manifold_vertex = true;
			//vertices[edges[i].second->vertices_.first].second->collaspe_weight += edges[i].second->faces_.size() - 2;
			//vertices[edges[i].second->vertices_.second].second->collaspe_weight += edges[i].second->faces_.size() - 2;
		}
	}

#if 0
	// real boundary edge and vertex, saved vertex
	for (unsigned i = 0; i != vertices.size(); i++)
	{
		if (vertices[i].first && vertices[i].second->fake_boundary_vertex)
		{
			//vertices[i].second->collaspe_weight += 1.0;

			set<unsigned>::iterator si = vertices[i].second->edges_.begin();
			vector<unsigned> fake_boundary_edge_number;
			for (; si != vertices[i].second->edges_.end(); si++)
				if (edges[*si].second->fake_boundary_edge)
					fake_boundary_edge_number.push_back(*si);

			if (fake_boundary_edge_number.size() >= 2)
				vertices[i].second->boundary_vertex = true;

			set<unsigned> face_set;
			unsigned face_number = 0;
			for (unsigned j = 0; j < fake_boundary_edge_number.size(); j++)
			{
				edges[fake_boundary_edge_number[j]].second->boundary_edge = true;

				face_number += edges[fake_boundary_edge_number[j]].second->faces_.size();
				face_set.insert(edges[fake_boundary_edge_number[j]].second->faces_.begin(), 
					edges[fake_boundary_edge_number[j]].second->faces_.end());
			}

			if (face_number != face_set.size())
			{
				vertices[i].second->collaspe_weight += 10.0 * edges[fake_boundary_edge_number[0]].second->sphere.center.Dot(edges[fake_boundary_edge_number[1]].second->sphere.center) 
					/ edges[fake_boundary_edge_number[0]].second->sphere.center.Length() / edges[fake_boundary_edge_number[1]].second->sphere.center.Length();

				vertices[i].second->saved_vertex = true;
			}
		}
	}

	// judge the bounding points of each boundary point
	for(std::set<unsigned>::iterator sit = boundary_vertexes.begin(); sit != boundary_vertexes.end(); sit++)
		for(std::set<unsigned>::iterator sit2 = boundary_vertexes.begin(); sit2 != boundary_vertexes.end(); sit2++)
		{
			if (*sit == *sit2)
				continue;
			Vector3d tempvec = vertices[*sit].second->sphere.center - vertices[*sit2].second->sphere.center;

			//if (tempvec.Length() <= vertices[*sit].second->sphere.radius) 
			if (tempvec.Length() <= vertices[*sit].second->sphere.radius) 
			{
				vertices[*sit].second->bound_point_vec.push_back(*sit2);
				vertices[*sit].second->boundVec += tempvec;
			}
		}

		for(std::set<unsigned>::iterator sit = boundary_vertexes.begin(); sit != boundary_vertexes.end(); sit++)
		{
			bound_vector.push_back(pair<unsigned, unsigned>(*sit, vertices[*sit].second->bound_point_vec.size()));
			if (vertices[*sit].second->bound_point_vec.size() == 0)
				boundVec_vector.push_back(pair<unsigned, double>(*sit, 0));
			else
				boundVec_vector.push_back(pair<unsigned, double>(*sit, vertices[*sit].second->boundVec.Length() / vertices[*sit].second->bound_point_vec.size()));
			//boundVec_vector.push_back(pair<unsigned, unsigned>(*sit, vertices[*sit].second->boundVec.Length()));
		}

		sort(bound_vector.begin(), bound_vector.end(), cmpByValue<unsigned>());
		sort(boundVec_vector.begin(), boundVec_vector.end(), cmpByBiggerValue<double>());

#endif
}

unsigned SlabMesh::GetSavedPointNumber()
{
	//unsigned count = vertices.size();
	//for (unsigned i = 0; i < count; i++)
	//{
	//	if (vertices[i].first)
	//	{
	//		if (vertices[i].second->edges_.size() == 1 && vertices[i].second->faces_.size() == 0)
	//		{
	//			if (vertices[i].second->sphere.radius > 1.0e-3)
	//			{
	//				InsertSavedPoint(vertices[i].second->index);
	//			}

	//		}
	//	}
	//}
	//return count;

	unsigned count = edges.size();
	for (unsigned i = 0; i < count; i++)
	{
		if (edges[i].first)
		{
			if (edges[i].second->faces_.size() == 0)
			{
				if (vertices[edges[i].second->vertices_.first].second->edges_.size() == 1)
				{
					InsertSavedPoint(edges[i].second->vertices_.first);
					//if (vertices[edges[i].second->vertices_.first].second->sphere.radius > 1.0e-3)
					//	InsertSavedPoint(edges[i].second->vertices_.first);
					//else if (vertices[edges[i].second->vertices_.second].second->sphere.radius > 1.0e-3)
					//	InsertSavedPoint(edges[i].second->vertices_.second);
				}else if (vertices[edges[i].second->vertices_.second].second->edges_.size() == 1)
				{
					InsertSavedPoint(edges[i].second->vertices_.second);
					//if (vertices[edges[i].second->vertices_.second].second->sphere.radius > 1.0e-3)
					//	InsertSavedPoint(edges[i].second->vertices_.second);
					//else if (vertices[edges[i].second->vertices_.first].second->sphere.radius > 1.0e-3)
					//	InsertSavedPoint(edges[i].second->vertices_.first);
				}
			}
		}
	}
	return count;
}

unsigned SlabMesh::GetConnectPointNumber()
{
	unsigned count = 0;
	for (unsigned i = 0; i < edges.size(); i++)
	{
		if (edges[i].first && edges[i].second->faces_.size() == 0)
		{
			if (vertices[edges[i].second->vertices_.first].second->edges_.size() <= 3)
				count++;

			if (vertices[edges[i].second->vertices_.second].second->edges_.size() <= 3)
				count++;
		}
	}
	return count;
}

void SlabMesh::InsertSavedPoint(unsigned vid)
{
	if (vertices[vid].second->saved_vertex)
		return;

	unsigned vid_tgt;
	InsertVertex(new SlabVertex, vid_tgt);

	std::vector< std::set<unsigned> > tri_vec;
	for(std::set<unsigned>::iterator si = vertices[vid].second->faces_.begin();
		si != vertices[vid].second->faces_.end(); si ++)
		if(!faces[*si].second->HasVertex(vid_tgt))
		{
			std::set<unsigned> vset = faces[*si].second->vertices_;
			vset.erase(vid);
			vset.insert(vid_tgt);
			tri_vec.push_back(vset);
		}

		std::vector< std::pair<unsigned,unsigned> > edge_vec;
		for(std::set<unsigned>::iterator si = vertices[vid].second->edges_.begin();
			si != vertices[vid].second->edges_.end(); si ++)
			if(!edges[*si].second->HasVertex(vid_tgt))
			{
				std::pair<unsigned, unsigned> vp = edges[*si].second->vertices_;
				if(vp.first == vid)
					vp.first = vid_tgt;
				if(vp.second == vid)
					vp.second = vid_tgt;
				edge_vec.push_back(vp);
			}

			vertices[vid_tgt].second->saved_vertex = true;
			vertices[vid_tgt].second->sphere = vertices[vid].second->sphere;
			vertices[vid_tgt].second->bplist = vertices[vid].second->bplist;

			vertices[vid_tgt].second->slab_A = vertices[vid].second->slab_A;
			vertices[vid_tgt].second->slab_b = vertices[vid].second->slab_b;
			vertices[vid_tgt].second->slab_c = vertices[vid].second->slab_c;

			vertices[vid_tgt].second->mean_square_error = vertices[vid].second->mean_square_error;
			vertices[vid_tgt].second->related_face = vertices[vid].second->related_face;

			DeleteVertex(vid);

			for(unsigned i = 0; i < tri_vec.size(); i ++)
				InsertFace(tri_vec[i]);

			for(unsigned i = 0; i < edge_vec.size(); i ++)
			{
				unsigned neweid;
				InsertEdge(edge_vec[i].first, edge_vec[i].second, neweid);
			}

			for (std::set<unsigned>::iterator si = vertices[vid_tgt].second->edges_.begin(); si != vertices[vid_tgt].second->edges_.end(); si ++)
			{
				EvaluateEdgeCollapseCost(*si);
				if (edges[*si].second->collapse_cost != DBL_MAX)
				{
					edge_collapses_queue.push(EdgeInfo(*si, edges[*si].second->collapse_cost));
					ComputeEdgeCone(*si);
				}
			}

			return;
}

// �ж��Ƿ����������η�ת���
bool SlabMesh::Contractible(unsigned vid_src1, unsigned vid_src2, Vector3d &v_tgt)
{
	if( !vertices[vid_src1].first || !vertices[vid_src2].first )
		return false;

	set<unsigned> fs1;
	set<unsigned> fs2;
	fs1 = vertices[vid_src1].second->faces_;
	fs2 = vertices[vid_src2].second->faces_;

	for (std::set<unsigned>::iterator si = fs1.begin(); si != fs1.end(); si++)
	{
		if(faces[*si].first)
		{
			if ( faces[*si].second->vertices_.find(vid_src1) != faces[*si].second->vertices_.end() 
				&& faces[*si].second->vertices_.find(vid_src2) != faces[*si].second->vertices_.end())
				continue;

			Vector3d pp[3], pa[3];
			unsigned count = 0;
			for(std::set<unsigned>::iterator si2 = faces[*si].second->vertices_.begin();
				si2 != faces[*si].second->vertices_.end(); si2 ++)
			{
				pp[count] = vertices[*si2].second->sphere.center;
				pa[count++] = (*si2 != vid_src1) ? vertices[*si2].second->sphere.center : v_tgt;
			}
			Vector3d pnorm = TriangleNormal(pp[0],pp[1],pp[2]);
			Vector3d anorm = TriangleNormal(pa[0],pa[1],pa[2]);

			//double angle = VectorAngle(pnorm, anorm);
			//if (angle > Wm4::Math<double>::PI * 2.0 / 3.0)
			//	return false;

			if(pnorm.Dot(anorm) < 0)
				return false;
		}
	}

	for (std::set<unsigned>::iterator si = fs2.begin(); si != fs2.end(); si++)
	{
		if(faces[*si].first)
		{
			if ( faces[*si].second->vertices_.find(vid_src1) != faces[*si].second->vertices_.end() 
				&& faces[*si].second->vertices_.find(vid_src2) != faces[*si].second->vertices_.end())
				continue;

			Vector3d pp[3], pa[3];
			unsigned count = 0;
			for(std::set<unsigned>::iterator si2 = faces[*si].second->vertices_.begin();
				si2 != faces[*si].second->vertices_.end(); si2 ++)
			{
				pp[count] = vertices[*si2].second->sphere.center;
				pa[count++] = (*si2 != vid_src2) ? vertices[*si2].second->sphere.center : v_tgt;
			}
			Vector3d pnorm = TriangleNormal(pp[0],pp[1],pp[2]);
			Vector3d anorm = TriangleNormal(pa[0],pa[1],pa[2]);

			//double angle = VectorAngle(pnorm, anorm);
			//if (angle > Wm4::Math<double>::PI * 2.0 / 3.0)
			//	return false;

			if(pnorm.Dot(anorm) < 0)
				return false;
		}
	}

	return true;
}

bool SlabMesh::MinCostBoundaryEdgeCollapse(unsigned & eid)
{
	//merge 2 vertices of the edge first, then move the combined vertex to the preferred point and resize it.
	unsigned v1, v2;
	v1 = edges[eid].second->vertices_.first;
	v2 = edges[eid].second->vertices_.second;
	Wm4::Matrix4d A = edges[eid].second->slab_A;
	Wm4::Vector4d b = edges[eid].second->slab_b;
	double c = edges[eid].second->slab_c;
	Sphere sphere = edges[eid].second->sphere;

	if (prevent_inversion == true)
	{
		if (!Contractible(v1, v2, sphere.center))
			return false;
	}


	set<unsigned> temp_bplist;
	for (set<unsigned>::iterator it = vertices[v1].second->bplist.begin(); it != vertices[v1].second->bplist.end(); it++)
		temp_bplist.insert(*it);
	for (set<unsigned>::iterator it = vertices[v2].second->bplist.begin(); it != vertices[v2].second->bplist.end(); it++)
		temp_bplist.insert(*it);	

	unsigned temp_related_face = vertices[v1].second->related_face + vertices[v2].second->related_face;
	double temp_mean_squre_error = edges[eid].second->collapse_cost < 0 ? 0 : edges[eid].second->collapse_cost / temp_related_face;
	max_mean_squre_error = max(temp_mean_squre_error, max_mean_squre_error);

	unsigned former_edge_number = edges.size();
	unsigned vid_tgt;
	if(MergeVertices(v1, v2, vid_tgt)){
		vertices[vid_tgt].second->slab_A = A;
		vertices[vid_tgt].second->slab_b = b;
		vertices[vid_tgt].second->slab_c = c;
		vertices[vid_tgt].second->sphere = sphere;
		vertices[vid_tgt].second->related_face = temp_related_face;
		vertices[vid_tgt].second->mean_square_error = temp_mean_squre_error;
		vertices[vid_tgt].second->bplist = temp_bplist;

		switch(boundary_compute_scale)
		{
		case 1:			
			// ������ӵı��ж�����
			for (unsigned i = former_edge_number; i < edges.size(); i++)
			{
				if (edges[i].second->faces_.size() <= 1)
				{
					edges[i].second->fake_boundary_edge = true;
					vertices[edges[i].second->vertices_.first].second->fake_boundary_vertex = true;
					vertices[edges[i].second->vertices_.second].second->fake_boundary_vertex = true;
					vertices[edges[i].second->vertices_.first].second->boundary_edge_vec.insert(i);
					vertices[edges[i].second->vertices_.second].second->boundary_edge_vec.insert(i);
				}
			}
			//for (std::set<unsigned>::iterator si = vertices[vid_tgt].second->edges_.begin(); si != vertices[vid_tgt].second->edges_.end(); si ++)
			//{
			//	if (edges[*si].second->faces_.size() <= 1)
			//	{
			//		unsigned fir = edges[*si].second->vertices_.first;
			//		unsigned sec = edges[*si].second->vertices_.second;
			//		vertices[fir].second->fake_boundary_vertex = true;
			//		vertices[sec].second->fake_boundary_vertex = true;
			//	}
			//}
			break;
		case 2:
			// ������ӵı��ж�����
			for (unsigned i = former_edge_number; i < edges.size(); i++)
			{
				if (edges[i].second->faces_.size() <= 1)
				{
					edges[i].second->fake_boundary_edge = true;
					vertices[edges[i].second->vertices_.first].second->fake_boundary_vertex = true;
					vertices[edges[i].second->vertices_.second].second->fake_boundary_vertex = true;
					vertices[edges[i].second->vertices_.first].second->boundary_edge_vec.insert(i);
					vertices[edges[i].second->vertices_.second].second->boundary_edge_vec.insert(i);
				}
			}
			//for (std::set<unsigned>::iterator si = vertices[vid_tgt].second->edges_.begin(); si != vertices[vid_tgt].second->edges_.end(); si ++)
			//{
			//	if (edges[*si].second->faces_.size() <= 1)
			//	{
			//		unsigned fir = edges[*si].second->vertices_.first;
			//		unsigned sec = edges[*si].second->vertices_.second;
			//		vertices[fir].second->fake_boundary_vertex = true;
			//		vertices[sec].second->fake_boundary_vertex = true;
			//	}
			//}
			break;
		case 3:
			break;
		default:
			break;
		}

		if (compute_hausdorff == true)
		{
			double temp_sum_haus_dis = meanhausdorff_distance * pmesh->pVertexList.size();
			for (set<unsigned>::iterator it = temp_bplist.begin(); it != temp_bplist.end(); it++)
			{
				unsigned temp_ind = *it;
				Vector3d bou_ver(pmesh->pVertexList[temp_ind]->point()[0], pmesh->pVertexList[temp_ind]->point()[1], pmesh->pVertexList[temp_ind]->point()[2]);

				temp_sum_haus_dis -= pmesh->pVertexList[temp_ind]->slab_hausdorff_dist;

				double min_dis = DBL_MAX;
				unsigned min_index = -1;
				for (unsigned j = 0; j < vertices.size(); j++)
				{
					if (vertices[j].first)
					{
						Sphere ma_ver = vertices[j].second->sphere;
						double temp_length = abs((bou_ver - ma_ver.center).Length() - ma_ver.radius);
						if (temp_length >= 0 && temp_length < min_dis)
						{
							min_dis = temp_length;
							min_index = j;
						}
					}
				}

				if (min_index != -1)
				{	
					double temp_near_dis = NearestPoint(bou_ver, min_index);
					min_dis = min(temp_near_dis, min_dis);

					maxhausdorff_distance = max(maxhausdorff_distance, min_dis);
					pmesh->pVertexList[temp_ind]->slab_hansdorff_index = min_index;
					pmesh->pVertexList[temp_ind]->slab_hausdorff_dist = min_dis;

					temp_sum_haus_dis += min_dis;
				}
				else
				{
					pmesh->pVertexList[temp_ind]->slab_hansdorff_index = vid_tgt;
					double temp_len = abs((vertices[vid_tgt].second->sphere.center - bou_ver).Length() - vertices[vid_tgt].second->sphere.radius);
					if (temp_len >= 0)
					{
						pmesh->pVertexList[temp_ind]->slab_hausdorff_dist = temp_len;
						maxhausdorff_distance = max(maxhausdorff_distance,temp_len);
						temp_sum_haus_dis += temp_len;
					}
				}
			}
			meanhausdorff_distance = temp_sum_haus_dis / pmesh->pVertexList.size();
		}

		for (std::set<unsigned>::iterator si = vertices[vid_tgt].second->edges_.begin(); si != vertices[vid_tgt].second->edges_.end(); si ++)
		{
			unsigned fir = edges[*si].second->vertices_.first;
			unsigned sec = edges[*si].second->vertices_.second;

			switch(boundary_compute_scale)
			{
			case 1:
				if (!vertices[fir].second->fake_boundary_vertex || !vertices[sec].second->fake_boundary_vertex)
					continue;
				break;
			case 2:
				if (!vertices[fir].second->fake_boundary_vertex && !vertices[sec].second->fake_boundary_vertex)
					//if (vertices[fir].second->boundary_edge_vec.size() < 2 && vertices[sec].second->boundary_edge_vec.size() < 2)
					continue;
				break;
			case 3:
				if (!vertices[fir].second->fake_boundary_vertex && !vertices[sec].second->fake_boundary_vertex)
					continue;
				break;
			default:
				break;
			}
			EvaluateEdgeHausdorffCost(*si);
			//ReEvaluateEdgeHausdorffCost(*si);
			boundary_edge_collapses_queue.push(EdgeInfo(*si, edges[*si].second->collapse_cost));
		}
	}

	return true;
}

bool SlabMesh::MinCostEdgeCollapse(unsigned & eid){
	//merge 2 vertices of the edge first, then move the combined vertex to the preferred point and resize it.
	unsigned v1, v2;
	v1 = edges[eid].second->vertices_.first;
	v2 = edges[eid].second->vertices_.second;
	Wm4::Matrix4d A = edges[eid].second->slab_A;
	Wm4::Vector4d b = edges[eid].second->slab_b;
	double c = edges[eid].second->slab_c;
	Sphere sphere = edges[eid].second->sphere;
	double hyperbolic_weight = vertices[v1].second->hyperbolic_weight + vertices[v2].second->hyperbolic_weight;

	//// ���ںϲ��ᷢ�����˸ı�ıߣ���������кϲ�

	if (!edges[eid].second->topo_contractable)
		return false;
	int stopcheck;
	cout << "checking...." << endl;
	//cin >> stopcheck;
	
	// ��������˷�ת�Ĵ���ʽ
	if (prevent_inversion == true)
	{
		cout << "merging..." << endl;
		cin >> stopcheck;

		// ��������תʱ��ѡȡû������ת�ķ�ʽ���кϲ�
		// 
		if (!Contractible(v1, v2, sphere.center))
		{
			Wm4::Vector4d lamdar;
			double coll_cost = 0.0;

			int count = 0;		
			double *collapse_costs = new double[3];
			Sphere *min_sphere = new Sphere[3];
			Vector4d min_vertex;
			int min_index = 0;
			if (Contractible(v1, v2, vertices[v1].second->sphere.center))
			{
				cout << "checking ... " << endl;
				if (vertices[v1].second->is_selected == 1) {
					cout << "I catch you at v1" << endl;
					int stops;
					cin >> stops;
				}
				min_sphere[count] = vertices[v1].second->sphere;
				min_vertex = Vector4d(min_sphere[count].center.X(), min_sphere[count].center.Y(), min_sphere[count].center.Z(), min_sphere[count].radius);
				collapse_costs[count] = 0.5 * (min_vertex * A).Dot(min_vertex) - b.Dot(min_vertex) + c;
				count++;
			}
			if (Contractible(v1, v2, vertices[v2].second->sphere.center))
			{
				cout << "checking ... " << endl;
				if (vertices[v2].second->is_selected == 1) {
					cout << "I catch you at v2" << endl;
					int stops2;
					cin >> stops2;
				}
				min_sphere[count] = vertices[v2].second->sphere;
				min_vertex = Vector4d(min_sphere[count].center.X(), min_sphere[count].center.Y(), min_sphere[count].center.Z(), min_sphere[count].radius);
				collapse_costs[count] = 0.5 * (min_vertex * A).Dot(min_vertex) - b.Dot(min_vertex) + c;
				count++;
			}
			if (Contractible(v1, v2, (vertices[v1].second->sphere.center + vertices[v2].second->sphere.center) / 2.0))
			{
				cout << "checking ... " << endl;
				if (vertices[v2].second->is_selected == 1) {
					cout << "I catch you at mid" << endl;
					int stops2;
					cin >> stops2;
				}
				min_sphere[count] = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
				min_vertex = Vector4d(min_sphere[count].center.X(), min_sphere[count].center.Y(), min_sphere[count].center.Z(), min_sphere[count].radius);
				collapse_costs[count] = 0.5 * (min_vertex * A).Dot(min_vertex) - b.Dot(min_vertex) + c;
				count++;
			}
			cout << "finally, you entered here" << endl;
			if (count == 1)
			{
				lamdar = Vector4d(min_sphere[0].center.X(), min_sphere[0].center.Y(), min_sphere[0].center.Z(), min_sphere[0].radius);
				coll_cost = collapse_costs[0];
			}else if (count == 2)
			{
				min_index = collapse_costs[0] > collapse_costs[1] ? 1 : 0;
				lamdar = Vector4d(min_sphere[min_index].center.X(), min_sphere[min_index].center.Y(), min_sphere[min_index].center.Z(), min_sphere[min_index].radius);
				coll_cost = collapse_costs[min_index];
			}else if (count == 3)
			{
				if (collapse_costs[0] >= collapse_costs[1]) min_index = 1;
				min_index = collapse_costs[min_index] > collapse_costs[2] ? 2 : min_index;
				lamdar = Vector4d(min_sphere[min_index].center.X(), min_sphere[min_index].center.Y(), min_sphere[min_index].center.Z(), min_sphere[min_index].radius);
				coll_cost = collapse_costs[min_index];
			}else
				coll_cost += 1e9;
			delete [] collapse_costs;
			delete [] min_sphere;

			edges[eid].second->qem_error = coll_cost;
			//coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight 
			//	* edges[eid].second->hyperbolic_weight	* edges[eid].second->hyperbolic_weight
			//	* edges[eid].second->hyperbolic_weight	* edges[eid].second->hyperbolic_weight;

			coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight;

			edges[eid].second->collapse_cost = coll_cost;

			edges[eid].second->sphere.center = Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z());
			edges[eid].second->sphere.radius = lamdar.W();

			edge_collapses_queue.push(EdgeInfo(eid, coll_cost));

			return false;
		}
	}

	cout << "jumped" << endl;
	// �����Ǳ߽�߽��м򻯻����ڲ��߽��м� --> ֻ�������ж�
	if (edges[eid].second->faces_.size() <= 1)
		simplified_boundary_edges++;
	else
		simplified_inside_edges++;
	// ÿ�μ���1000����֮������򻯽��
	if (simplified_boundary_edges + simplified_inside_edges == 1000)
	{
		ExportSimplifyResult();
		simplified_inside_edges = simplified_boundary_edges = 0;
	}

	set<unsigned> temp_bplist;
	if (compute_hausdorff)
	{
		for (set<unsigned>::iterator it = vertices[v1].second->bplist.begin(); it != vertices[v1].second->bplist.end(); it++)
			temp_bplist.insert(*it);
		for (set<unsigned>::iterator it = vertices[v2].second->bplist.begin(); it != vertices[v2].second->bplist.end(); it++)
			temp_bplist.insert(*it);
	}

	unsigned temp_related_face = vertices[v1].second->related_face + vertices[v2].second->related_face;
	double temp_mean_squre_error = edges[eid].second->collapse_cost < 0 ? 0 : edges[eid].second->collapse_cost / temp_related_face;
	max_mean_squre_error = max(temp_mean_squre_error, max_mean_squre_error);

	unsigned vid_tgt;
	if(MergeVertices(v1, v2, vid_tgt)){  
		// checking the merging...
		cout << "inside the merge_verices. " << endl;
		cout << v1 << " -- " << v2 << " -- " << vid_tgt << endl;

		// int jj;
		// cin >> jj;
		vertices[vid_tgt].second->slab_A = A;
		vertices[vid_tgt].second->slab_b = b;
		vertices[vid_tgt].second->slab_c = c;
		vertices[vid_tgt].second->sphere = sphere;
		vertices[vid_tgt].second->related_face = temp_related_face;
		vertices[vid_tgt].second->mean_square_error = temp_mean_squre_error;
		vertices[vid_tgt].second->hyperbolic_weight = hyperbolic_weight;

		// ����������Ϣ
		InitialTopologyProperty(vid_tgt);

		for (std::set<unsigned>::iterator si = vertices[vid_tgt].second->edges_.begin(); si != vertices[vid_tgt].second->edges_.end(); si ++)
		{
			unsigned fir = edges[*si].second->vertices_.first;
			unsigned sec = edges[*si].second->vertices_.second;

			EvaluateEdgeCollapseCost(*si);
			ComputeEdgeCone(*si);
			edge_collapses_queue.push(EdgeInfo(*si, edges[*si].second->collapse_cost));
		}

		if (compute_hausdorff)
		{
			double temp_sum_haus_dis = meanhausdorff_distance * pmesh->pVertexList.size();
			for (set<unsigned>::iterator it = temp_bplist.begin(); it != temp_bplist.end(); it++)
			{
				unsigned temp_ind = *it;
				Vector3d bou_ver(pmesh->pVertexList[temp_ind]->point()[0], pmesh->pVertexList[temp_ind]->point()[1], pmesh->pVertexList[temp_ind]->point()[2]);
				bou_ver /= pmesh->bb_diagonal_length;

				temp_sum_haus_dis -= pmesh->pVertexList[temp_ind]->slab_hausdorff_dist;

				double min_dis = DBL_MAX;
				unsigned min_index = -1;
				for (unsigned j = 0; j < vertices.size(); j++)
				{
					if (vertices[j].first)
					{
						Sphere ma_ver = vertices[j].second->sphere;
						double temp_length = abs((bou_ver - ma_ver.center).Length() - ma_ver.radius);
						if (temp_length >= 0 && temp_length < min_dis)
						{
							min_dis = temp_length;
							min_index = j;
						}
					}
				}

				if (min_index != -1)
				{	
					double temp_near_dis = NearestPoint(bou_ver, min_index);
					min_dis = min(temp_near_dis, min_dis);

					vertices[min_index].second->bplist.insert(temp_ind);
					maxhausdorff_distance = max(maxhausdorff_distance, min_dis);
					pmesh->pVertexList[temp_ind]->slab_hansdorff_index = min_index;
					pmesh->pVertexList[temp_ind]->slab_hausdorff_dist = min_dis;

					temp_sum_haus_dis += min_dis;
				}
				else
				{
					pmesh->pVertexList[temp_ind]->slab_hansdorff_index = vid_tgt;
					double temp_len = abs((vertices[vid_tgt].second->sphere.center - bou_ver).Length() - vertices[vid_tgt].second->sphere.radius);
					if (temp_len >= 0)
					{
						pmesh->pVertexList[temp_ind]->slab_hausdorff_dist = temp_len;
						maxhausdorff_distance = max(maxhausdorff_distance,temp_len);
						temp_sum_haus_dis += temp_len;
					}
				}
			}
			meanhausdorff_distance = temp_sum_haus_dis / pmesh->pVertexList.size();
		}
	}

	return true;
}

void SlabMesh::EvaluateEdgeCollapseCost(unsigned eid){
	//Modified by Zhiyang ...
// this is key part of cost assigning...
	if (!edges[eid].first)
		return ;

	unsigned v1, v2;
	v1 = edges[eid].second->vertices_.first;
	v2 = edges[eid].second->vertices_.second; 

	double weight = vertices[v1].second->hyperbolic_weight + vertices[v2].second->hyperbolic_weight;

	// set the hyperbolic weight to the related edge
	switch(hyperbolic_weight_type)
	{
	case 1:
		edges[eid].second->hyperbolic_weight = GetHyperbolicLength(eid);
		break;
	case 2:
		edges[eid].second->hyperbolic_weight = GetHyperbolicLength(eid);
		break;
	case 3:
		edges[eid].second->hyperbolic_weight = GetRatioHyperbolicEuclid(eid);
		break;
	default:
		break;
	}


	//double w1 = 1e-5, w2 = 1e-5;
	double w1 = 1.0, w2 = 1.0;
	//// �Բ�ͬratio�ı߽���ӳ�䴦��С��0.2�Ĳ�������(0.2,1)ӳ�䵽(2, 10)
	//w1 = edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight;
	//w2 = w1;
	//if (edges[eid].second->hyperbolic_weight >= 0.2)
	//{
	//	w1 = (edges[eid].second->hyperbolic_weight - 0.2) / 0.8 * (10 - 2) + 2;
	//	w2 = w1;
	//}

	Wm4::Matrix4d A1 = vertices[v1].second->slab_A;
	Wm4::Matrix4d A2 = vertices[v2].second->slab_A; 
	Wm4::Matrix4d add_A1 = vertices[v1].second->add_A * w1;
	Wm4::Matrix4d add_A2 = vertices[v2].second->add_A * w2;

	Wm4::Vector4d b1 = vertices[v1].second->slab_b;
	Wm4::Vector4d b2 = vertices[v2].second->slab_b;	
	Wm4::Vector4d add_b1 = vertices[v1].second->add_b * w1;
	Wm4::Vector4d add_b2 = vertices[v2].second->add_b * w2;

	double c1 = vertices[v1].second->slab_c;
	double c2 = vertices[v2].second->slab_c;
	double add_c1 = vertices[v1].second->add_c * w1;
	double add_c2 = vertices[v2].second->add_c * w2;

	edges[eid].second->slab_A = A1 + A2;
	edges[eid].second->slab_b = b1 + b2;
	edges[eid].second->slab_c = c1 + c2;

	Matrix4d inverse_A_matrix = edges[eid].second->slab_A.Inverse();
	Wm4::Vector4d lamdar;
	double coll_cost = 0.0;

	if ((vertices[v1].second->saved_vertex && !vertices[v2].second->saved_vertex) ||
		(vertices[v2].second->saved_vertex && !vertices[v1].second->saved_vertex))
	{
		cout << "This fixed vertex." << endl; // This is fixed!
		// one of the two vertices is fixed sphere.
		Sphere mid_sphere; 
		if (vertices[v1].second->saved_vertex)
			mid_sphere = vertices[v1].second->sphere;
		else
			mid_sphere = vertices[v2].second->sphere;

		lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
	}	
	else if ((vertices[v1].second->saved_vertex && vertices[v2].second->saved_vertex))
	{
		cout << "stop here!" << endl;
		cout << "the two edges are both saved_vertex..." << endl;
		cout << "inside the vertices" << endl;
		cout << "v " << vertices[v1].second->sphere.center << endl;
		cout << "v " << vertices[v2].second->sphere.center << endl;

		//cin >> sss;
		
		
		//edges[eid].second->collapse_cost = DBL_MAX;
		//GetBestBoundaryPoint(eid);
		//read it here
		double collapse_costs[3];
		Sphere min_sphere[3];
		Vector4d min_vertex;
		int min_index = 0;

		min_sphere[0] = vertices[v1].second->sphere;
		min_vertex = Vector4d(min_sphere[0].center.X(), min_sphere[0].center.Y(), min_sphere[0].center.Z(), min_sphere[0].radius);
		collapse_costs[0] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
			- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
		min_sphere[1] = vertices[v2].second->sphere;
		min_vertex = Vector4d(min_sphere[1].center.X(), min_sphere[1].center.Y(), min_sphere[1].center.Z(), min_sphere[1].radius);
		collapse_costs[1] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
			- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
		min_sphere[2] = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
		min_vertex = Vector4d(min_sphere[2].center.X(), min_sphere[2].center.Y(), min_sphere[2].center.Z(), min_sphere[2].radius);
		collapse_costs[2] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
			- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;

		if (collapse_costs[0] >= collapse_costs[1]) min_index = 1;

		min_index = collapse_costs[min_index] > collapse_costs[2] ? 2 : min_index;

		//edges[eid].second->collapse_cost = collapse_costs[min_index];
		edges[eid].second->collapse_cost = 99999999;
		edges[eid].second->sphere.center = min_sphere[min_index].center;
		edges[eid].second->sphere.radius = min_sphere[min_index].radius;

		return;
	}
	else
	{
		if (inverse_A_matrix != Matrix4d() || edges[eid].second->faces_.size() == 0)
		{
			// add the boundary preserving.
			//if (edges[eid].second->hyperbolic_weight >= 0.1)
			//{
			edges[eid].second->slab_A = A1 + A2 + add_A1 + add_A2;
			edges[eid].second->slab_b = b1 + b2 + add_b1 + add_b2;
			edges[eid].second->slab_c = c1 + c2 + add_c1 + add_c2;
			//}
			inverse_A_matrix = edges[eid].second->slab_A.Inverse();
			if (inverse_A_matrix != Matrix4d())
			{
				lamdar = inverse_A_matrix * edges[eid].second->slab_b;
				if (lamdar.W() < 0)
				{
					Sphere mid_sphere = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
					lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
				}
			}
			else
			{
				// it's now calculate as the middle of the spheres
				Sphere mid_sphere = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
				lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
			}
		}
		else
		{
			//// it's now calculate as the middle of the spheres
			//Sphere mid_sphere = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
			//lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);

			double collapse_costs[3];
			Sphere min_sphere[3];
			Vector4d min_vertex;
			int min_index = 0;

			min_sphere[0] = vertices[v1].second->sphere;
			min_vertex = Vector4d(min_sphere[0].center.X(), min_sphere[0].center.Y(), min_sphere[0].center.Z(), min_sphere[0].radius);
			collapse_costs[0] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
				- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
			min_sphere[1] = vertices[v2].second->sphere;
			min_vertex = Vector4d(min_sphere[1].center.X(), min_sphere[1].center.Y(), min_sphere[1].center.Z(), min_sphere[1].radius);
			collapse_costs[1] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
				- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
			min_sphere[2] = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
			min_vertex = Vector4d(min_sphere[2].center.X(), min_sphere[2].center.Y(), min_sphere[2].center.Z(), min_sphere[2].radius);
			collapse_costs[2] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
				- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;

			if (collapse_costs[0] >= collapse_costs[1]) min_index = 1;
			min_index = collapse_costs[min_index] > collapse_costs[2] ? 2 : min_index;

			coll_cost = collapse_costs[min_index];
			switch(hyperbolic_weight_type)
			{
			case 1:
				coll_cost = coll_cost * edges[eid].second->hyperbolic_weight;
				break;
			case 2:
				if (weight <= 1e-12)
					coll_cost = 0.0;
				else
					coll_cost = collapse_costs[min_index] / weight;
				break;
			case 3:
				//coll_cost = edges[eid].second->hyperbolic_weight;
				edges[eid].second->qem_error = coll_cost;
				coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight;

				//coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight 
				//	* edges[eid].second->hyperbolic_weight	* edges[eid].second->hyperbolic_weight
				//	* edges[eid].second->hyperbolic_weight	* edges[eid].second->hyperbolic_weight;

				//coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight;
				//coll_cost = (coll_cost + edges[eid].second->hyperbolic_weight) * edges[eid].second->hyperbolic_weight;
				break;
			default:
				break;
			}

			edges[eid].second->collapse_cost = coll_cost;
			edges[eid].second->sphere.center = min_sphere[min_index].center;
			edges[eid].second->sphere.radius = min_sphere[min_index].radius;

			return;
		}
	}
	// set the first case.?
	coll_cost = 0.5 * (lamdar * edges[eid].second->slab_A).Dot(lamdar) 
		- edges[eid].second->slab_b.Dot(lamdar) + edges[eid].second->slab_c;

	// ��������תʱ��ѡȡû������ת�ķ�ʽ���кϲ�
	if (!Contractible(v1, v2, Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z())))
	{
		int count = 0;		
		double *collapse_costs = new double[3];
		Sphere *min_sphere = new Sphere[3];
		Vector4d min_vertex;
		int min_index = 0;
		if (!Contractible(v1, v2, vertices[v1].second->sphere.center))
		{
			min_sphere[count] = vertices[v1].second->sphere;
			min_vertex = Vector4d(min_sphere[count].center.X(), min_sphere[count].center.Y(), min_sphere[count].center.Z(), min_sphere[count].radius);
			collapse_costs[count] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
				- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
			count++;
		}
		if (!Contractible(v1, v2, vertices[v2].second->sphere.center))
		{
			min_sphere[count] = vertices[v2].second->sphere;
			min_vertex = Vector4d(min_sphere[count].center.X(), min_sphere[count].center.Y(), min_sphere[count].center.Z(), min_sphere[count].radius);
			collapse_costs[count] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
				- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
			count++;
		}
		if (!Contractible(v1, v2, (vertices[v1].second->sphere.center + vertices[v2].second->sphere.center) / 2.0))
		{
			min_sphere[count] = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
			min_vertex = Vector4d(min_sphere[count].center.X(), min_sphere[count].center.Y(), min_sphere[count].center.Z(), min_sphere[count].radius);
			collapse_costs[count] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
				- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
			count++;
		}

		if (count == 1)
		{
			lamdar = Vector4d(min_sphere[0].center.X(), min_sphere[0].center.Y(), min_sphere[0].center.Z(), min_sphere[0].radius);
			//coll_cost = collapse_costs[0];
			//coll_cost = 0.0;
		}else if (count == 2)
		{
			min_index = collapse_costs[0] > collapse_costs[1] ? 1 : 0;
			//min_index = min_sphere[0].radius > min_sphere[1].radius ? 0 : 1;
			lamdar = Vector4d(min_sphere[min_index].center.X(), min_sphere[min_index].center.Y(), min_sphere[min_index].center.Z(), min_sphere[min_index].radius);
			//coll_cost = collapse_costs[min_index];
			//coll_cost = 0.0;
		}else if (count == 3)
		{
			if (collapse_costs[0] >= collapse_costs[1]) min_index = 1;
			min_index = collapse_costs[min_index] > collapse_costs[2] ? 2 : min_index;
			//if (min_sphere[0].radius >= min_sphere[1].radius) min_index = 0;
			//min_index = min_sphere[min_index].radius > min_sphere[2].radius ? min_index : 2;
			lamdar = Vector4d(min_sphere[min_index].center.X(), min_sphere[min_index].center.Y(), min_sphere[min_index].center.Z(), min_sphere[min_index].radius);
			//coll_cost = collapse_costs[min_index];
			//coll_cost = 0.0;
		}else
			coll_cost += 1e9;
		delete [] collapse_costs;
		delete [] min_sphere;
	}

	switch(hyperbolic_weight_type)
	{
	case 1:
		coll_cost = coll_cost * edges[eid].second->hyperbolic_weight;
		break;
	case 2:
		if (weight <= 1e-12)
			coll_cost = 0.0;
		else
			coll_cost = coll_cost / weight;
		break;
	case 3:
		//coll_cost = coll_cost * edges[eid].second->hyperbolic_weight;
		//coll_cost = edges[eid].second->hyperbolic_weight;
		edges[eid].second->qem_error = coll_cost;
		coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight;

		//coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight 
		//	* edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight
		//	* edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight;

		//coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight;
		//coll_cost = (coll_cost + edges[eid].second->hyperbolic_weight) * edges[eid].second->hyperbolic_weight;
		break;
	default:
		break; 
	}

	edges[eid].second->collapse_cost = coll_cost;
	edges[eid].second->sphere.center = Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z());
	edges[eid].second->sphere.radius = lamdar.W();
}

void SlabMesh::EvaluateEdgeHausdorffCost(unsigned eid)
{
	if (!edges[eid].first)
		return ;

	unsigned v1, v2;
	v1 = edges[eid].second->vertices_.first;
	v2 = edges[eid].second->vertices_.second;

	double w1 = 1.0, w2 = 1.0;

	Wm4::Matrix4d A1 = vertices[v1].second->slab_A;
	Wm4::Matrix4d A2 = vertices[v2].second->slab_A;
	Wm4::Matrix4d add_A1 = vertices[v1].second->add_A * w1;
	Wm4::Matrix4d add_A2 = vertices[v2].second->add_A * w2;

	Wm4::Vector4d b1 = vertices[v1].second->slab_b;
	Wm4::Vector4d b2 = vertices[v2].second->slab_b;	
	Wm4::Vector4d add_b1 = vertices[v1].second->add_b * w1;
	Wm4::Vector4d add_b2 = vertices[v2].second->add_b * w2;

	double c1 = vertices[v1].second->slab_c;
	double c2 = vertices[v2].second->slab_c;
	double add_c1 = vertices[v1].second->add_c * w1;
	double add_c2 = vertices[v2].second->add_c * w2;

	edges[eid].second->slab_A = A1 + A2;
	edges[eid].second->slab_b = b1 + b2;
	edges[eid].second->slab_c = c1 + c2;

	Matrix4d inverse_A_matrix = edges[eid].second->slab_A.Inverse();
	Wm4::Vector4d lamdar;
	double coll_cost;

	if (inverse_A_matrix != Matrix4d())
	{
		// add the boundary preserving.
		edges[eid].second->slab_A = A1 + A2 + add_A1 + add_A2;
		edges[eid].second->slab_b = b1 + b2 + add_b1 + add_b2;
		edges[eid].second->slab_c = c1 + c2 + add_c1 + add_c2;

		inverse_A_matrix = edges[eid].second->slab_A.Inverse();
		if (inverse_A_matrix != Matrix4d())
		{
			lamdar = inverse_A_matrix * edges[eid].second->slab_b;
			if (lamdar.W() < 0)
			{
				Sphere mid_sphere = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
				lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
			}
		}
		else
		{
			// it's now calculate as the middle of the spheres
			Sphere mid_sphere = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
			lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
		}
	}
	else
	{
		// it's now calculate as the middle of the spheres
		Sphere mid_sphere = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
		lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
	}

	set<unsigned> temp_bplist;
	for (set<unsigned>::iterator it = vertices[v1].second->bplist.begin(); it != vertices[v1].second->bplist.end(); it++)
		temp_bplist.insert(*it);
	for (set<unsigned>::iterator it = vertices[v2].second->bplist.begin(); it != vertices[v2].second->bplist.end(); it++)
		temp_bplist.insert(*it);

	double max_hausdorff = 0;
	for (set<unsigned>::iterator it = temp_bplist.begin(); it != temp_bplist.end(); it++)
	{
		unsigned temp_ind = *it;
		Vector3d bou_ver(pmesh->pVertexList[temp_ind]->point()[0], pmesh->pVertexList[temp_ind]->point()[1], pmesh->pVertexList[temp_ind]->point()[2]);

		//double min_dis = DBL_MIN;
		//unsigned min_index = -1;
		//for (unsigned j = 0; j < vertices.size(); j++)
		//{
		//	if (vertices[j].first && (j != v1 && j != v2))
		//	{
		//		Sphere ma_ver = vertices[j].second->sphere;
		//		double temp_length = abs((bou_ver - ma_ver.center).Length() - ma_ver.radius);
		//		if (temp_length >= 0 && temp_length < min_dis)
		//		{
		//			min_dis = temp_length;
		//			min_index = j;
		//		}
		//	}
		//}

		//double temp_near_dis = NearestPoint(bou_ver, min_index);
		//min_dis = min(temp_near_dis, min_dis);

		// ��������С���Ƿ��������ɵĵ�
		double temp_length = abs((bou_ver - Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z())).Length() - lamdar.W());
		//min_dis = max(temp_length, min_dis);

		max_hausdorff = max(temp_length, max_hausdorff);
	}
	//set<unsigned> neighbors_v, tdneighbors_v;
	//GetNeighborVertices(v1, neighbors_v);
	//GetNeighborVertices(v2, tdneighbors_v);
	//neighbors_v.insert(tdneighbors_v.begin(), tdneighbors_v.end());
	//neighbors_v.erase(v1);
	//neighbors_v.erase(v2);

	//set< set<unsigned> > neighbors_f;
	//for(set<unsigned>::iterator si = vertices[v1].second->faces_.begin();
	//	si != vertices[v1].second->faces_.end(); si ++)
	//{
	//	if(!faces[*si].second->HasVertex(v2))
	//	{
	//		set<unsigned> vset = faces[*si].second->vertices_;
	//		vset.erase(v1);
	//		neighbors_f.insert(vset);
	//	}
	//}
	//for(set<unsigned>::iterator si = vertices[v2].second->faces_.begin();
	//	si != vertices[v2].second->faces_.end(); si ++)
	//{
	//	if(!faces[*si].second->HasVertex(v1))
	//	{
	//		set<unsigned> vset = faces[*si].second->vertices_;
	//		vset.erase(v2);
	//		neighbors_f.insert(vset);
	//	}
	//}
	//max_hausdorff = EvaluateVertexDistanceErrorEnvelope(lamdar, neighbors_v, neighbors_f, temp_bplist);

	coll_cost = max_hausdorff;

	double temp_coll_cost = 0.5 * (lamdar * edges[eid].second->slab_A).Dot(lamdar) 
		- edges[eid].second->slab_b.Dot(lamdar) + edges[eid].second->slab_c;

	edges[eid].second->collapse_cost = coll_cost;
	edges[eid].second->sphere.center = Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z());
	edges[eid].second->sphere.radius = lamdar.W();
}

void SlabMesh::ReEvaluateEdgeHausdorffCost(unsigned eid)
{
	if (!edges[eid].first)
		return ;

	unsigned v[2];
	v[0] = edges[eid].second->vertices_.first;
	v[1] = edges[eid].second->vertices_.second;

	Wm4::Matrix4d A[2];
	Wm4::Vector4d b[2];
	double c[2]  = {0, 0};

	for(unsigned i = 0; i < 2; i++)
	{
		SlabVertex sv = *vertices[v[i]].second;
		std::set<unsigned> fset = sv.faces_;
		Vector4d C1(sv.sphere.center.X(), sv.sphere.center.Y(), sv.sphere.center.Z(), sv.sphere.radius);

		for (set<unsigned>::iterator si = fset.begin(); si != fset.end(); si++)
		{
			SlabFace sf = *faces[*si].second;

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

			A[i] += temp_A1;
			A[i] += temp_A2;
			b[i] += temp_b1;
			b[i] += temp_b2;
			c[i] += temp_c1;
			c[i] += temp_c2;
		}
	}

	edges[eid].second->slab_A = A[0] + A[1];
	edges[eid].second->slab_b = b[0] + b[1];
	edges[eid].second->slab_c = c[0] + c[1];

	Matrix4d inverse_A_matrix = edges[eid].second->slab_A.Inverse();
	Wm4::Vector4d lamdar;
	double coll_cost;

	if (inverse_A_matrix != Matrix4d())
	{
		lamdar = inverse_A_matrix * edges[eid].second->slab_b;
		if (lamdar.W() < 0)
		{
			Sphere mid_sphere = (vertices[v[0]].second->sphere + vertices[v[1]].second->sphere) * 0.5;
			lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
		}
	}
	else
	{
		// it's now calculate as the middle of the spheres
		Sphere mid_sphere = (vertices[v[0]].second->sphere + vertices[v[1]].second->sphere) * 0.5;
		lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
	}

	set<unsigned> temp_bplist;
	for (set<unsigned>::iterator it = vertices[v[0]].second->bplist.begin(); it != vertices[v[0]].second->bplist.end(); it++)
		temp_bplist.insert(*it);
	for (set<unsigned>::iterator it = vertices[v[1]].second->bplist.begin(); it != vertices[v[1]].second->bplist.end(); it++)
		temp_bplist.insert(*it);

	double max_hausdorff = 0;
	for (set<unsigned>::iterator it = temp_bplist.begin(); it != temp_bplist.end(); it++)
	{
		unsigned temp_ind = *it;
		Vector3d bou_ver(pmesh->pVertexList[temp_ind]->point()[0], pmesh->pVertexList[temp_ind]->point()[1], pmesh->pVertexList[temp_ind]->point()[2]);

		// ��������С���Ƿ��������ɵĵ�
		double temp_length = abs((bou_ver - Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z())).Length() - lamdar.W());
		//min_dis = max(temp_length, min_dis);

		max_hausdorff = max(temp_length, max_hausdorff);
	}

	coll_cost = max_hausdorff;

	double temp_coll_cost = 0.5 * (lamdar * edges[eid].second->slab_A).Dot(lamdar) 
		- edges[eid].second->slab_b.Dot(lamdar) + edges[eid].second->slab_c;

	edges[eid].second->collapse_cost = coll_cost;
	edges[eid].second->sphere.center = Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z());
	edges[eid].second->sphere.radius = lamdar.W();
}

vector<string> SlabMesh::my_split(const string &str, const string &pattern)
{
	//just for quick test...
	//�����ַ������и�
	//const char* convert to char*
	char * strc = new char[strlen(str.c_str()) + 1];
	strcpy(strc, str.c_str());
	vector<string> resultVec;
	char* tmpStr = strtok(strc, pattern.c_str());
	while (tmpStr != NULL)
	{
		resultVec.push_back(string(tmpStr));
		tmpStr = strtok(NULL, pattern.c_str());
	}

	delete[] strc;

	return resultVec;
}
void SlabMesh::readMA_ball_diff_radius(string objpath, vector<vector<double> >& vset, vector<vector<int> >& eset, vector<vector<int> >& fset) {
	/* ��ȡ��MA�ļ�
	* �ڶ�ȡ֮ǰ�ǵ�Ҫ�Ƚ���λ�����
	* */
	string line;
	fstream f;
	f.open(objpath, ios::in);
	if (!f.is_open())
		std::cout << "�ļ��򿪳���" << endl;
	int v_counter = 1;
	int f_counter = 1;
	while (!f.eof()) {
		getline(f, line);//�õ�obj�ļ���һ�У���Ϊһ���ַ���
		vector <string> parameters;
		string tailMark = " ";
		string ans = "";
		line = line.append(tailMark);
		std::cout << line << endl;
		
		if (line[0] != 'v' && line[0] != 'f' && line[0] != 'e')
			continue;
		for (int i = 0; i < line.length(); i++) {
			char ch = line[i];
			if (ch != ' ') {
				ans += ch;
			}
			else {
				if (ans != "") {
					parameters.push_back(ans); //ȡ���ַ����е�Ԫ�أ��Կո��з�
					ans = "";
				}
			}
		}
		if (parameters[0] == "v") {   //����Ƕ���Ļ�
									  //�����
			vector<double> Point;
			v_counter++;
			Point.push_back(atof(parameters[1].c_str()));
			Point.push_back(atof(parameters[2].c_str()));
			Point.push_back(atof(parameters[3].c_str()));
			Point.push_back(fabs(atof(parameters[4].c_str())));//�뾶��Ϣ
			vset.push_back(Point);
		}
		else if (parameters[0] == "e") {
			vector<int> Edge;   //ʹ��int��Ϊindex
			Edge.push_back(atof(parameters[1].c_str()));
			Edge.push_back(atof(parameters[2].c_str()));
			eset.push_back(Edge);
			//�����
		}
		else if (parameters[0] == "f") {   //�������Ļ�����Ŷ��������
										   //������
			vector<int> vIndexSets;          //��ʱ��ŵ�ļ���
			for (int i = 1; i < 4; i++) {
				string x = parameters[i];
				string ans = "";
				for (int j = 0; j < x.length(); j++) {   //������/��
					char ch = x[j];
					if (ch != '/')
						ans += ch;
					else
						break;
				}
				vector <string> res = my_split(ans, "/");
				int index = atof(res[0].c_str());
				index--;//��Ϊ����������obj�ļ����Ǵ�1��ʼ�ģ������Ǵ�ŵĶ���vector�Ǵ�0��ʼ�ģ����Ҫ��1
				vIndexSets.push_back(index);
			}
			fset.push_back(vIndexSets);
		}

	}
	f.close();
	int vert_number = vset.size();
	int face_number = fset.size();
	int edge_number = eset.size();
	std::cout << "vnumber" << vert_number << endl;
	std::cout << "fnumber" << face_number << endl;
	std::cout << "enumber" << edge_number << endl;
	//������еĵ������ܵ�һ��obj�ļ���
	//����������ǵĽ��
	

}
double SlabMesh::p_distance(double x1, double y1, double z1, double x2, double y2, double z2) {
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));

}
void SlabMesh::Simplify_with_Selected_Pole(int threshold, vector<vector<double> >& selected_pole){

	// Their setting...
	/*
	Zhiyang: 
	1. Here we first introduce our selected poles.
	2. Get the corresponding vertex index.
	3. Do not remove these points...
	*/


	vector<vector<double> > vset;

	vector<vector<int> >eset;
	vector<vector<int> > fset;

	cout << "simplifying..." << endl;



	int preserve_num = 30;
    readMA_ball_diff_radius("C://Users//frank//Desktop//cover_pics//01kitten//01kitten_selected_inner_points.txt", vset, eset, fset); // 49

	std::ofstream fsout("C://Users//frank//Desktop//all_poles.obj");
	for (int vid = 0; vid < vertices.size(); vid++) {
		fsout << "v " << vertices[vid].second->sphere.center << endl;
	}
	fsout.close();

	// then we iterate in vertices set and label them...
	for (int select_id = 0;select_id < vset.size(); select_id++) {
		double min_dis = 999999;
		int min_idx = -1;
		for (int vid = 0; vid < vertices.size(); vid++) {
			vertices[vid].second->is_non_manifold;

			vertices[vid].second->sphere.center;
			double temp_dis = p_distance(vset[select_id][0],
				vset[select_id][1],
				vset[select_id][2],
				vertices[vid].second->sphere.center[0],
				vertices[vid].second->sphere.center[1],
				vertices[vid].second->sphere.center[2]);
			if (temp_dis < min_dis 
				&& vertices[vid].second->is_selected==0 
				//&& vertices[vid].second->is_pole // No, it does not need to be pole...
				) {
				min_dis = temp_dis;
				min_idx = vid;
			}
		}
		vertices[min_idx].second->is_selected = 1;
		vertices[min_idx].second->saved_vertex = 1;
		vertices[min_idx].second->is_pole = 1;
		vector<double> mPoint;
		mPoint.push_back(vertices[min_idx].second->sphere.center[0]);
		mPoint.push_back(vertices[min_idx].second->sphere.center[1]);
		mPoint.push_back(vertices[min_idx].second->sphere.center[2]);
		selected_pole.push_back(mPoint);

	}
	int check_count = 0;
	
	cout << "how many: " << check_count << " / "<< vset.size() << endl;

	// this is for checking...
	// double checked.
	/*
	std::ofstream fsout2("C://Users//frank//Desktop//nearest_poles.obj");
	for (int vid = 0; vid < vertices.size(); vid++) {
		if (vertices[vid].second->is_selected == 1) {
			fsout2 << "v " << vertices[vid].second->sphere.center << endl;
		}
	}
	fsout2.close();
	for (int vid = 0; vid < vertices.size(); vid++) {
		if (vertices[vid].second->is_selected == 1) {
			cout << "v " << vertices[vid].second->sphere.center << endl;
			cout << "is_pole: " << vertices[vid].second->is_pole << endl;
			cout << "is_non_manifold: " << vertices[vid].second->is_non_manifold << endl;
			cout << "is_disk: " << vertices[vid].second->is_disk << endl;
			cout << "is_boundary: " << vertices[vid].second->is_boundary << endl;
			cout << "is_selected" << vertices[vid].second->is_selected << endl;
			cout << "-----------------------------" << endl;
			
		}
	}*/

	//int sss;



		// find the nearest one...
	


	// ���򻯵�С��50������ʱ������������˵�ı߽��кϲ�
	//
	if (numVertices <= preserve_num)
	{
		if (initial_boundary_preserve == false)
		{
			initial_boundary_preserve = true;
			InitialTopologyProperty();
			for (int i = 0; i < vertices.size(); i++)
			{
				if (vertices[i].first)
				{
					set<unsigned> fir_edges = vertices[i].second->edges_;
					for (set<unsigned>::iterator si = fir_edges.begin(); si != fir_edges.end(); si++)
					{
						int index = edges[*si].second->vertices_.first == i ? 
							edges[*si].second->vertices_.second : edges[*si].second->vertices_.first;

						if (vertices[index].second->edges_.size() == 1 && vertices[index].second->faces_.size() == 0)
						{
							edges[*si].second->topo_contractable = false;
						}
					}
				}
			}
		}
	}

	// this is for selected poles:
	
		
	// This is the key simplification process.
	int deleteSphereNum = 0;
	if (!boundary_edge_collapses_queue.empty())
	{
		while (deleteSphereNum < threshold && numVertices > 1 && !boundary_edge_collapses_queue.empty())
		{
			EdgeInfo topEdge = boundary_edge_collapses_queue.top();
			boundary_edge_collapses_queue.pop();
			unsigned eid = topEdge.edge_num;
			if(edges[eid].first && ValidVertex(edges[eid].second->vertices_.first) && ValidVertex(edges[eid].second->vertices_.second))
			{
				std::cout << edges[eid].second->vertices_.first << "  " << edges[eid].second->vertices_.second << endl;
				std::cout << "Under the first condition --------------------------------------" << endl;
				cout << "This is the first time we meet the situation... have a double check..." << endl;
				cout << "Slabmesh.cpp Line 2155" << endl;
				int stophere;
				cin >> stophere;
				if (MinCostBoundaryEdgeCollapse(eid)) 
					deleteSphereNum ++;    
			} 
		} 
	}else
	{ // second condition...
		while (deleteSphereNum < threshold && numVertices > 1 && !edge_collapses_queue.empty())
		{
			//if (maxhausdorff_distance / pmesh->bb_diagonal_length >= end_multi)
			//	break;

			//if (sqrt(max_mean_squre_error) / pmesh->bb_diagonal_length >= end_multi)
			//	break;

			EdgeInfo topEdge = edge_collapses_queue.top();
			edge_collapses_queue.pop(); 
			unsigned eid = topEdge.edge_num;
			if(edges[eid].first && ValidVertex(edges[eid].second->vertices_.first) && ValidVertex(edges[eid].second->vertices_.second)
				&& (vertices[edges[eid].second->vertices_.first].second->saved_vertex == 0 || vertices[edges[eid].second->vertices_.second].second->saved_vertex == 0)
			)
			{

				cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
				cout << "Under the second condition --------------------------------------" << endl;

					
				// read the points in...

				//

				//cout << edges[eid].second->vertices_.first << endl;
				//cout << vertices[edges[eid].second->vertices_.first].first << endl;
				//cout << vertices[edges[eid].second->vertices_.first].second->sphere.center << endl;
				//if (deleteSphereNum == 36002) {
				//	break;
				//}
				if (MinCostEdgeCollapse(eid))
				{
					deleteSphereNum++;
					cout <<"hhh deleting in Simplify_Select_Pole..."<< deleteSphereNum << " / "<< vertices.size() << endl;
				}
			}

			//if (maxhausdorff_distance / pmesh->bb_diagonal_length >= start_multi)
			//	GetSavedPointNumber();

			//if (sqrt(max_mean_squre_error) / pmesh->bb_diagonal_length >= start_multi)
			//	GetSavedPointNumber();
		}
	}
}

void SlabMesh::Simplify(int threshold) {


	/*
	
	
	This method is not
	
	NOt used!
	NOt used!
	NOt used!
	NOt used!
	NOt used!
	NOt used!
	NOt used!
	NOt used!
	NOt used!
	NOt used!
	NOt used!
	NOt used!
	NOt used!
	NOt used!
	NOt used!

	*/




	vector<vector<double> > vset;

	vector<vector<int> >eset;
	vector<vector<int> > fset;



	// then we iterate in vertices set and label them...
	for (int select_id = 0;select_id < vset.size(); select_id++) {
		double min_dis = 999999;
		int min_idx = -1;
		for (int vid = 0; vid < vertices.size(); vid++) {
			vertices[vid].second->is_non_manifold;

			vertices[vid].second->sphere.center;
			double temp_dis = p_distance(vset[select_id][0],
				vset[select_id][1],
				vset[select_id][2],
				vertices[vid].second->sphere.center[0],
				vertices[vid].second->sphere.center[1],
				vertices[vid].second->sphere.center[2]);
			if (temp_dis < min_dis && vertices[vid].second->is_selected == 0) {
				min_dis = temp_dis;
				min_idx = vid;
			}
		}
		vertices[min_idx].second->is_selected = 1;
		

	}
	int check_count = 0;

	cout << "how many: " << check_count << " / " << vset.size() << endl;

	// this is for checking...
	// double checked.
	std::ofstream fsout2("C://Users//frank//Desktop//nearest_poles.obj");
	for (int vid = 0; vid < vertices.size(); vid++) {
		if (vertices[vid].second->is_selected == 1) {
			fsout2 << "v " << vertices[vid].second->sphere.center << endl;
		}
	}
	fsout2.close();



	// find the nearest one...



	// ���򻯵�С��50������ʱ������������˵�ı߽��кϲ�
	//
	if (numVertices <= 20)
	{
		if (initial_boundary_preserve == false)
		{
			initial_boundary_preserve = true;
			InitialTopologyProperty();
			for (int i = 0; i < vertices.size(); i++)
			{
				if (vertices[i].first)
				{
					set<unsigned> fir_edges = vertices[i].second->edges_;
					for (set<unsigned>::iterator si = fir_edges.begin(); si != fir_edges.end(); si++)
					{
						int index = edges[*si].second->vertices_.first == i ?
							edges[*si].second->vertices_.second : edges[*si].second->vertices_.first;

						if (vertices[index].second->edges_.size() == 1 && vertices[index].second->faces_.size() == 0)
						{
							edges[*si].second->topo_contractable = false;
						}
					}
				}
			}
		}
	}
	// This is the key simplification process.
	int deleteSphereNum = 0;
	if (!boundary_edge_collapses_queue.empty())
	{
		while (deleteSphereNum < threshold && numVertices > 1 && !boundary_edge_collapses_queue.empty())
		{
			EdgeInfo topEdge = boundary_edge_collapses_queue.top();
			boundary_edge_collapses_queue.pop();
			unsigned eid = topEdge.edge_num;
			if (edges[eid].first && ValidVertex(edges[eid].second->vertices_.first) && ValidVertex(edges[eid].second->vertices_.second))
			{
				std::cout << edges[eid].second->vertices_.first << "  " << edges[eid].second->vertices_.second << endl;
				std::cout << "Under the first condition --------------------------------------" << endl;
				int stophere;
				cin >> stophere;
				if (MinCostBoundaryEdgeCollapse(eid))
					deleteSphereNum++;
			}
		}
	}
	else
	{ // second condition...
		while (deleteSphereNum < threshold && numVertices > 1 && !edge_collapses_queue.empty())
		{
			//if (maxhausdorff_distance / pmesh->bb_diagonal_length >= end_multi)
			//	break;

			//if (sqrt(max_mean_squre_error) / pmesh->bb_diagonal_length >= end_multi)
			//	break;

			EdgeInfo topEdge = edge_collapses_queue.top();
			edge_collapses_queue.pop();
			unsigned eid = topEdge.edge_num;
			if (edges[eid].first && ValidVertex(edges[eid].second->vertices_.first) && ValidVertex(edges[eid].second->vertices_.second)
				//&& vertices[edges[eid].second->vertices_.first].second->is_selected == 0
				//&& vertices[edges[eid].second->vertices_.second].second->is_selected == 0)
				)
			{

				cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
				cout << "Under the second condition --------------------------------------" << endl;


				// read the points in...

				//

				//cout << edges[eid].second->vertices_.first << endl;
				//cout << vertices[edges[eid].second->vertices_.first].first << endl;
				//cout << vertices[edges[eid].second->vertices_.first].second->sphere.center << endl;
				//if (deleteSphereNum == 36000) {
					//break;
				//}
				if (MinCostEdgeCollapse(eid))
				{
					deleteSphereNum++;
					cout << deleteSphereNum << " / " << vertices.size() << endl;
					
				}
			}

			//if (maxhausdorff_distance / pmesh->bb_diagonal_length >= start_multi)
			//	GetSavedPointNumber();

			//if (sqrt(max_mean_squre_error) / pmesh->bb_diagonal_length >= start_multi)
			//	GetSavedPointNumber();
		}
	}
}




void SlabMesh::initCollapseQueue(){

	// first initial the edges with fake boundary edge.
	for (int i = 0; i < numEdges; i++)
	{ 
		if (edges[i].first)
		{
			EvaluateEdgeCollapseCost(i);
			edge_collapses_queue.push(EdgeInfo(i, edges[i].second->collapse_cost));
		}
	}
}

void SlabMesh::initBoundaryCollapseQueue()
{
	for (int i = 0; i < edges.size(); i ++)
	{
		if (edges[i].first)
		{
			unsigned fir = edges[i].second->vertices_.first;
			unsigned sec = edges[i].second->vertices_.second;

			switch(boundary_compute_scale)
			{
			case 1:
				if (!vertices[fir].second->fake_boundary_vertex || !vertices[sec].second->fake_boundary_vertex)
					continue;
				break;
			case 2:
				if (!vertices[fir].second->fake_boundary_vertex && !vertices[sec].second->fake_boundary_vertex)
					//if (vertices[fir].second->boundary_edge_vec.size() < 2 && vertices[sec].second->boundary_edge_vec.size() < 2)
					continue;
				break;
			case 3:
				if (!vertices[fir].second->fake_boundary_vertex && !vertices[sec].second->fake_boundary_vertex)
					continue;
				break;
			default:
				break;
			}

			//EvaluateEdgeCollapseCost(i);
			EvaluateEdgeHausdorffCost(i);
			boundary_edge_collapses_queue.push(EdgeInfo(i, edges[i].second->collapse_cost));
		}
	}
}

double SlabMesh::NearestPoint(Vector3d point, unsigned vid)
{
	set<unsigned> near_faces = vertices[vid].second->faces_;
	set<unsigned> near_edges = vertices[vid].second->edges_;

	double mind = DBL_MAX;
	// calculation of related faces
	for (set<unsigned>::iterator si = near_faces.begin(); si != near_faces.end(); si++)
	{
		if (!faces[*si].first)
			continue;
		SlabFace sf = *faces[*si].second;
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
	for (set<unsigned>::iterator si = near_edges.begin(); si != near_edges.end(); si++)
	{
		if (!edges[*si].first)
			continue;
		SlabEdge se = *edges[*si].second;
		if (se.valid_cone == false)
			continue;

		SlabVertex v[2];
		Vector3d tfp;
		double td;
		double tr;
		v[0] = *vertices[se.vertices_.first].second;
		v[1] = *vertices[se.vertices_.second].second;
		Vector3d v0 = v[0].sphere.center;
		Vector3d v1 = v[1].sphere.center;
		double t((point-v0).Dot(v1-v0) / (v1-v0).SquaredLength());
		if( (t >= 0.0) && (t <= 1.0) )
		{
			tfp = (1.0-t)*v0 + t*v1;
			td = (point-tfp).Length();
			tr = (1.0-t)*v[0].sphere.radius + t*v[1].sphere.radius;

			if (abs(td - tr) < mind) mind = abs(td - tr);
		}
	}

	return mind;
}

// simple method, do not add any plane to preserve the boundary 
void SlabMesh::PreservBoundaryMethodOne()
{
	// �����е�boundary_edge���ϱ߽籣��
	for (unsigned i = 0; i < edges.size(); i++)
	{
		if (!edges[i].first || !edges[i].second->fake_boundary_edge)
			continue;

		SlabEdge se = *(edges[i].second);
		if (edges[i].second->faces_.size() == 0)
		{
			unsigned ver_index[2];
			ver_index[0] = se.vertices_.first;
			ver_index[1] = se.vertices_.second;
			Sphere ver[3];
			ver[0]= vertices[ver_index[0]].second->sphere;
			ver[1] = vertices[ver_index[1]].second->sphere;
			ver[2].center = (ver[0].center + ver[1].center);
			ver[2].radius = (ver[0].radius + ver[1].radius) / 2.0;

			SimpleTriangle st[2];
			Wm4::Vector3d pos[3];
			double radius[3];
			for(int count = 0; count < 3; count++)
			{
				pos[count] = ver[count].center;
				radius[count] = ver[count].radius;
			}
			if(TriangleFromThreeSpheres(pos[0],radius[0],pos[1],radius[1],pos[2],radius[2],st[0],st[1]))
			{
				if (st[0].normal == Vector3d(0., 0., 0.) || st[1].normal == Vector3d(0., 0., 0.))
					continue;

				// �ӵ�һ��slab�е�����ƽ��
				Vector4d normal1(st[0].normal.X(), st[0].normal.Y(), st[0].normal.Z(), 1.0);
				Vector4d normal2(st[1].normal.X(), st[1].normal.Y(), st[1].normal.Z(), 1.0);
				// compute the matrix of A
				Matrix4d temp_A1, temp_A2;
				temp_A1.MakeTensorProduct(normal1, normal1);
				temp_A2.MakeTensorProduct(normal2, normal2);
				temp_A1 *= 2.0;
				temp_A2 *= 2.0;
				for (int i = 0; i < 2; i++)
				{
					Vector4d C1(ver[i].center.X(), ver[i].center.Y(), ver[i].center.Z(), ver[i].radius);
					// compute the matrix of b
					double normal_mul_point1 = normal1.Dot(C1);
					double normal_mul_point2 = normal2.Dot(C1);
					Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
					Wm4::Vector4d temp_b2 = normal2 * 2 * normal_mul_point2;
					//compute c
					double temp_c1 = normal_mul_point1 * normal_mul_point1;
					double temp_c2 = normal_mul_point2 * normal_mul_point2;
					vertices[ver_index[i]].second->add_A += temp_A1;
					vertices[ver_index[i]].second->add_A += temp_A2;
					vertices[ver_index[i]].second->add_b += temp_b1;
					vertices[ver_index[i]].second->add_b += temp_b2;
					vertices[ver_index[i]].second->add_c += temp_c1;
					vertices[ver_index[i]].second->add_c += temp_c2;
				}

				// �ӵڶ���slab�е�����ƽ��
				Vector3d ver1_to_ver2 = ver[0].center - ver[1].center;
				Vector3d t1 = ver1_to_ver2.Cross(st[0].normal);
				Vector3d t2 = ver1_to_ver2.Cross(st[1].normal);
				Vector4d tnormal1(t1.X(), t1.Y(), t1.Z(), 1.0);
				Vector4d tnormal2(t2.X(), t2.Y(), t2.Z(), 1.0);
				// compute the matrix of A
				temp_A1.MakeTensorProduct(tnormal1, tnormal1);
				temp_A2.MakeTensorProduct(tnormal2, tnormal2);
				temp_A1 *= 2.0;
				temp_A2 *= 2.0;
				for (int i = 0; i < 2; i++)
				{
					Vector4d C1(ver[i].center.X(), ver[i].center.Y(), ver[i].center.Z(), ver[i].radius);
					// compute the matrix of b
					double normal_mul_point1 = tnormal1.Dot(C1);
					double normal_mul_point2 = tnormal2.Dot(C1);
					Wm4::Vector4d temp_b1 = tnormal1 * 2 * normal_mul_point1;
					Wm4::Vector4d temp_b2 = tnormal2 * 2 * normal_mul_point2;
					//compute c
					double temp_c1 = normal_mul_point1 * normal_mul_point1;
					double temp_c2 = normal_mul_point2 * normal_mul_point2;
					vertices[ver_index[i]].second->add_A += temp_A1;
					vertices[ver_index[i]].second->add_A += temp_A2;
					vertices[ver_index[i]].second->add_b += temp_b1;
					vertices[ver_index[i]].second->add_b += temp_b2;
					vertices[ver_index[i]].second->add_c += temp_c1;
					vertices[ver_index[i]].second->add_c += temp_c2;
				}
			}

			//// ���ڱ�¶�����ĵ��ټ�һ������ƽ��
			//if (vertices[ver_index[0]].second->edges_.size() == 1)
			//{
			//	Vector3d ver1_to_ver2 = ver[0].center - ver[1].center;
			//	Vector4d normal1(ver1_to_ver2.X(), ver1_to_ver2.Y(), ver1_to_ver2.Z(), 1.0);
			//	// compute the matrix of A
			//	Matrix4d temp_A1;
			//	temp_A1.MakeTensorProduct(normal1, normal1);
			//	temp_A1 *= 2.0;
			//	Vector4d C1(ver[0].center.X(), ver[0].center.Y(), ver[0].center.Z(), ver[0].radius);
			//	// compute the matrix of b
			//	double normal_mul_point1 = normal1.Dot(C1);
			//	Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
			//	//compute c
			//	double temp_c1 = normal_mul_point1 * normal_mul_point1;
			//	vertices[ver_index[0]].second->add_A += temp_A1 * 100.0;
			//	vertices[ver_index[0]].second->add_b += temp_b1 * 100.0;
			//	vertices[ver_index[0]].second->add_c += temp_c1 * 100.0;
			//}
			//if (vertices[ver_index[1]].second->edges_.size() == 1)
			//{
			//	Vector3d ver2_to_ver1 = ver[1].center - ver[0].center;
			//	Vector4d normal1(ver2_to_ver1.X(), ver2_to_ver1.Y(), ver2_to_ver1.Z(), 1.0);
			//	// compute the matrix of A
			//	Matrix4d temp_A1;
			//	temp_A1.MakeTensorProduct(normal1, normal1);
			//	temp_A1 *= 2.0;
			//	Vector4d C1(ver[1].center.X(), ver[1].center.Y(), ver[1].center.Z(), ver[1].radius);
			//	// compute the matrix of b
			//	double normal_mul_point1 = normal1.Dot(C1);
			//	Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
			//	//compute c
			//	double temp_c1 = normal_mul_point1 * normal_mul_point1;
			//	vertices[ver_index[1]].second->add_A += temp_A1 * 100.0;
			//	vertices[ver_index[1]].second->add_b += temp_b1 * 100.0;
			//	vertices[ver_index[1]].second->add_c += temp_c1 * 100.0;
			//}
		}
	}
}

// just add one normal to all the Boundary Vertex Sphere with two Boundary edges.
void SlabMesh::PreservBoundaryMethodTwo()
{
	// for each boundary vertex sphere, add a normal
	for (unsigned i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].first && vertices[i].second->boundary_edge_vec.size() == 2)
		{
			Vector3d add_normal(0.0, 0.0, 0.0);
			int valid_edge = 0;
			for (set<unsigned>::iterator vi = vertices[i].second->boundary_edge_vec.begin(); vi != vertices[i].second->boundary_edge_vec.end(); vi++)
			{
				if (edges[*vi].second->faces_.size() == 0)
					continue;
				unsigned face_num = *(edges[*vi].second->faces_.begin());			
				SlabFace sf = *faces[face_num].second;
				if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) || 
					sf.st[1].normal == Vector3d(0., 0., 0.))
					continue;

				valid_edge++;
				unsigned ver_index = edges[*vi].second->vertices_.first == i ? edges[*vi].second->vertices_.second 
					: edges[*vi].second->vertices_.first;
				Vector3d temp_norm = vertices[i].second->sphere.center - vertices[ver_index].second->sphere.center;
				temp_norm.Normalize();
				add_normal += temp_norm;
			}
			if (valid_edge == 0)
				continue;

			add_normal /= valid_edge;
			add_normal.Normalize();

			vertices[i].second->boundVec = add_normal;

			// ȷ�����ӵ����Ȩ�ش�С
			//Vector3d boudary_vec[2];
			//for (int index = 0; index < 2; index++)
			//{
			//	set<unsigned>::iterator vi = vertices[i].second->boundary_edge_vec.begin() + index;
			//	unsigned ver_index = edges[*vi].second->vertices_.first == i ? edges[*vi].second->vertices_.second 
			//		: edges[*vi].second->vertices_.first;
			//	Vector3d temp_norm = vertices[i].second->sphere.center - vertices[ver_index].second->sphere.center;
			//	boudary_vec[index] = temp_norm;
			//}
			//vertices[i].second->collaspe_weight = sin(VectorAngle(boudary_vec[0], boudary_vec[1]));

			// �ж�����߽��������͹�����ĵ㻹�ǰ���ȥ�ĵ�
			bool boundary_vertex = false;
			for (auto si = vertices[i].second->edges_.begin(); si != vertices[i].second->edges_.end(); si++)
			{
				unsigned temp_ind = edges[*si].second->vertices_.first == i ? edges[*si].second->vertices_.second 
					: edges[*si].second->vertices_.first;
				Vector3d temp_vec = vertices[temp_ind].second->sphere.center - vertices[i].second->sphere.center;
				double temp_angle = acos(temp_vec.Dot(add_normal) / temp_vec.Length());
				boundary_vertex = temp_angle < Wm4::Math<double>::PI / 2.0 ? true : false;
				if (boundary_vertex == true)
				{
					vertices[i].second->boundary_vertex = true;
					break;
				}
			}

			if (boundary_vertex == false)
			{
				Vector4d normal(add_normal.X(), add_normal.Y(), add_normal.Z(), 1.0);
				Matrix4d temp_A;
				temp_A.MakeTensorProduct(normal, normal);
				temp_A *= 2.0;
				Sphere ve = vertices[i].second->sphere;
				Vector4d C1(ve.center.X(), ve.center.Y(), ve.center.Z(), ve.radius);
				double normal_mul_point = normal.Dot(C1);
				Wm4::Vector4d temp_b = normal * 2 * normal_mul_point;
				double temp_c = normal_mul_point * normal_mul_point;

				vertices[i].second->add_A += temp_A;
				vertices[i].second->add_b += temp_b;
				vertices[i].second->add_c += temp_c;

				vertices[i].second->related_face ++;
			}
			else
			{
				add_normal = add_normal * (-1.0);
				Vector4d normal(add_normal.X(), add_normal.Y(), add_normal.Z(), 1.0);
				Matrix4d temp_A;
				temp_A.MakeTensorProduct(normal, normal);
				temp_A *= 2.0;
				Sphere ve = vertices[i].second->sphere;
				Vector4d C1(ve.center.X(), ve.center.Y(), ve.center.Z(), ve.radius);
				double normal_mul_point = normal.Dot(C1);
				Wm4::Vector4d temp_b = normal * 2 * normal_mul_point;
				double temp_c = normal_mul_point * normal_mul_point;

				vertices[i].second->add_A += temp_A;
				vertices[i].second->add_b += temp_b;
				vertices[i].second->add_c += temp_c;

				vertices[i].second->related_face ++;
			}
		}
	}
}

// add a plane to the Possible Spike Vertex, and add a slab to the Possible Boundary Vertex.
void SlabMesh::PreservBoundaryMethodThree()
{
	// for each boundary vertex sphere, add a normal
	for (unsigned i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].first && vertices[i].second->boundary_edge_vec.size() == 2)
		{
			Vector3d add_normal(0.0, 0.0, 0.0);
			int valid_edge = 0;
			for (set<unsigned>::iterator vi = vertices[i].second->boundary_edge_vec.begin(); vi != vertices[i].second->boundary_edge_vec.end(); vi++)
			{
				if (edges[*vi].second->faces_.size() == 0)
					continue;
				unsigned face_num = *(edges[*vi].second->faces_.begin());			
				SlabFace sf = *faces[face_num].second;
				if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) || 
					sf.st[1].normal == Vector3d(0., 0., 0.))
					continue;

				valid_edge++;
				unsigned ver_index = edges[*vi].second->vertices_.first == i ? edges[*vi].second->vertices_.second 
					: edges[*vi].second->vertices_.first;
				Vector3d temp_norm = vertices[i].second->sphere.center - vertices[ver_index].second->sphere.center;
				temp_norm.Normalize();
				add_normal += temp_norm;
			}
			if (valid_edge == 0)
				continue;

			add_normal /= valid_edge;
			add_normal.Normalize();

			vertices[i].second->boundVec = add_normal;

			// ȷ�����ӵ����Ȩ�ش�С
			//Vector3d boudary_vec[2];
			//for (int index = 0; index < 2; index++)
			//{
			//	auto vi = vertices[i].second->boundary_edge_vec.begin() + index;
			//	unsigned ver_index = edges[*vi].second->vertices_.first == i ? edges[*vi].second->vertices_.second 
			//		: edges[*vi].second->vertices_.first;
			//	Vector3d temp_norm = vertices[i].second->sphere.center - vertices[ver_index].second->sphere.center;
			//	boudary_vec[index] = temp_norm;
			//}
			//vertices[i].second->collaspe_weight = sin(VectorAngle(boudary_vec[0], boudary_vec[1]));

			// �ж�����߽��������͹�����ĵ㻹�ǰ���ȥ�ĵ�
			bool boundary_vertex = false;
			for (auto si = vertices[i].second->edges_.begin(); si != vertices[i].second->edges_.end(); si++)
			{
				unsigned temp_ind = edges[*si].second->vertices_.first == i ? edges[*si].second->vertices_.second 
					: edges[*si].second->vertices_.first;
				Vector3d temp_vec = vertices[temp_ind].second->sphere.center - vertices[i].second->sphere.center;
				double temp_angle = acos(temp_vec.Dot(add_normal) / temp_vec.Length());
				boundary_vertex = temp_angle < Wm4::Math<double>::PI / 2.0 ? true : false;
				if (boundary_vertex == true)
				{
					vertices[i].second->boundary_vertex = true;
					break;
				}
			}

			if (boundary_vertex == false)
			{
				Vector4d normal(add_normal.X(), add_normal.Y(), add_normal.Z(), 1.0);
				Matrix4d temp_A;
				temp_A.MakeTensorProduct(normal, normal);
				temp_A *= 2.0;
				Sphere ve = vertices[i].second->sphere;
				Vector4d C1(ve.center.X(), ve.center.Y(), ve.center.Z(), ve.radius);
				double normal_mul_point = normal.Dot(C1);
				Wm4::Vector4d temp_b = normal * 2 * normal_mul_point;
				double temp_c = normal_mul_point * normal_mul_point;

				vertices[i].second->add_A += temp_A;
				vertices[i].second->add_b += temp_b;
				vertices[i].second->add_c += temp_c;

				vertices[i].second->related_face ++;
			}
		}
	}

	// �����е�boundary_edge���ϱ߽籣��
	for (unsigned i = 0; i < edges.size(); i++)
	{
		if (edges[i].first && edges[i].second->faces_.size() == 1)
		{
			unsigned face_num = *(edges[i].second->faces_.begin());

			unsigned ver_index[2];
			ver_index[0] = edges[i].second->vertices_.first;
			ver_index[1] = edges[i].second->vertices_.second;
			Sphere ver[2];
			ver[0]= vertices[ver_index[0]].second->sphere;
			ver[1] = vertices[ver_index[1]].second->sphere;
			Vector3d ver1_to_ver2 = ver[0].center - ver[1].center;

			SlabFace sf = *faces[face_num].second;
			if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) || 
				sf.st[1].normal == Vector3d(0., 0., 0.))
				continue;

			Vector3d temp_normal1(sf.st[0].normal.X(), sf.st[0].normal.Y(), sf.st[0].normal.Z());
			Vector3d temp_normal2(sf.st[1].normal.X(), sf.st[1].normal.Y(), sf.st[1].normal.Z());
			Vector3d t1 = ver1_to_ver2.Cross(temp_normal1);
			Vector3d t2 = ver1_to_ver2.Cross(temp_normal2);

			// ����boundary_edge����һ��slab���б߽籣��
			Vector4d normal1(t1.X(), t1.Y(), t1.Z(), 1.0);
			Vector4d normal2(t2.X(), t2.Y(), t2.Z(), 1.0);

			// compute the matrix of A
			Matrix4d temp_A1, temp_A2;
			temp_A1.MakeTensorProduct(normal1, normal1);
			temp_A2.MakeTensorProduct(normal2, normal2);
			temp_A1 *= 2.0;
			temp_A2 *= 2.0;

			for (int i = 0; i < 2; i++)
			{	
				//if (vertices[ver_index[i]].second->boundary_vertex == false)
				//	continue;

				Vector4d C1(ver[i].center.X(), ver[i].center.Y(), ver[i].center.Z(), ver[i].radius);

				// compute the matrix of b
				double normal_mul_point1 = normal1.Dot(C1);
				double normal_mul_point2 = normal2.Dot(C1);
				Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
				Wm4::Vector4d temp_b2 = normal2 * 2 * normal_mul_point2;

				//compute c
				double temp_c1 = normal_mul_point1 * normal_mul_point1;
				double temp_c2 = normal_mul_point2 * normal_mul_point2;

				vertices[ver_index[i]].second->add_A += temp_A1;
				vertices[ver_index[i]].second->add_A += temp_A2;
				vertices[ver_index[i]].second->add_b += temp_b1;
				vertices[ver_index[i]].second->add_b += temp_b2;
				vertices[ver_index[i]].second->add_c += temp_c1;
				vertices[ver_index[i]].second->add_c += temp_c2;
			}
		}
	}
}

void SlabMesh::PreservBoundaryMethodFour()
{
	// �����е�boundary_edge���ϱ߽籣��
	for (unsigned i = 0; i < edges.size(); i++)
	{
		if (!edges[i].first || !edges[i].second->fake_boundary_edge)
			continue;

		SlabEdge se = *(edges[i].second);
		if (se.faces_.size() == 1)
		{
			unsigned face_num = *(se.faces_.begin());

			set<unsigned> sv = faces[face_num].second->vertices_;
			Vector3d face_normal = faces[face_num].second->normal;

			unsigned ver_index[3];
			ver_index[0] = se.vertices_.first;
			sv.erase(ver_index[0]);
			ver_index[1] = se.vertices_.second;
			sv.erase(ver_index[1]);
			ver_index[2] = *(sv.begin());
			Sphere ver[2];
			ver[0]= vertices[ver_index[0]].second->sphere;
			ver[1] = vertices[ver_index[1]].second->sphere;

			Vector3d v1v2 = vertices[ver_index[0]].second->sphere.center - vertices[ver_index[1]].second->sphere.center;
			Vector3d v1v3 = vertices[ver_index[2]].second->sphere.center - vertices[ver_index[0]].second->sphere.center;
			Vector3d temp_nor = face_normal.Cross(v1v2);

			double temp_angle = acos(temp_nor.Dot(v1v3) / temp_nor.Length() / v1v3.Length());
			bool dir = temp_angle > Wm4::Math<double>::PI / 2.0 ? true : false;

			if (dir == false)
				temp_nor *= -1;

			// ����boundary_edge����һ��ƽ����б߽籣��
			Vector4d normal1(temp_nor.X(), temp_nor.Y(), temp_nor.Z(), 1.0);
			// compute the matrix of A
			Matrix4d temp_A1;
			temp_A1.MakeTensorProduct(normal1, normal1);
			temp_A1 *= 2.0;

			// �Բ�ͬratio�ı߽���ӳ�䴦��С��0.2�Ĳ�������(0.2,1)ӳ�䵽(2, 10)
			double ratio = GetRatioHyperbolicEuclid(i);
			double w1 = 1.0;
			//w1 = 0.02 * ratio * ratio * ratio * ratio * ratio * ratio;
			w1 =  0.1 * ratio * ratio;
			//w1 = 3 * ratio * ratio;
			//w1 = 1 * ratio * ratio * ratio;
			//if (ratio >= 0.3)
			//{
			//	w1 = (ratio - 0.3) / 0.7 * (3 - 1) + 1;
			//}

			for (int i = 0; i < 2; i++)
			{	
				Vector4d C1(ver[i].center.X(), ver[i].center.Y(), ver[i].center.Z(), ver[i].radius);

				// compute the matrix of b
				double normal_mul_point1 = normal1.Dot(C1);
				Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;

				//compute c
				double temp_c1 = normal_mul_point1 * normal_mul_point1;

				vertices[ver_index[i]].second->add_A += temp_A1 * bound_weight * w1;
				vertices[ver_index[i]].second->add_b += temp_b1 * bound_weight * w1;
				vertices[ver_index[i]].second->add_c += temp_c1 * bound_weight * w1;
			}
		}else if (edges[i].second->faces_.size() == 0)
		{
			unsigned ver_index[2];
			ver_index[0] = se.vertices_.first;
			ver_index[1] = se.vertices_.second;
			Sphere ver[3];
			ver[0]= vertices[ver_index[0]].second->sphere;
			ver[1] = vertices[ver_index[1]].second->sphere;
			ver[2].center = (ver[0].center + ver[1].center);
			ver[2].radius = (ver[0].radius + ver[1].radius) / 2.0;

			SimpleTriangle st[2];
			Wm4::Vector3d pos[3];
			double radius[3];
			for(int count = 0; count < 3; count++)
			{
				pos[count] = ver[count].center;
				radius[count] = ver[count].radius;
			}
			if(TriangleFromThreeSpheres(pos[0],radius[0],pos[1],radius[1],pos[2],radius[2],st[0],st[1]))
			{
				if (st[0].normal == Vector3d(0., 0., 0.) || st[1].normal == Vector3d(0., 0., 0.))
					continue;

				// �ӵ�һ��slab�е�����ƽ��
				Vector4d normal1(st[0].normal.X(), st[0].normal.Y(), st[0].normal.Z(), 1.0);
				Vector4d normal2(st[1].normal.X(), st[1].normal.Y(), st[1].normal.Z(), 1.0);
				// compute the matrix of A
				Matrix4d temp_A1, temp_A2;
				temp_A1.MakeTensorProduct(normal1, normal1);
				temp_A2.MakeTensorProduct(normal2, normal2);
				temp_A1 *= 2.0;
				temp_A2 *= 2.0;
				for (int i = 0; i < 2; i++)
				{
					Vector4d C1(ver[i].center.X(), ver[i].center.Y(), ver[i].center.Z(), ver[i].radius);
					// compute the matrix of b
					double normal_mul_point1 = normal1.Dot(C1);
					double normal_mul_point2 = normal2.Dot(C1);
					Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
					Wm4::Vector4d temp_b2 = normal2 * 2 * normal_mul_point2;
					//compute c
					double temp_c1 = normal_mul_point1 * normal_mul_point1;
					double temp_c2 = normal_mul_point2 * normal_mul_point2;
					vertices[ver_index[i]].second->add_A += temp_A1;
					vertices[ver_index[i]].second->add_A += temp_A2;
					vertices[ver_index[i]].second->add_b += temp_b1;
					vertices[ver_index[i]].second->add_b += temp_b2;
					vertices[ver_index[i]].second->add_c += temp_c1;
					vertices[ver_index[i]].second->add_c += temp_c2;
				}

				// �ӵڶ���slab�е�����ƽ��
				Vector3d ver1_to_ver2 = ver[0].center - ver[1].center;
				Vector3d t1 = ver1_to_ver2.Cross(st[0].normal);
				Vector3d t2 = ver1_to_ver2.Cross(st[1].normal);
				Vector4d tnormal1(t1.X(), t1.Y(), t1.Z(), 1.0);
				Vector4d tnormal2(t2.X(), t2.Y(), t2.Z(), 1.0);
				// compute the matrix of A
				temp_A1.MakeTensorProduct(tnormal1, tnormal1);
				temp_A2.MakeTensorProduct(tnormal2, tnormal2);
				temp_A1 *= 2.0;
				temp_A2 *= 2.0;
				for (int i = 0; i < 2; i++)
				{
					Vector4d C1(ver[i].center.X(), ver[i].center.Y(), ver[i].center.Z(), ver[i].radius);
					// compute the matrix of b
					double normal_mul_point1 = tnormal1.Dot(C1);
					double normal_mul_point2 = tnormal2.Dot(C1);
					Wm4::Vector4d temp_b1 = tnormal1 * 2 * normal_mul_point1;
					Wm4::Vector4d temp_b2 = tnormal2 * 2 * normal_mul_point2;
					//compute c
					double temp_c1 = normal_mul_point1 * normal_mul_point1;
					double temp_c2 = normal_mul_point2 * normal_mul_point2;
					vertices[ver_index[i]].second->add_A += temp_A1;
					vertices[ver_index[i]].second->add_A += temp_A2;
					vertices[ver_index[i]].second->add_b += temp_b1;
					vertices[ver_index[i]].second->add_b += temp_b2;
					vertices[ver_index[i]].second->add_c += temp_c1;
					vertices[ver_index[i]].second->add_c += temp_c2;
				}
			}

			double ratio = GetRatioHyperbolicEuclid(i);
			double w1 = 1.0;
			//w1 = 0.02 * ratio * ratio * ratio * ratio * ratio * ratio;
			w1 = 0.1 * ratio * ratio;
			//w1 = 1 * ratio * ratio * ratio * ratio * ratio * ratio;
			//w1 = 1 * ratio * ratio * ratio;
			// ���ڱ�¶�����ĵ��ټ�һ������ƽ��
			if (vertices[ver_index[0]].second->edges_.size() == 1)
			{

				Vector3d ver1_to_ver2 = ver[0].center - ver[1].center;
				Vector4d normal1(ver1_to_ver2.X(), ver1_to_ver2.Y(), ver1_to_ver2.Z(), 1.0);
				// compute the matrix of A
				Matrix4d temp_A1;
				temp_A1.MakeTensorProduct(normal1, normal1);
				temp_A1 *= 2.0;
				Vector4d C1(ver[0].center.X(), ver[0].center.Y(), ver[0].center.Z(), ver[0].radius);
				// compute the matrix of b
				double normal_mul_point1 = normal1.Dot(C1);
				Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
				//compute c
				double temp_c1 = normal_mul_point1 * normal_mul_point1;
				//vertices[ver_index[0]].second->add_A += (temp_A1 * 100.0);
				//vertices[ver_index[0]].second->add_b += (temp_b1 * 100.0);
				//vertices[ver_index[0]].second->add_c += (temp_c1 * 100.0);
				vertices[ver_index[0]].second->add_A += temp_A1 * bound_weight * w1;
				vertices[ver_index[0]].second->add_b += temp_b1 * bound_weight * w1;
				vertices[ver_index[0]].second->add_c += temp_c1 * bound_weight * w1;
			}
			if (vertices[ver_index[1]].second->edges_.size() == 1)
			{
				Vector3d ver2_to_ver1 = ver[1].center - ver[0].center;
				Vector4d normal1(ver2_to_ver1.X(), ver2_to_ver1.Y(), ver2_to_ver1.Z(), 1.0);
				// compute the matrix of A
				Matrix4d temp_A1;
				temp_A1.MakeTensorProduct(normal1, normal1);
				temp_A1 *= 2.0;
				Vector4d C1(ver[1].center.X(), ver[1].center.Y(), ver[1].center.Z(), ver[1].radius);
				// compute the matrix of b
				double normal_mul_point1 = normal1.Dot(C1);
				Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
				//compute c
				double temp_c1 = normal_mul_point1 * normal_mul_point1;
				//vertices[ver_index[1]].second->add_A += (temp_A1 * 100.0);
				//vertices[ver_index[1]].second->add_b += (temp_b1 * 100.0);
				//vertices[ver_index[1]].second->add_c += (temp_c1 * 100.0);
				vertices[ver_index[1]].second->add_A += temp_A1 * bound_weight * w1;
				vertices[ver_index[1]].second->add_b += temp_b1 * bound_weight * w1;
				vertices[ver_index[1]].second->add_c += temp_c1 * bound_weight * w1;
			}
		}
	} 
}

void SlabMesh::clear()
{
	for (unsigned i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].first)
		{
			vertices[i].second->slab_A.MakeZero();
			vertices[i].second->add_A.MakeZero();
			vertices[i].second->slab_b = Wm4::Vector4<double>::ZERO;
			vertices[i].second->add_b = Wm4::Vector4<double>::ZERO;
			vertices[i].second->slab_c = 0.0;
			vertices[i].second->add_c = 0.0;
			vertices[i].second->fake_boundary_vertex = false;
			vertices[i].second->boundary_vertex = false;
			vertices[i].second->boundary_edge_vec.clear();
			vertices[i].second->mean_square_error = 0.0;
			vertices[i].second->related_face = 0;
		}
	}

	for (unsigned i = 0; i < edges.size(); i++)
	{
		if (edges[i].first)
		{
			edges[i].second->slab_A.MakeZero();
			edges[i].second->slab_b = Wm4::Vector4<double>::ZERO;
			edges[i].second->slab_c = 0.0;
			edges[i].second->fake_boundary_edge = false;
			edges[i].second->non_manifold_edge = false;
		}
	}

	max_mean_squre_error = 0.0;

	DistinguishVertexType();
}

void SlabMesh::RecomputerVertexType()
{
	for (unsigned i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].first)
		{
			vertices[i].second->fake_boundary_vertex = false;
			vertices[i].second->boundary_vertex = false;
			vertices[i].second->boundary_edge_vec.clear();
		}
	}

	for (unsigned i = 0; i < edges.size(); i++)
	{
		if (edges[i].first)
		{
			edges[i].second->fake_boundary_edge = false;
			edges[i].second->non_manifold_edge = false;
		}
	}

	DistinguishVertexType();
}

void SlabMesh::computebb()
{
	m_min[0] = 1e20;
	m_min[1] = 1e20;
	m_min[2] = 1e20;
	m_max[0] = -1e20;
	m_max[1] = -1e20;
	m_max[2] = -1e20;
	//check the size of the model...

	for (unsigned i = 0; i < vertices.size(); i++)
	{
		if (!vertices[i].first)
			continue;

		Vector3d ver = vertices[i].second->sphere.center;




		if (ver[0] < m_min[0])
			m_min[0] = ver[0]; 
		if (ver[1] < m_min[1])
			m_min[1] = ver[1]; 
		if (ver[2] < m_min[2])
			m_min[2] = ver[2];

		if (ver[0] > m_max[0])
			m_max[0] = ver[0]; 
		if (ver[1] > m_max[1])
			m_max[1] = ver[1]; 
		if (ver[2] > m_max[2])
			m_max[2] = ver[2];
	}
}

void SlabMesh::GetEnvelopeSet(const Vector4d & lamder, const set<unsigned> & neighbor_v, const set< std::set<unsigned> > & adj_faces, vector<Sphere> & sph_vec, vector<Cone> & con_vec, vector<SimpleTriangle> & st_vec)
{
	Vector3d ps(lamder.X(), lamder.Y(), lamder.Z());
	double rs = lamder.W();

	sph_vec.push_back(Sphere(ps,rs));

	for(std::set<unsigned>::iterator si = neighbor_v.begin(); si != neighbor_v.end(); si ++)
	{
		//sph_vec.push_back(Sphere(vertices[*si].second->pos, vertices[*si].second->radius));
		Cone newc = Cone(ps,rs,vertices[*si].second->sphere.center, vertices[*si].second->sphere.radius);
		if(newc.type != 1)
			con_vec.push_back(newc);
	}

	for(std::set< std::set<unsigned> >::iterator si = adj_faces.begin(); si != adj_faces.end(); si ++)
	{
		Vector3d cen[2];
		double rad[2];
		unsigned count = 0;
		for(std::set<unsigned>::iterator si2 = (*si).begin(); si2 != (*si).end(); si2 ++, count ++)
		{
			cen[count] = vertices[*si2].second->sphere.center;
			rad[count] = vertices[*si2].second->sphere.radius;
		}

		SimpleTriangle st[2];
		if(TriangleFromThreeSpheres(cen[0], rad[0], cen[1], rad[1], ps, rs, st[0], st[1]))
		{
			st_vec.push_back(st[0]);
			st_vec.push_back(st[1]);
		}

	}

	return;
}


double SlabMesh::EvaluateVertexDistanceErrorEnvelope(Vector4d & lamdar, set<unsigned> & neighbor_vertices, set< set<unsigned> > & neighbor_faces, set<unsigned> & bplist)
{

    // this is for calculating the error of recon.
	bool valid_cone = true;
	vector<Sphere> sph_vec;
	vector<Cone> con_vec;
	vector<SimpleTriangle> st_vec;
	GetEnvelopeSet(lamdar, neighbor_vertices, neighbor_faces, sph_vec, con_vec, st_vec);

	double maxerror(0.0);
	for(set<unsigned>::iterator si = bplist.begin(); si != bplist.end(); si ++)
	{
		Vector3d p(pmesh->pVertexList[*si]->point()[0], pmesh->pVertexList[*si]->point()[1], pmesh->pVertexList[*si]->point()[2]);
		double tempdist;
		Vector3d tempfp;

		double mindist(1e20);
		for(unsigned i = 0; i < sph_vec.size(); i ++)
		{
			sph_vec[i].ProjectOntoSphere(p,tempfp,tempdist);
			tempdist = fabs(tempdist);
			mindist = min(tempdist,mindist);
		}

		tempdist = abs((p - Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z())).Length() - lamdar.W());
		mindist = min(tempdist,mindist);

		for(unsigned i = 0; i < con_vec.size(); i ++)
		{
			con_vec[i].ProjectOntoCone(p,tempfp,tempdist);
			//if(tempdist < -2.*error_threshold)
			//	valid_cone = false;
			if(tempdist < 0)
				tempdist = -tempdist;
			mindist = min(tempdist,mindist);
		}

		for(unsigned i = 0; i < st_vec.size(); i ++)
		{
			st_vec[i].ProjectOntoSimpleTriangle(p,tempfp,tempdist);
			mindist = min(tempdist, mindist);
		}
		maxerror = max(maxerror, mindist);
	}
	//vertices[vid].second->v_evaluated_distance_error_envelope = maxerror;
	return maxerror;
}

double SlabMesh::GetHyperbolicLength(unsigned eid)
{
	double hyperbolic_weight;
	unsigned v1 = edges[eid].second->vertices_.first;
	unsigned v2 = edges[eid].second->vertices_.second;
	double edge_length = Vector3d(vertices[v1].second->sphere.center - vertices[v2].second->sphere.center).Length();
	double r1 = vertices[v1].second->sphere.radius;
	double r2 = vertices[v2].second->sphere.radius;
	if (r1 <= r2)
		hyperbolic_weight = edge_length - (r2 - r1);
	else
		hyperbolic_weight = edge_length - (r1 - r2);
	hyperbolic_weight = max(hyperbolic_weight, 0.0); 
	return hyperbolic_weight;
}

double SlabMesh::GetRatioHyperbolicEuclid(unsigned eid)
{
	double hyperbolic_distance;
	unsigned v1 = edges[eid].second->vertices_.first;
	unsigned v2 = edges[eid].second->vertices_.second;
	double edge_length = Vector3d(vertices[v1].second->sphere.center - vertices[v2].second->sphere.center).Length();
	double r1 = vertices[v1].second->sphere.radius;
	double r2 = vertices[v2].second->sphere.radius;
	if (r1 <= r2)
		hyperbolic_distance = edge_length - (r2 - r1);
	else
		hyperbolic_distance = edge_length - (r1 - r2);
	hyperbolic_distance = max(hyperbolic_distance, 0.0); 

	if (edge_length == 0.0)
		return 0.0;

	return hyperbolic_distance / edge_length;
}

void SlabMesh::ExportSimplifyResult()
{
	//std::ofstream f_result_out;
	//f_result_out.open("Result.txt", ios::app);

	//f_result_out << simplified_boundary_edges << "\t" << simplified_inside_edges << "\t" << maxhausdorff_distance << endl;
}


void SlabMesh::Export(std::string fname){
	fname += "___v_";
	fname += std::to_string(static_cast<long long>(numVertices));
	fname += "___e_";
	fname += std::to_string(static_cast<long long>(numEdges));
	fname += "___f_";
	fname += std::to_string(static_cast<long long>(numFaces));

	AdjustStorage();

	std::string maname = fname;
	maname += ".ma";

	std::ofstream fout(maname);

	//	GraphVertexIterator gvi,gvi_end;

	fout << numVertices << " " << numEdges << " " << numFaces << std::endl;

	//fout << num_vertices(*g) << " " << num_edges(*g) << " " << g->tris.size() << std::endl;

	for(unsigned i = 0; i < vertices.size(); i ++)
		//fout << "v " << vertices[i].second->sphere.center << " " << vertices[i].second->sphere.radius << std::endl;
		fout << "v " << setiosflags(ios::fixed) << setprecision(15) << (vertices[i].second->sphere.center * pmesh->bb_diagonal_length) << " " << (vertices[i].second->sphere.radius * pmesh->bb_diagonal_length) << std::endl;

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
	fout.close();
}


void SlabMesh::Export_OBJ(std::string fname) {
	fname += "___v_";
	fname += std::to_string(static_cast<long long>(numVertices));
	fname += "___e_";
	fname += std::to_string(static_cast<long long>(numEdges));
	fname += "___f_";
	fname += std::to_string(static_cast<long long>(numFaces));

	AdjustStorage();

	std::string maname = fname;
	maname += ".obj";

	std::ofstream fout(maname);

	//	GraphVertexIterator gvi,gvi_end;


	//fout << num_vertices(*g) << " " << num_edges(*g) << " " << g->tris.size() << std::endl;

	for (unsigned i = 0; i < vertices.size(); i++)
		fout << "v " << setiosflags(ios::fixed) << setprecision(15) << (vertices[i].second->sphere.center * pmesh->bb_diagonal_length) << " " << (vertices[i].second->sphere.radius * pmesh->bb_diagonal_length) << std::endl;

	for (unsigned i = 0; i < edges.size(); i++)
		fout << "l " << edges[i].second->vertices_.first+1 << " " << edges[i].second->vertices_.second+1 << std::endl;
	for (unsigned i = 0; i < faces.size(); i++)
	{
		fout << "f";
		for (std::set<unsigned>::iterator si = faces[i].second->vertices_.begin();
			si != faces[i].second->vertices_.end(); si++)
			fout << " " << *si+1;
		fout << std::endl;
	}
	fout.close();
}




void SlabMesh::InitialTopologyProperty(unsigned vid) {
	if (numVertices > 50)
		return;

	set<unsigned> fir_faces = vertices[vid].second->faces_;
	set<unsigned> fir_edges = vertices[vid].second->edges_;
	for (set<unsigned>::iterator si = fir_edges.begin(); si != fir_edges.end(); si++)
	{
		int index = edges[*si].second->vertices_.first == vid ? 
			edges[*si].second->vertices_.second : edges[*si].second->vertices_.first;

		set<unsigned> sec_edges = vertices[index].second->edges_;
		sec_edges.erase(*si);
		for (set<unsigned>::iterator si2 = sec_edges.begin(); si2 != sec_edges.end(); si2++)
		{
			int index2 = edges[*si2].second->vertices_.first == index ? 
				edges[*si2].second->vertices_.second : edges[*si2].second->vertices_.first;

			set<unsigned> third_edges = vertices[index2].second->edges_;
			third_edges.erase(*si2);
			for (set<unsigned>::iterator si3 = third_edges.begin(); si3 != third_edges.end(); si3++)
			{
				int index3 = edges[*si3].second->vertices_.first == index2 ? 
					edges[*si3].second->vertices_.second : edges[*si3].second->vertices_.first;

				if (index3 == vid)
				{
					// �������γ��˻�·��������滹��hole
					bool is_hole = true;
					for (set<unsigned>::iterator fi = fir_faces.begin(); fi != fir_faces.end(); fi++) 
					{
						set<unsigned> ver = faces[*fi].second->vertices_;
						if (ver.find(index) != ver.end() && ver.find(index2) != ver.end())
						{
							is_hole = false;
							break;
						}
					}

					if (is_hole)
					{
						if (edges[*si].second->faces_.size() <= 1 && edges[*si2].second->faces_.size() <= 1 
							&& edges[*si3].second->faces_.size() <= 1)
						{
							edges[*si].second->topo_contractable = false;
							edges[*si2].second->topo_contractable = false;
							edges[*si3].second->topo_contractable = false; 
						}
					}
				}
			}
		}
	}

	// ���򻯵�С��50������ʱ������������˵�ı߽��кϲ�
	if (numVertices <= 50)
	{
		for (set<unsigned>::iterator si = fir_edges.begin(); si != fir_edges.end(); si++)
		{
			int index = edges[*si].second->vertices_.first == vid ? 
				edges[*si].second->vertices_.second : edges[*si].second->vertices_.first;

			if (vertices[index].second->edges_.size() == 1 && vertices[index].second->faces_.size() == 0)
			{
				edges[*si].second->topo_contractable = false;
			}
		}
	}
}

void SlabMesh::InitialTopologyProperty() {
	for (int i = 0 ;i < vertices.size(); i++)
	{
		if(vertices[i].first)
		{
			InitialTopologyProperty(i);
		}
	}
}