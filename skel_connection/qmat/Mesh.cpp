#include "Mesh.h"

#include <CGAL/centroid.h>

#include <Eigen/Dense>

double Triangulation::TetCircumRadius(const Tetrahedron & tet)
{
	return (to_wm4(tet.vertex(0))-to_wm4(CGAL::circumcenter(tet))).Length();
}

double Triangulation::TetLargestEdgeLength(const Tetrahedron & tet)
{
	Vector3d tv0(to_wm4(tet.vertex(0)));
	Vector3d tv1(to_wm4(tet.vertex(1)));
	Vector3d tv2(to_wm4(tet.vertex(2)));
	Vector3d tv3(to_wm4(tet.vertex(3)));

	double length01,length02,length03,length12,length13,length23;
	//double maxlength;
	length01 = (tv0-tv1).Length();
	length02 = (tv0-tv2).Length();
	length03 = (tv0-tv3).Length();

	length12 = (tv1-tv2).Length();
	length13 = (tv1-tv3).Length();

	length23 = (tv2-tv3).Length();
	return std::max(std::max(std::max(length01,length02),std::max(length03,length12)),std::max(length13,length23));
}

double Triangulation::TetBoundingSphereRadius(const Tetrahedron & tet)
{
	double tv0coor[3];
	double tv1coor[3];
	double tv2coor[3];
	double tv3coor[3];

	tv0coor[0] = tet.vertex(0)[0]; tv0coor[1] = tet.vertex(0)[1]; tv0coor[2] = tet.vertex(0)[2];
	tv1coor[0] = tet.vertex(1)[0]; tv1coor[1] = tet.vertex(1)[1]; tv1coor[2] = tet.vertex(1)[2];
	tv2coor[0] = tet.vertex(2)[0]; tv2coor[1] = tet.vertex(2)[1]; tv2coor[2] = tet.vertex(2)[2];
	tv3coor[0] = tet.vertex(3)[0]; tv3coor[1] = tet.vertex(3)[1]; tv3coor[2] = tet.vertex(3)[2];

	MSPoint pts[4];
	pts[0] = MSPoint(3,tv0coor,tv0coor+3);
	pts[1] = MSPoint(3,tv1coor,tv1coor+3);
	pts[2] = MSPoint(3,tv2coor,tv2coor+3);
	pts[3] = MSPoint(3,tv3coor,tv3coor+3);

	Min_sphere minsphere(pts, pts+3);
		
	return std::sqrt(minsphere.squared_radius());
}


MPMesh::MPMesh()
{
	domain = NULL;
	m_max_v_gaussiancurvature = 1.;
	m_min_v_gaussiancurvature = 0.;
	m_max_v_meancurvature = 1.;
	m_min_v_meancurvature = 0.;
	m_max_v_abs_gaussiancurvature = 1.;
	m_min_v_abs_gaussiancurvature = 0.;
	m_max_v_ratio_principalcurvatures = 1.;
	m_min_v_ratio_principalcurvatures = 0.;
	m_max_v_maxcurvature = 1.;
	m_min_v_maxcurvature = 0.;
	m_max_v_mincurvature = 1.;
	m_min_v_mincurvature = 0.;

	m_max_f_gaussiancurvature = 1.;
	m_min_f_gaussiancurvature = 0.;
	m_max_f_meancurvature = 1.;
	m_min_f_meancurvature = 0.;
	m_max_f_abs_gaussiancurvature = 1.;
	m_min_f_abs_gaussiancurvature = 0.;
	m_max_f_ratio_principalcurvatures = 1.;
	m_min_f_ratio_principalcurvatures = 0.;
	m_max_f_maxcurvature = 1.;
	m_min_f_maxcurvature = 0.;
	m_max_f_mincurvature = 1.;
	m_min_f_mincurvature = 0.;

	m_max_gaussiancurvature_Meyer = 1e20;
	m_max_gaussiancurvature_Meyer = -1e20;

	m_density_policy = SQRTK1K2;
	m_metric_policy = K1_K2;

	number_of_bad_vertices = 0;
}

// compute the bounding box
void MPMesh::computebb()
{
	m_min[0] = 1e20;
	m_min[1] = 1e20;
	m_min[2] = 1e20;
	m_max[0] = -1e20;
	m_max[1] = -1e20;
	m_max[2] = -1e20;

	MPMesh::Vertex_iterator pVertex = vertices_begin();
	for(; pVertex != vertices_end(); pVertex ++)
	{
		if(pVertex->point()[0] < m_min[0])
			m_min[0] = pVertex->point()[0];
		if(pVertex->point()[1] < m_min[1])
			m_min[1] = pVertex->point()[1];
		if(pVertex->point()[2] < m_min[2])
			m_min[2] = pVertex->point()[2];

		if(pVertex->point()[0] > m_max[0])
			m_max[0] = pVertex->point()[0];
		if(pVertex->point()[1] > m_max[1])
			m_max[1] = pVertex->point()[1];
		if(pVertex->point()[2] > m_max[2])
			m_max[2] = pVertex->point()[2];	
	}

	bb_diagonal_length = sqrt((m_max[0] - m_min[0]) * (m_max[0] - m_min[0]) + (m_max[1] - m_min[1]) * (m_max[1] - m_min[1])
						+(m_max[2] - m_min[2]) * (m_max[2] - m_min[2]));
	bb_diagonal_length = 1; // this is modified by Zhiyang. for error checking...

}

void MPMesh::GenerateVertexList()
{
	pVertexList.clear();
	Vertex_iterator pVertex;
	int idx = 0;
	pVertex = vertices_begin();
	for(; pVertex != vertices_end(); pVertex ++, idx ++)
	{
		pVertexList.push_back(pVertex);
		pVertex->id = idx;
	}
	/*
	double scale_factor = 1.0 / bb_diagonal_length; 
	// 将顶点归一化
	for(unsigned i = 0; i < pVertexList.size(); i ++)
	{
		pVertexList[i]->point()[0] *= scale_factor;
		pVertexList[i]->point()[1] *= scale_factor;
		pVertexList[i]->point()[2] *= scale_factor;
	}*/
}

void MPMesh::GenerateFaceList()
{
	pFaceList.clear();
	Facet_iterator pFacet;
	int idx = 0;
	pFacet = facets_begin();
	for(; pFacet != facets_end(); pFacet ++, idx ++)
	{
		pFaceList.push_back(pFacet);
		pFacet->id = idx;
	}
}

void MPMesh::GenerateList()
{
	GenerateVertexList();
	GenerateFaceList();
}

void MPMesh::GenerateRandomColor()
{
	srand( (unsigned int)(time(NULL)) );
	for(unsigned int i = 0; i < pVertexList.size(); i ++)
	{
		pVertexList[i]->color[0] = rand() / (double) RAND_MAX;
		pVertexList[i]->color[1] = rand() / (double) RAND_MAX;
		pVertexList[i]->color[2] = rand() / (double) RAND_MAX;
	}
	
	for(unsigned int i = 0; i < pFaceList.size(); i ++)
	{
		pFaceList[i]->color[0] = rand() / (double) RAND_MAX;
		pFaceList[i]->color[1] = rand() / (double) RAND_MAX;
		pFaceList[i]->color[2] = rand() / (double) RAND_MAX;
	}
}

// compute the matrix of A and b for sphere mesh
void MPMesh::compute_sphere_matrix()
{
	for (int i = 0; i < pFaceList.size(); i ++)
	{
		MPMesh::Point p;
		compute_facet_center(pFaceList[i], p);
		//Wm4::Vector3d normal(pFaceList[i]->normal.x() / sqrt(pFaceList[i]->normal.squared_length()),
		//	pFaceList[i]->normal.y() / sqrt(pFaceList[i]->normal.squared_length()), 
		//	pFaceList[i]->normal.z() / sqrt(pFaceList[i]->normal.squared_length()));
		Wm4::Vector3d normal(pFaceList[i]->normal.x(), pFaceList[i]->normal.y(), pFaceList[i]->normal.z());
		Wm4::Vector3d point(p.x(), p.y(), p.z());

		// compute the matrix of A
		Wm4::Matrix3d normal_mul_normal;
		normal_mul_normal = normal_mul_normal.MakeTensorProduct(normal, normal);
		double A00 = normal_mul_normal.GetRow(0).X();
		double A01 = normal_mul_normal.GetRow(0).Y();
		double A02 = normal_mul_normal.GetRow(0).Z();
		double A03 = normal.X();
		double A10 = normal_mul_normal.GetRow(1).X();
		double A11 = normal_mul_normal.GetRow(1).Y();
		double A12 = normal_mul_normal.GetRow(1).Z();
		double A13 = normal.Y();
		double A20 = normal_mul_normal.GetRow(2).X();
		double A21 = normal_mul_normal.GetRow(2).Y();
		double A22 = normal_mul_normal.GetRow(2).Z();
		double A23 = normal.Z();
		double A30 = normal.X();
		double A31 = normal.Y();
		double A32 = normal.Z();
		double A33 = 1.0;
		Wm4::Matrix4d temp_A(2 * A00, 2 * A01, 2 * A02, 2 * A03, 
			2 * A10, 2 * A11, 2 * A12, 2 * A13, 
			2 * A20, 2 * A21, 2 * A22, 2 * A23, 
			2 * A30, 2 * A31, 2 * A32, 2 * A33);

		//Wm4::Vector4d temp_normal1(pFaceList[i]->normal.x(), pFaceList[i]->normal.y(), pFaceList[i]->normal.z(), 1.0);
		//Matrix4d temp_A1, temp_A2;
		//temp_A1.MakeTensorProduct(temp_normal1, temp_normal1);
		//temp_A1 *= 2.0;

		// compute the matrix of b
		double normal_mul_point = normal.Dot(point);
		Wm4::Vector4d temp_b(2 * normal_mul_point * normal.X(), 2 * normal_mul_point * normal.Y(), 
			2 * normal_mul_point * normal.Z(), 2 * normal_mul_point);

		//compute c
		double temp_c = normal_mul_point * normal_mul_point;

		Halfedge_around_facet_circulator pHalfedge = pFaceList[i]->facet_begin();
		Halfedge_around_facet_circulator end = pHalfedge;
		// compute the area of the triangle
		Point triangle_vertex[3];
		int index = 0;
		CGAL_For_all(pHalfedge, end)
		{
			triangle_vertex[index] = pHalfedge->vertex()->point();
			index++;
		}
		Wm4::Vector3d ab(triangle_vertex[0].x() - triangle_vertex[1].x(), 
			triangle_vertex[0].y() - triangle_vertex[1].y(),
			triangle_vertex[0].z() - triangle_vertex[1].z());
		Wm4::Vector3d ac(triangle_vertex[0].x() - triangle_vertex[2].x(), 
			triangle_vertex[0].y() - triangle_vertex[2].y(),
			triangle_vertex[0].z() - triangle_vertex[2].z());
		double area = ab.Cross(ac).Length() / 2.0;

		// traverse the vertex in the face
		pHalfedge = pFaceList[i]->facet_begin();
		CGAL_For_all(pHalfedge, end)
		{
			//pHalfedge->vertex()->A += temp_A * area / 3.0;
			//pHalfedge->vertex()->b += temp_b * area / 3.0;
			//pHalfedge->vertex()->c += temp_c * area / 3.0;
			pHalfedge->vertex()->A += temp_A;
			pHalfedge->vertex()->b += temp_b;
			pHalfedge->vertex()->c += temp_c;
		}
	}
}

void MPMesh::compute_normals_per_facet()
{
	std::for_each( facets_begin(), facets_end(), Facet_normal() );
}

void MPMesh::compute_normals_per_vertex()
{
	std::for_each( vertices_begin(), vertices_end(), Vertex_normal() );
}

void MPMesh::compute_normals()
{
	compute_normals_per_facet();
	compute_normals_per_vertex();
}

void MPMesh::copybb(MPMesh * pmesh)
{
	m_min[0] = pmesh->m_min[0];
	m_min[1] = pmesh->m_min[1];
	m_min[2] = pmesh->m_min[2];
	
	m_max[0] = pmesh->m_max[0];
	m_max[1] = pmesh->m_max[1];
	m_max[2] = pmesh->m_max[2];
} 

void MPMesh::tag_facets(const int tag)
{
	for( Facet_iterator pFace = facets_begin(); pFace != facets_end(); pFace ++)
		pFace->tag = tag;
}

void MPMesh::tag_halfedges(const int tag)
{
	for( Halfedge_iterator pHalfedge = halfedges_begin(); pHalfedge != halfedges_end(); pHalfedge ++)
		pHalfedge->tag = tag;
}

void MPMesh::tag_vertices(const int tag)
{
	for( Vertex_iterator pVertex = vertices_begin(); pVertex != vertices_end(); pVertex ++)
		pVertex->tag = tag;
}

void MPMesh::compute_facet_center(Facet_handle pFace, Point & center)
{
	Halfedge_around_facet_circulator pHalfedge = pFace->facet_begin();
	Halfedge_around_facet_circulator end = pHalfedge;
	Vector vec(0.0, 0.0, 0.0);
	int degree = 0;
	CGAL_For_all(pHalfedge, end)
	{
		vec = vec + (pHalfedge->vertex()->point() - CGAL::ORIGIN);
		degree ++;
	}
	center = CGAL::ORIGIN + (vec / (simple_kernel::FT) degree);
}

Halfedge_handle MPMesh::get_border_halfedge_tag(int tag)
{
	for(Halfedge_iterator pHalfedge = halfedges_begin(); pHalfedge != halfedges_end(); pHalfedge ++)
		if( (pHalfedge->is_border()) && (pHalfedge->tag == tag) )
			return pHalfedge;
	return NULL;
}

Facet_handle MPMesh::get_facet_tag(const int tag)
{
	for(Facet_iterator pFace = facets_begin(); pFace != facets_end(); pFace ++)
		if(pFace->tag == tag)
			return pFace;
	return NULL;
}

void MPMesh::tag_component(Facet_handle pSeedFacet, const int tag_free, const int tag_done)
{
	pSeedFacet->tag = tag_done;
	std::list<Facet_handle> facets;

	facets.push_front(pSeedFacet);
	while( !facets.empty() )
	{
		Facet_handle pFacet = facets.front();
		facets.pop_front();
		pFacet->tag = tag_done;
		Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
		Halfedge_around_facet_circulator end = pHalfedge;
		CGAL_For_all(pHalfedge, end)
		{
			Facet_handle pNFacet = pHalfedge->opposite()->facet();
			if( (pNFacet != NULL) && (pNFacet->tag == tag_free) )
			{
				facets.push_front(pNFacet);
				pNFacet->tag = tag_done;
			}
		}
	}
}

unsigned int MPMesh::nb_boundaries()
{
	unsigned int nb = 0;
	tag_halfedges(0);
	Halfedge_handle seed_halfedge = NULL;
	while( ( seed_halfedge = get_border_halfedge_tag(0) ) != NULL )
	{
		nb ++;
		seed_halfedge->tag = 1;
		Vertex_handle seed_vertex = seed_halfedge->prev()->vertex();
		Halfedge_handle current_halfedge = seed_halfedge;
		Halfedge_handle next_halfedge;
		do
		{
			next_halfedge = current_halfedge->next();
			next_halfedge->tag = 1;
			current_halfedge = next_halfedge;
		} while ( next_halfedge->prev()->vertex() != seed_vertex );
	}
	return nb;
}

unsigned int MPMesh::nb_components()
{
	unsigned int nb = 0;
	tag_facets(0);
	Facet_handle seed_facet = NULL;
	while( ( seed_facet = get_facet_tag(0) ) != NULL )
	{
		nb ++;
		tag_component(seed_facet, 0, 1);
	}
	return nb;
}

bool MPMesh::is_simple_watertight()
{
	if( (nb_components() == 1) && (nb_boundaries() == 0) )
		return true;
	return false;
}

int MPMesh::genus()
{
	int c = (int)(nb_components());
	int b = (int)(nb_boundaries());
	int v = (int)(size_of_vertices());
	int e = (int)(size_of_halfedges()) / 2;
	int f = (int)(size_of_facets());
	return genus(c,v,f,e,b);
}

int MPMesh::genus(int c, int v, int f, int e, int b)
{
	return (2*c+e-b-f-v)/2;
}

void MPMesh::compute_components_boundaries_genus()
{
	m_nb_components = nb_components();
	m_nb_boundaries = nb_boundaries();
	m_genus = genus();
}

unsigned int MPMesh::get_nb_components()
{
	return m_nb_components;
}

unsigned int MPMesh::get_nb_boundaries()
{
	return m_nb_boundaries;
}

int MPMesh::get_genus()
{
	return m_genus;
}

Wm4::Vector3d MPMesh::GetCentroid(Facet_iterator pFace)
{
	Halfedge_around_facet_circulator hafc = pFace->facet_begin();
	Wm4::Vector3d centroid = Wm4::Vector3d::ZERO;
	int degree = (int)(CGAL::circulator_size( pFace->facet_begin() ));
	do
	{
		centroid += to_wm4( hafc->vertex()->point() );
	} while( ++hafc != pFace->facet_begin() );
	centroid /= degree;
	return centroid;
}

double MPMesh::GetFaceLargestAngle(Facet_handle pFace)
{
	Halfedge_around_facet_circulator hafc;
	hafc = pFace->facet_begin();
	Wm4::Vector3d v0, v1, v2;
	v0 = to_wm4( hafc->vertex()->point() );
	v1 = to_wm4( hafc->next()->vertex()->point() );
	v2 = to_wm4( hafc->next()->next()->vertex()->point() );

	double a, b, c;
	a = (v0-v1).Length();
	b = (v1-v2).Length();
	c = (v2-v0).Length();
	
	if(a*b*c == 0)
		return 180.;
	
	double anglea, angleb, anglec;
	anglea = 180. * acos( (b*b+c*c-a*a) / (2.*b*c) ) / MY_PI;
	angleb = 180. * acos( (a*a+c*c-b*b) / (2.*a*c) ) / MY_PI;
	anglec = 180. * acos( (a*a+b*b-c*c) / (2.*a*b) ) / MY_PI;

	if( (anglea >= angleb) && (anglea >= anglec) )
		return anglea;
	if( (angleb >= anglea) && (angleb >= anglec) )
		return angleb;
	return anglec;
	
}

double MPMesh::GetFaceSmallestAngle(Facet_handle pFace)
{
	Halfedge_around_facet_circulator hafc;
	hafc = pFace->facet_begin();
	Wm4::Vector3d v0, v1, v2;
	v0 = to_wm4( hafc->vertex()->point() );
	v1 = to_wm4( hafc->next()->vertex()->point() );
	v2 = to_wm4( hafc->next()->next()->vertex()->point() );

	double a, b, c;
	a = (v0-v1).Length();
	b = (v1-v2).Length();
	c = (v2-v0).Length();
	
	if(a*b*c == 0)
		return 180.;
	
	double anglea, angleb, anglec;

	anglea = 180. * acos( (b*b+c*c-a*a) / (2.*b*c) ) / MY_PI;
	angleb = 180. * acos( (a*a+c*c-b*b) / (2.*a*c) ) / MY_PI;
	anglec = 180. * acos( (a*a+b*b-c*c) / (2.*a*b) ) / MY_PI;

	if( (anglea <= angleb) && (anglea <= anglec) )
		return anglea;
	if( (angleb <= anglea) && (angleb <= anglec) )
		return angleb;
	return anglec;
}

double MPMesh::GetArea(Facet_handle pFace)
{
	Wm4::Vector3d cent = GetCentroid(pFace);

	double a(0.0);
	Halfedge_around_facet_circulator h = pFace->facet_begin();
	do
	{
		a += 0.5 * fabs( ( (to_wm4(h->next()->vertex()->point())-cent).Cross(to_wm4(h->vertex()->point())-cent) ).Length() );
	} while(++h != pFace->facet_begin());
	return a;
}

bool MPMesh::inside(const Wm4::Vector3d & p)
{
	// fast culling
	if(!inside_boundingbox(p))
		return false;

	if(is_simple_watertight() == false)
		return false;

	int insideCount = 0;

	int NumberRays = 1;
	
	for(int i = 0; i < NumberRays; i ++)
	{
		Wm4::Vector3d dir = RandomPointonSphere();
		bool odd = false;

		for(int j = 0; j < pFaceList.size(); j ++)
		{
			if( FastNoIntersect(p, dir, to_wm4(pFaceList[j]->facet_begin()->vertex()->point()), to_wm4(pFaceList[j]->normal)) )
				continue;

			if(RayTriangleIntersect(p,dir,to_wm4(pFaceList[j]->facet_begin()->vertex()->point()),
				to_wm4(pFaceList[j]->facet_begin()->next()->vertex()->point()),
				to_wm4(pFaceList[j]->facet_begin()->next()->next()->vertex()->point())))
				odd = !odd;
		}

		if(odd)
			insideCount ++;
	}

	return insideCount > NumberRays/2;
}

bool MPMesh::inside_boundingbox(const Wm4::Vector3d & p)
{
	if( (p[0] < m_min[0]) || (p[0] > m_max[0]) || (p[1] < m_min[1]) || (p[1] > m_max[1]) || (p[2] < m_min[2]) || (p[2] > m_max[2]) )
		return false;
	return true;
}

void MPMesh::EvaluateCurvatureInterval()
{
	double gc, mc, absgc, pcr;
	m_max_v_gaussiancurvature = -1e20;
	m_min_v_gaussiancurvature = 1e20;
	m_max_v_meancurvature = -1e20;
	m_min_v_meancurvature = 1e20;
	m_max_v_abs_gaussiancurvature = -1e20;
	m_min_v_abs_gaussiancurvature = 1e20;
	m_max_v_ratio_principalcurvatures = -1e20;
	m_min_v_ratio_principalcurvatures = 1e20;
	m_max_v_maxcurvature = -1e20;
	m_min_v_maxcurvature = 1e20;
	m_max_v_mincurvature = -1e20;
	m_min_v_mincurvature = 1e20;
	for(unsigned int i = 0; i < pVertexList.size(); i ++)
	{
		gc = pVertexList[i]->maxcurvature * pVertexList[i]->mincurvature;
		mc = 0.5 * (pVertexList[i]->maxcurvature + pVertexList[i]->mincurvature);
		absgc = fabs(gc);
		pcr = min(fabs(pVertexList[i]->mincurvature),fabs(pVertexList[i]->maxcurvature)) / max(fabs(pVertexList[i]->mincurvature),fabs(pVertexList[i]->maxcurvature));
		/*
		if( fabs(pVertexList[i]->maxcurvature) > fabs(pVertexList[i]->mincurvature) )
			//pcr = 1.0 - fabs(pVertexList[i]->mincurvature)/fabs(pVertexList[i]->maxcurvature);
			pcr = fabs(pVertexList[i]->mincurvature)/fabs(pVertexList[i]->maxcurvature);
		else
			//pcr = 1.0 - fabs(pVertexList[i]->maxcurvature)/fabs(pVertexList[i]->mincurvature);
			pcr = fabs(pVertexList[i]->maxcurvature)/fabs(pVertexList[i]->mincurvature);
		*/

		if(gc > m_max_v_gaussiancurvature)
			m_max_v_gaussiancurvature = gc;
		if(gc < m_min_v_gaussiancurvature)
			m_min_v_gaussiancurvature = gc;
		if(mc > m_max_v_meancurvature)
			m_max_v_meancurvature = mc;
		if(mc < m_min_v_meancurvature)
			m_min_v_meancurvature = mc;
		if(absgc > m_max_v_abs_gaussiancurvature)
			m_max_v_abs_gaussiancurvature = absgc;
		if(absgc < m_min_v_abs_gaussiancurvature)
			m_min_v_abs_gaussiancurvature = absgc;
		if(pcr > m_max_v_ratio_principalcurvatures)
			m_max_v_ratio_principalcurvatures = pcr;
		if(pcr < m_min_v_ratio_principalcurvatures)
			m_min_v_ratio_principalcurvatures = pcr;

		if(pVertexList[i]->maxcurvature > m_max_v_maxcurvature)
			m_max_v_maxcurvature = pVertexList[i]->maxcurvature;
		if(pVertexList[i]->maxcurvature < m_min_v_maxcurvature)
			m_min_v_maxcurvature = pVertexList[i]->maxcurvature;

		if(pVertexList[i]->mincurvature > m_max_v_mincurvature)
			m_max_v_mincurvature = pVertexList[i]->mincurvature;
		if(pVertexList[i]->mincurvature < m_min_v_mincurvature)
			m_min_v_mincurvature = pVertexList[i]->mincurvature;
	}

	m_max_f_gaussiancurvature = -1e20;
	m_min_f_gaussiancurvature = 1e20;
	m_max_f_meancurvature = -1e20;
	m_min_f_meancurvature = 1e20;
	m_max_f_abs_gaussiancurvature = -1e20;
	m_min_f_abs_gaussiancurvature = 1e20;
	m_max_f_ratio_principalcurvatures = -1e20;
	m_min_f_ratio_principalcurvatures = 1e20;
	m_max_f_maxcurvature = -1e20;
	m_min_f_maxcurvature = 1e20;
	m_max_f_mincurvature = -1e20;
	m_min_f_mincurvature = 1e20;
	for(unsigned int i = 0; i < pFaceList.size(); i ++)
	{
		gc = pFaceList[i]->maxcurvature * pFaceList[i]->mincurvature;
		mc = 0.5 * (pFaceList[i]->maxcurvature + pFaceList[i]->mincurvature);
		absgc = fabs(gc);
		pcr = min(fabs(pFaceList[i]->maxcurvature),fabs(pFaceList[i]->mincurvature)) / max(fabs(pFaceList[i]->maxcurvature),fabs(pFaceList[i]->mincurvature));
		/*
		if( fabs(pFaceList[i]->maxcurvature) > fabs(pFaceList[i]->mincurvature) )
			//pcr = 1.0 - fabs(pFaceList[i]->mincurvature)/fabs(pFaceList[i]->maxcurvature);
			pcr = fabs(pFaceList[i]->mincurvature)/fabs(pFaceList[i]->maxcurvature);
		else
			//pcr = 1.0 - fabs(pFaceList[i]->maxcurvature)/fabs(pFaceList[i]->mincurvature);
			pcr = fabs(pFaceList[i]->maxcurvature)/fabs(pFaceList[i]->mincurvature);
		*/

		if(gc > m_max_f_gaussiancurvature)
			m_max_f_gaussiancurvature = gc;
		if(gc < m_min_f_gaussiancurvature)
			m_min_f_gaussiancurvature = gc;
		if(mc > m_max_f_meancurvature)
			m_max_f_meancurvature = mc;
		if(mc < m_min_f_meancurvature)
			m_min_f_meancurvature = mc;
		if(absgc > m_max_f_abs_gaussiancurvature)
			m_max_f_abs_gaussiancurvature = absgc;
		if(absgc < m_min_f_abs_gaussiancurvature)
			m_min_f_abs_gaussiancurvature = absgc;
		if(pcr > m_max_f_ratio_principalcurvatures)
			m_max_f_ratio_principalcurvatures = pcr;
		if(pcr < m_min_f_ratio_principalcurvatures)
			m_min_f_ratio_principalcurvatures = pcr;

		if(pFaceList[i]->maxcurvature > m_max_f_maxcurvature)
			m_max_f_maxcurvature = pFaceList[i]->maxcurvature;
		if(pFaceList[i]->maxcurvature < m_min_f_maxcurvature)
			m_min_f_maxcurvature = pFaceList[i]->maxcurvature;

		if(pFaceList[i]->mincurvature > m_max_f_mincurvature)
			m_max_f_mincurvature = pFaceList[i]->mincurvature;
		if(pFaceList[i]->mincurvature < m_min_f_mincurvature)
			m_min_f_mincurvature = pFaceList[i]->mincurvature;
	}
}

void MPMesh::ComputeFaceCurvatures(double MinThreshold, double RatioThreshold)
{
	std::cout << "Computing Face Curvatures..." << std::endl;
	for(unsigned int i = 0; i < pFaceList.size(); i ++)
	{
		Halfedge_around_facet_circulator pHalfedge = pFaceList[i]->facet_begin();

		pFaceList[i]->normdir = to_wm4(pFaceList[i]->normal);
		Eigen::Matrix3d maxmat(Eigen::Matrix3d::Zero()), minmat(Eigen::Matrix3d::Zero());
		double maxcur(0), mincur(0);		
		int count(0);
		do
		{
			Eigen::MatrixXd maxdir(3,1), mindir(3,1);
			maxdir << pHalfedge->vertex()->maxdir[0], pHalfedge->vertex()->maxdir[1], pHalfedge->vertex()->maxdir[2];
			mindir << pHalfedge->vertex()->mindir[0], pHalfedge->vertex()->mindir[1], pHalfedge->vertex()->mindir[2];
			maxmat += maxdir*maxdir.transpose();
			minmat += mindir*mindir.transpose();
			maxcur += pHalfedge->vertex()->maxcurvature;
			mincur += pHalfedge->vertex()->mincurvature;
			count ++;
		} while(++pHalfedge != pFaceList[i]->facet_begin());
		maxcur /= count;
		mincur /= count;

		Eigen::EigenSolver<Eigen::MatrixXd> esmax(maxmat), esmin(minmat);
		Eigen::VectorXcd esmaxval(esmax.eigenvalues()), esminval(esmin.eigenvalues());
		int maxindex, minindex;
		maxindex = MaxIndexOf3( esmaxval[0].real(), esmaxval[1].real(), esmaxval[2].real() );
		minindex = MaxIndexOf3( esminval[0].real(), esminval[1].real(), esminval[2].real() );
		
		Vector3d fmaxdir(esmax.eigenvectors().col(maxindex)[0].real(),esmax.eigenvectors().col(maxindex)[1].real(),esmax.eigenvectors().col(maxindex)[2].real()), fmindir(esmin.eigenvectors().col(minindex)[0].real(),esmin.eigenvectors().col(minindex)[1].real(),esmin.eigenvectors().col(minindex)[2].real());
		
		fmaxdir = fmaxdir - fmaxdir.Dot(pFaceList[i]->normdir)*pFaceList[i]->normdir;
		fmaxdir.Normalize();
		fmindir = fmindir - fmindir.Dot(pFaceList[i]->normdir)*pFaceList[i]->normdir;
		fmindir = fmindir - fmindir.Dot(fmaxdir)*fmaxdir;
		fmindir.Normalize();

		pFaceList[i]->maxdir = fmaxdir;
		pFaceList[i]->mindir = fmindir;
		pFaceList[i]->maxcurvature = maxcur;
		pFaceList[i]->mincurvature = mincur;

		
		if(fabs(pFaceList[i]->maxcurvature) < MinThreshold)
		{
			if(pFaceList[i]->maxcurvature > 0)
				pFaceList[i]->maxcurvature = MinThreshold;
			else
				pFaceList[i]->maxcurvature = -MinThreshold;
		}
		
		if(fabs(pFaceList[i]->mincurvature) < MinThreshold)
		{
			if(pFaceList[i]->mincurvature > 0)
				pFaceList[i]->mincurvature = MinThreshold;
			else
				pFaceList[i]->mincurvature = -MinThreshold;
		}

		double absmax, absmin;
		absmax = fabs(pFaceList[i]->maxcurvature);
		absmin = fabs(pFaceList[i]->mincurvature);
		double ratio = absmax / absmin;

		if(ratio < 1)
			ratio = 1. / ratio;

		if(ratio > RatioThreshold)
		{
			if(absmax > absmin)
				absmax = RatioThreshold * absmin;
			else
				absmin = RatioThreshold * absmax;
		}

		if(pFaceList[i]->maxcurvature > 0)
			pFaceList[i]->maxcurvature = absmax;
		else
			pFaceList[i]->maxcurvature = -absmax;

		if(pFaceList[i]->mincurvature > 0)
			pFaceList[i]->mincurvature = absmin;
		else
			pFaceList[i]->mincurvature = -absmin;
	}
	ComputeFaceDensity();
	std::cout << "Done." << std::endl;
}

void MPMesh::ComputeFaceDensity()
{
	std::cout << "Computing Face Density..." << std::endl;
	for(unsigned int i = 0; i < pFaceList.size(); i ++)
	{
		if(m_density_policy == ONE)
			pFaceList[i]->density = 1.;
		else if(m_density_policy == SQRTK1K2)
			pFaceList[i]->density = sqrt(fabs(pFaceList[i]->maxcurvature*pFaceList[i]->mincurvature));
		else if(m_density_policy == K1K2)
			pFaceList[i]->density = fabs(pFaceList[i]->maxcurvature*pFaceList[i]->mincurvature);
	}
	std::cout << "Done." << std::endl;
}

void MPMesh::ComputeVertexDensity()
{
	std::cout << "Computing Vertex Density..." << std::endl;
	Vertex_iterator pVertex;
	for(pVertex = vertices_begin(); pVertex != vertices_end(); pVertex ++)
	{
		if(m_density_policy == ONE)
			pVertex->density = 1.;
		else if(m_density_policy == SQRTK1K2)
			pVertex->density = sqrt(fabs(pVertex->maxcurvature*pVertex->mincurvature));
		else if(m_density_policy == K1K2)
			pVertex->density = fabs(pVertex->maxcurvature*pVertex->mincurvature);
	}
	std::cout << "Done." << std::endl;
}

void MPMesh::BuildVertexMetricTensors(double normal_coefficient)
{
	std::cout << "Building Vertex Metric Tensors..." << std::endl;
	for(Vertex_iterator pVertex = vertices_begin(); pVertex != vertices_end(); pVertex ++)
	{
		// build the metric tensor
		Matrix3d pf(pVertex->maxdir, pVertex->mindir, pVertex->normdir, true);
		Matrix3d pt;
		if(m_metric_policy == K1_K2)
			pt = Matrix3d(fabs(pVertex->maxcurvature), fabs(pVertex->mincurvature), normal_coefficient);
		else if(m_metric_policy == SQUAREK1_SQUAREK2)
			pt = Matrix3d(pow(pVertex->maxcurvature,2), pow(pVertex->mincurvature,2), normal_coefficient);
		//std::cout << pt[0][0] << " " << pt[1][1] << " " << pt[2][2] << std::endl << std::endl;
		pVertex->metric = pf * pt * (pf.Transpose());
	}
	std::cout << "Done." << std::endl;
}

void MPMesh::BuildFaceMetricTensors(double normal_coefficient)
{
	std::cout << "Building Face Metric Tensors..." << std::endl;
	for(unsigned int i = 0; i < pFaceList.size(); i ++)
	{
		// build the metric tensor
		Matrix3d pf(pFaceList[i]->maxdir, pFaceList[i]->mindir, pFaceList[i]->normdir, true);
		Matrix3d pt;
		if(m_metric_policy == K1_K2)
			pt = Matrix3d(fabs(pFaceList[i]->maxcurvature), fabs(pFaceList[i]->mincurvature), normal_coefficient);
		else if(m_metric_policy == SQUAREK1_SQUAREK2)
			pt = Matrix3d(pow(pFaceList[i]->maxcurvature,2), pow(pFaceList[i]->mincurvature,2), normal_coefficient);
		//std::cout << pt[0][0] << " " << pt[1][1] << " " << pt[2][2] << std::endl << std::endl;
		pFaceList[i]->metric = pf * pt * (pf.Transpose());
	}
	std::cout << "Done." << std::endl;
}

void MPMesh::LoadPrincipalCurvatures(string filename, double MinThreshold, double RatioThreshold, double IsotropicCoefficient)
{
	std::cout << "Loading Principal Curvatures..." << std::endl;
	ifstream fin(filename.c_str());
	string str;
	fin >> str >> str >> str >> str >> str;

	Vertex_iterator pVertex;
	int pVertexid;
	for(pVertex = vertices_begin(); pVertex != vertices_end(); pVertex ++)
	{
		// loading principal curvatures
		fin >> pVertexid;
		fin >> pVertex->maxcurvature;
		fin >> pVertex->maxdir[0] >> pVertex->maxdir[1] >> pVertex->maxdir[2];
		fin >> pVertex->mincurvature;
		fin >> pVertex->mindir[0] >> pVertex->mindir[1] >> pVertex->mindir[2];
		fin >> pVertex->normdir[0] >> pVertex->normdir[1] >> pVertex->normdir[2];

		if(fabs(pVertex->maxcurvature) < MinThreshold)
		{
			if(pVertex->maxcurvature > 0)
				pVertex->maxcurvature = MinThreshold;
			else
				pVertex->maxcurvature = -MinThreshold;
		}
		
		if(fabs(pVertex->mincurvature) < MinThreshold)
		{
			if(pVertex->mincurvature > 0)
				pVertex->mincurvature = MinThreshold;
			else
				pVertex->mincurvature = -MinThreshold;
		}

		double absmax, absmin;
		absmax = fabs(pVertex->maxcurvature);
		absmin = fabs(pVertex->mincurvature);
		double ratio = absmax / absmin;

		if(ratio < 1)
			ratio = 1. / ratio;

		if(ratio > RatioThreshold)
		{
			if(absmax > absmin)
				absmax = RatioThreshold * absmin;
			else
				absmin = RatioThreshold * absmax;
		}

		if(pVertex->maxcurvature > 0)
			pVertex->maxcurvature = absmax;
		else
			pVertex->maxcurvature = -absmax;

		if(pVertex->mincurvature > 0)
			pVertex->mincurvature = absmin;
		else
			pVertex->mincurvature = -absmin;


		if(pVertex->maxcurvature > 0)
			pVertex->maxcurvature += IsotropicCoefficient;
		else
			pVertex->maxcurvature -= IsotropicCoefficient;

		if(pVertex->mincurvature > 0)
			pVertex->mincurvature += IsotropicCoefficient;
		else
			pVertex->mincurvature -= IsotropicCoefficient;
	}
	fin.close();
	ComputeVertexDensity();
	std::cout << "Done." << std::endl;
}

double MPMesh::vertex_voronoi_area_Meyer(Vertex_handle vh)
{
	double area(0.0);
	
	Halfedge_around_vertex_circulator havc = vh->vertex_begin();
	Vector3d v = to_wm4(vh->point());
	do
	{
		if(havc->facet() != NULL)
		{
			Vector3d vp = to_wm4(havc->opposite()->vertex()->point());
			Vector3d va = to_wm4(havc->next()->vertex()->point());
		
			Triangle t(v,vp,va);
			area += t.voronoi_area_Meyer(0);
		}
		havc ++;
	}
	while( havc != vh->vertex_begin() );

	return area;
}

double MPMesh::vertex_voronoi_area_new(Vertex_handle vh)
{
	double area(0.0);
	
	Halfedge_around_vertex_circulator havc = vh->vertex_begin();
	Vector3d v = to_wm4(vh->point());
	do
	{
		if(havc->facet() != NULL)
		{
			Vector3d vp = to_wm4(havc->opposite()->vertex()->point());
			Vector3d va = to_wm4(havc->next()->vertex()->point());
		
			Triangle t(v,vp,va);
			area += t.voronoi_area_new(0);
		}
		havc ++;
	}
	while( havc != vh->vertex_begin() );

	return area;
}

double MPMesh::vertex_angle(Vertex_handle vh)
{
	double angle(0.0);
	
	Halfedge_around_vertex_circulator havc = vh->vertex_begin();
	Vector3d v = to_wm4(vh->point());
	do
	{
		if(havc->facet() != NULL)
		{
			Vector3d vp = to_wm4(havc->opposite()->vertex()->point());
			Vector3d va = to_wm4(havc->next()->vertex()->point());
		
			angle += angle_from_cotan(v,vp,va);
		}
		havc ++;
	}
	while( havc != vh->vertex_begin() );

	return angle;
}

bool MPMesh::vertex_gaussian_curvature_Meyer(Vertex_handle vh)
{
	double area = vertex_voronoi_area_Meyer(vh);
	
	double angledefect = 2. * Wm4::Mathd::PI - vertex_angle(vh);
	
	vh->gaussiancurvature_Meyer = angledefect / area;
	
	return true;
}

void MPMesh::ComputeVertexGaussianCurvature_Meyer()
{
	for(unsigned int i = 0; i < pVertexList.size(); i ++)
	{
		if(is_border(pVertexList[i]) == false)
			vertex_gaussian_curvature_Meyer(pVertexList[i]);
	}

	m_max_gaussiancurvature_Meyer = -1e20;
	m_min_gaussiancurvature_Meyer = 1e20;
	for(unsigned int i = 0; i < pVertexList.size(); i ++)
	{
		if(is_border(pVertexList[i]) == false)
		{
			if(pVertexList[i]->gaussiancurvature_Meyer > m_max_gaussiancurvature_Meyer)
				m_max_gaussiancurvature_Meyer = pVertexList[i]->gaussiancurvature_Meyer;
			if(pVertexList[i]->gaussiancurvature_Meyer < m_min_gaussiancurvature_Meyer)
				m_min_gaussiancurvature_Meyer = pVertexList[i]->gaussiancurvature_Meyer;
		}
	}
}

bool MPMesh::vertex_gaussian_curvature_new(Vertex_handle vh)
{
	double area = vertex_voronoi_area_new(vh);
	
	double angledefect = 2. * Wm4::Mathd::PI - vertex_angle(vh);
	
	vh->gaussiancurvature_new = angledefect / area;
	
	return true;
}

void MPMesh::ComputeVertexGaussianCurvature_new()
{
	for(unsigned int i = 0; i < pVertexList.size(); i ++)
	{
		if(is_border(pVertexList[i]) == false)
			vertex_gaussian_curvature_new(pVertexList[i]);
	}

	m_max_gaussiancurvature_new = -1e20;
	m_min_gaussiancurvature_new = 1e20;
	for(unsigned int i = 0; i < pVertexList.size(); i ++)
	{
		if(is_border(pVertexList[i]) == false)
		{
			if(pVertexList[i]->gaussiancurvature_new > m_max_gaussiancurvature_new)
				m_max_gaussiancurvature_new = pVertexList[i]->gaussiancurvature_Meyer;
			if(pVertexList[i]->gaussiancurvature_new < m_min_gaussiancurvature_new)
				m_min_gaussiancurvature_new = pVertexList[i]->gaussiancurvature_new;
		}
	}
}

Wm4::Vector3d MPMesh::vertex_cotan(Vertex_handle vh)
{
	Wm4::Vector3d sum(Wm4::Vector3d::ZERO);
	
	Halfedge_around_vertex_circulator havc = vh->vertex_begin();
	Vector3d v = to_wm4(vh->point());
	do
	{
		if(havc->facet() != NULL)
		{
			Vector3d vp = to_wm4(havc->opposite()->vertex()->point());
			Vector3d va = to_wm4(havc->next()->vertex()->point());
		
			Wm4::Vector3d vpv = v - vp;
			Wm4::Vector3d vav = v - va;
			
			double cotvp, cotva;
			cotvp = cotan(vp,v,va);
			cotva = cotan(va,v,vp);			

			sum += cotvp * vpv;
			sum += cotva * vav;
		}
		havc ++;
	}
	while( havc != vh->vertex_begin() );

	return sum;
}

bool MPMesh::vertex_mean_curvature_normal_Meyer(Vertex_handle vh)
{
	double area = vertex_voronoi_area_Meyer(vh);

	Wm4::Vector3d sum_cotan = vertex_cotan(vh);

	sum_cotan /= (2.*area);

	vh->normal_Meyer = sum_cotan;
	vh->normal_Meyer.Normalize();

	vh->meancurvature_Meyer =  sum_cotan.Length();
	
	return true;
}

void MPMesh::ComputeVertexMeanCurvatureNormal_Meyer()
{
	for(unsigned int i = 0; i < pVertexList.size(); i ++)
	{
		if(is_border(pVertexList[i]) == false)
			vertex_mean_curvature_normal_Meyer(pVertexList[i]);
	}

	m_max_meancurvature_Meyer = -1e20;
	m_min_meancurvature_Meyer = 1e20;
	for(unsigned int i = 0; i < pVertexList.size(); i ++)
	{
		if(is_border(pVertexList[i]) == false)
		{
			if(pVertexList[i]->meancurvature_Meyer > m_max_meancurvature_Meyer)
				m_max_meancurvature_Meyer = pVertexList[i]->meancurvature_Meyer;
			if(pVertexList[i]->meancurvature_Meyer < m_min_meancurvature_Meyer)
				m_min_meancurvature_Meyer = pVertexList[i]->meancurvature_Meyer;
		}
	}
}

bool MPMesh::vertex_principal_curvature_Meyer(Vertex_handle vh)
{
	double temp = vh->meancurvature_Meyer * vh->meancurvature_Meyer - vh->gaussiancurvature_Meyer;

	if (temp < 0.0) temp = 0.0;
	temp = sqrt (temp);

	vh->maxcurvature_Meyer = vh->meancurvature_Meyer + temp;
	vh->mincurvature_Meyer = vh->meancurvature_Meyer - temp;

	return true;
}

void MPMesh::ComputeVertexPrincipalCurvature_Meyer()
{
	for(unsigned int i = 0; i < pVertexList.size(); i ++)
	{
		if(is_border(pVertexList[i]) == false)
			vertex_principal_curvature_Meyer(pVertexList[i]);
	}

	m_max_maxcurvature_Meyer = -1e20;
	m_min_maxcurvature_Meyer = 1e20;
	m_max_mincurvature_Meyer = -1e20;
	m_min_mincurvature_Meyer = 1e20;
	for(unsigned int i = 0; i < pVertexList.size(); i ++)
	{
		if(is_border(pVertexList[i]) == false)
		{
			if(pVertexList[i]->maxcurvature_Meyer > m_max_maxcurvature_Meyer)
				m_max_maxcurvature_Meyer = pVertexList[i]->maxcurvature_Meyer;
			if(pVertexList[i]->maxcurvature_Meyer < m_min_maxcurvature_Meyer)
				m_min_maxcurvature_Meyer = pVertexList[i]->maxcurvature_Meyer;

			if(pVertexList[i]->mincurvature_Meyer > m_max_mincurvature_Meyer)
				m_max_mincurvature_Meyer = pVertexList[i]->mincurvature_Meyer;
			if(pVertexList[i]->mincurvature_Meyer < m_min_mincurvature_Meyer)
				m_min_mincurvature_Meyer = pVertexList[i]->mincurvature_Meyer;
		}
	}
}

void MPMesh::EstimateNormalCurvature_Meyer()
{
	ComputeVertexGaussianCurvature_Meyer();
	ComputeVertexGaussianCurvature_new();
	ComputeVertexMeanCurvatureNormal_Meyer();
	ComputeVertexPrincipalCurvature_Meyer();
}

void MPMesh::MeyerEstimationReverseOrientation()
{
	for(unsigned int i = 0; i < pVertexList.size(); i ++)
	{
		if(is_border(pVertexList[i]) == false)
		{
			pVertexList[i]->normal_Meyer = -pVertexList[i]->normal_Meyer;
			pVertexList[i]->meancurvature_Meyer = -pVertexList[i]->meancurvature_Meyer;

			Wm4::Vector3d oldmaxd, oldmind;
			oldmaxd = pVertexList[i]->maxdir_Meyer;
			oldmind = pVertexList[i]->mindir_Meyer;
			pVertexList[i]->maxdir_Meyer = oldmind;
			pVertexList[i]->mindir_Meyer = oldmaxd;

			double oldmaxc, oldminc;
			oldmaxc = pVertexList[i]->maxcurvature_Meyer; 
			oldminc = pVertexList[i]->mincurvature_Meyer; 
			pVertexList[i]->maxcurvature = -oldminc;
			pVertexList[i]->mincurvature = -oldmaxc;
			
		}
	}
}

unsigned MPMesh::NearestVertexId(Vector3d p)
{
	Triangulation::Vertex_handle vh = dt.nearest_vertex(Triangulation::Point(p[0],p[1],p[2]));
	return vh->info().id;
}

Vector3d MPMesh::NearestVertex(Vector3d p)
{
	Triangulation::Vertex_handle vh = dt.nearest_vertex(Triangulation::Point(p[0],p[1],p[2]));
	return to_wm4(vh->point());
}

Vector3d MPMesh::NearestPoint(Vector3d p)
{
	/*
	{
		Vertex_handle_t vht = dt.nearest_vertex(Point_t(p[0],p[1],p[2]));
		Vector3d fp(vht->point()[0], vht->point()[1], vht->point()[2]);
		return fp;
	}
	{
		Point_t fp(0., 0., 0.);
		fp = domain->project_on_surface(Point_t(p[0],p[1],p[2]));
		return Vector3d(fp[0], fp[1], fp[2]);
	}
	*/

	//unsigned vid = NearestVertexId(p);
	Vertex_handle_t vht = dt.nearest_vertex(Point_t(p[0],p[1],p[2]));
	unsigned vid = vht->info().id;

	MPMesh::Halfedge_around_vertex_circulator havc = pVertexList[vid]->vertex_begin();
	double mind(1e20);
	Vector3d fp;
	do
	{
		MPMesh::Face_handle fh = havc->facet(); 
		Vector3d v[3];
		MPMesh::Halfedge_around_facet_circulator hafc = fh->facet_begin();
		unsigned j = 0;
		do
		{
			v[j++] = to_wm4(hafc->vertex()->point());
			hafc ++;
		}while(hafc != fh->facet_begin());

		Vector3d tfp;
		double td;
		ProjectOntoTriangle(p, v[0], v[1], v[2], tfp, td);
		if(td < mind)
		{
			mind = td;
			fp = tfp;
		}

		havc ++;
	}while(havc != pVertexList[vid]->vertex_begin());
	return fp;
}

void MPMesh::computedt()
{
	// compute dt
	dt.clear();

	Vertex_iterator pVertex;
	int idx = 0;
	pVertex = vertices_begin();
	int my_counter = 0;
	for(; pVertex != vertices_end(); pVertex ++, idx ++)
	{
		Point_t p(pVertex->point()[0], pVertex->point()[1], pVertex->point()[2]);
		//Point_t p2(pVertex->point()[0]+0.001, pVertex->point()[1]-0.001, pVertex->point()[2]+0.001);
		Vertex_handle_t vh;
		//Vertex_handle_t vh2;
		vh = dt.insert(p);
		//vh2 = dt.insert(p2);
		vh->info().id = idx*2;
		//vh2->info().id = idx*2+1;

		my_counter++;
	}
	cout << "total number: " <<my_counter << endl;

	//for(unsigned int i = 0; i < pVertexList.size(); i ++)
	//{
	//	Point_t p(pVertexList[i]->point()[0],pVertexList[i]->point()[1],pVertexList[i]->point()[2]);
	//	Vertex_handle_t vh;
	//	vh = dt.insert(p);
	//	vh->info().id = pVertexList[i]->id;
	//}

	unsigned int fid(0);
	for(Finite_cells_iterator_t fci = dt.finite_cells_begin(); 
		fci != dt.finite_cells_end();
		fci ++)
	{
		fci->info().id = fid ++;
		Point_t cent = CGAL::circumcenter(dt.tetrahedron(fci));
		//Point_t fp = domain->project_on_surface(cent);
		//fci->info().dist_center_to_boundary = sqrt((cent-fp).squared_length());
		
		Wm4::Vector3d cent_wm4(cent[0], cent[1], cent[2]);
		if(!inside_boundingbox(cent_wm4))
			fci->info().inside = false;
		else
			fci->info().inside = (domain->is_in_domain_object()(cent) > 0);
	}
}

void MPMesh::computesimpledt()
{
	// compute dt
	dt.clear();
	for(unsigned int i = 0; i < pVertexList.size(); i ++)
	{
		Point_t p(pVertexList[i]->point()[0],pVertexList[i]->point()[1],pVertexList[i]->point()[2]);
		Vertex_handle_t vh;
		vh = dt.insert(p);
		vh->info().id = pVertexList[i]->id;
	}
}

void MPMesh::markpoles()
{
	for(Finite_cells_iterator_t fci = dt.finite_cells_begin(); fci != dt.finite_cells_end(); fci ++)
		fci->info().is_pole = false;

	for(Finite_vertices_iterator_t fvi = dt.finite_vertices_begin(); fvi != dt.finite_vertices_end(); fvi ++)
	{
		Vector3d p = to_wm4(fvi->point());
		std::vector<Cell_handle_t> fic;
		dt.finite_incident_cells(fvi,std::back_inserter(fic)); 
		double ld(0);
		Cell_handle_t ld_ch;
		bool found = false;
		for(unsigned i = 0; i < fic.size(); i ++)
		{
			if(fic[i]->info().inside)
			{
				Vector3d cp = to_wm4(fic[i]->circumcenter());
				double td = (p-cp).SquaredLength();
				if(td > ld)
				{
					ld = td;
					ld_ch = fic[i];
					found = true;
				}
			}
		}
		if(found)
		{
			ld_ch->info().is_pole = true;
			ld_ch->info().pole_bplist.insert(fvi->info().id);
		}
	}

	for(Finite_cells_iterator_t fci = dt.finite_cells_begin(); fci != dt.finite_cells_end(); fci ++)
	{
		double radius = dt.TetCircumRadius(dt.tetrahedron(fci));
		if(fci->info().dist_center_to_boundary < 0.8*radius)
			fci->info().is_pole = false;
	}
}

int MPMesh::LocalFlipCount(Vertex_handle vh)
{
	if(is_border(vh) == true)
		return -1;

	Halfedge_around_vertex_circulator havc = vh->vertex_begin();
	Wm4::Vector3d v = to_wm4(vh->point());

	std::vector<int> signv;
	do
	{
		Wm4::Vector3d oppv_v = to_wm4(havc->next()->opposite()->next()->vertex()->point()) - v;

		if(oppv_v.Dot(to_wm4(havc->facet()->normal)) > 0)
			signv.push_back(1);
		else
			signv.push_back(-1);
		havc ++;
	}
	while(havc != vh->vertex_begin());

	int count(0);
	for(int i = 0; i < signv.size(); i ++)
	{
		if(signv[i]*signv[(i+1)%(signv.size())] == -1)
			count ++;
	}

	return count;

}






