#include "medialaxissimplification3d.h"
#include <ctime>

// construction function
MedialAxisSimplification3D::MedialAxisSimplification3D(QWidget *parent, Qt::WindowFlags flags)
	: QMainWindow(parent, flags), m_isSimplified(false)
{
	ui.setupUi(this);

	m_pThreeDimensionalShape = new ThreeDimensionalShape;
	m_pCentralWidget = new QWidget;
	setCentralWidget(m_pCentralWidget);

	m_pGLWidget = new GLWidget(0);
	m_pCentralArea = new QScrollArea;
	m_pCentralArea->setWidget(m_pGLWidget);
	m_pCentralArea->setWidgetResizable(true);
	m_pCentralArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	m_pCentralArea->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	m_pCentralArea->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	m_pCentralArea->setMinimumSize(200,200);

	createActions();
	//createMenus();

	m_pCentralLayout = new QGridLayout;
	m_pCentralLayout->addWidget(m_pCentralArea, 0, 0);
	m_pCentralWidget->setLayout(m_pCentralLayout);
	statusBar()->showMessage(tr("Ready"));
	setWindowTitle(tr("Medial Axis Simplification 3D -- Computing"));
	setMinimumSize(800,600);

	//m_pThreeDimensionalShape->start_multiple = 0.02;
	//m_pThreeDimensionalShape->end_multiple = 0.05;
	m_pThreeDimensionalShape->slab_initial = false;
	slab_initial = false;
	k = 0.00001;
}

MedialAxisSimplification3D::~MedialAxisSimplification3D()
{
}

void MedialAxisSimplification3D::updateviews()
{
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::openmeshfile()
{

	// Zhiyang
	// Modifying...

	QString filename = QFileDialog::getOpenFileName(this, tr("Select a 3D model to open"), NULL, tr("3D model(*.off)"));

	if(!filename.isEmpty())
	{
		QDir qd(filename);
		QString ext = filename.right(3).toLower();
		std::cout << "Loading file " << filename.toLocal8Bit().constData() << std::endl;
		ThreeDimensionalShape * pThreeDimensionalShape = new ThreeDimensionalShape;
		bool suc = false;
		std::ifstream stream( filename.toLatin1() );
		if(stream)
		{
			stream >> pThreeDimensionalShape->input;

			// compute the properties of the input mesh
			pThreeDimensionalShape->input.computebb();		// bounding box
			pThreeDimensionalShape->input.GenerateList();		// generate vertex and triangle list
			pThreeDimensionalShape->input.GenerateRandomColor();	// color of vertex and triangle
			//pThreeDimensionalShape->input.compute_normals();		// normal of vertex and triangle
			//pThreeDimensionalShape->input.compute_sphere_matrix();  // compute the related matrix

			std::ifstream streampol(  filename.toLatin1() );
			Polyhedron pol;
			streampol >> pol;

			// set the non manifold mesh
			Mesh_domain * pdom;
			pdom = new Mesh_domain(pol);
			pThreeDimensionalShape->input.domain = pdom;

			// Computing DT and MA 
			pThreeDimensionalShape->input.computedt();
			//pThreeDimensionalShape->input.markpoles();
			pThreeDimensionalShape->ComputeInputNMM();
			//pThreeDimensionalShape->PruningSlabMesh();

			suc = true;
		}
		if(suc)
		{
			if(m_pThreeDimensionalShape != NULL)
			{
				delete m_pThreeDimensionalShape;
				m_pThreeDimensionalShape = NULL;
			}
			m_pThreeDimensionalShape = pThreeDimensionalShape;

			m_pThreeDimensionalShape->input_nmm.domain = m_pThreeDimensionalShape->input.domain;
			m_pThreeDimensionalShape->input_nmm.pmesh = &(m_pThreeDimensionalShape->input);
			m_pThreeDimensionalShape->slab_mesh.pmesh = &(m_pThreeDimensionalShape->input);
			m_pThreeDimensionalShape->slab_mesh.type = 1;
			loadconfigurefile();

			// show the input mesh in the dialog.
			m_pGLWidget->set3DShape(m_pThreeDimensionalShape);
			statusBar()->showMessage(filename + tr(" is loaded successfully.") );
			setWindowTitle( tr("Medial Axis Simplification 3D - ") + filename );
			ui.actionShow_MAT->setChecked(true);
		}
		else
		{
			if(m_pThreeDimensionalShape != NULL)
			{
				delete m_pThreeDimensionalShape;
				m_pThreeDimensionalShape = NULL;
			}
		}
	}
}

void MedialAxisSimplification3D::showmeansquareerror(bool _val)
{
	m_pGLWidget->m_display_mean_squre_error = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::openmeshfile2(){

	QString filename = QFileDialog::getOpenFileName(this, tr("Select a 3D model to open"), NULL, tr("3D model(*.off)"));
	m_offFileName = filename;
	if(!filename.isEmpty())
	{
		QDir qd(filename);
		QString prefix = filename.left(filename.size() - 4);
		std::cout << "Loading file " << filename.toLocal8Bit().constData() << std::endl;
		cout << "At medialaxissimplification3d.cpp -> openmeshfile2()!" << endl;
		ThreeDimensionalShape * pThreeDimensionalShape = new ThreeDimensionalShape;
		bool suc = false;
		std::ifstream stream( filename.toLatin1() );

		if(stream)
		{
			stream >> pThreeDimensionalShape->input;
			// compute the properties of the input mesh
			pThreeDimensionalShape->input.computebb();		// bounding box
			pThreeDimensionalShape->input.GenerateList();		// generate vertex and triangle list
			pThreeDimensionalShape->input.GenerateRandomColor();	// color of vertex and triangle
			pThreeDimensionalShape->input.compute_normals();		// normal of vertex and triangle
			pThreeDimensionalShape->input_nmm.meshname = filename.toLocal8Bit().constData();

			std::ifstream streampol(  filename.toLatin1() );
			Polyhedron pol;
			streampol >> pol;

			// set the non manifold mesh
			Mesh_domain * pdom;
			pdom = new Mesh_domain(pol);			
			suc = true;
		}
		if(suc)
		{
			if(m_pThreeDimensionalShape != NULL)
			{
				delete m_pThreeDimensionalShape;
				m_pThreeDimensionalShape = NULL;
			}
			m_pThreeDimensionalShape = pThreeDimensionalShape;			
			//ui.actionShow_Edge->setChecked(true);

			ui.actionShow_Face->setChecked(true);
			ui.actionReverse_Orientation->setChecked(true);
			

			importVP(prefix); // load the view point...
			bool re = importMA(prefix);
			if (re == false)
				return;

			m_pThreeDimensionalShape->slab_mesh.k = k;
			initialize();
			long ti = m_pThreeDimensionalShape->LoadSlabMesh();
			slab_initial = true;


			cout << "ZY: finish loading both mesh and MA..." << endl;

			// set back
			ui.actionSet_k_value->setChecked(false);			
			
			// show the input mesh in the dialog.
			m_pGLWidget->set3DShape(m_pThreeDimensionalShape);
			statusBar()->showMessage(filename + tr(" is loaded successfully.") );
			setWindowTitle( tr("Medial Axis Simplification 3D - ") + filename );

			m_isSimplified = false;
		}
		else
		{
			if(m_pThreeDimensionalShape != NULL)
			{
				delete m_pThreeDimensionalShape;
				m_pThreeDimensionalShape = NULL;
			}
		}
	}
}

void MedialAxisSimplification3D::showface(bool _val)
{
	m_pGLWidget->m_display_faces = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::showedge(bool _val)
{
	m_pGLWidget->m_display_edges = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::showvertex(bool _val)
{
	m_pGLWidget->m_display_vertices = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::showfacetnormal(bool _val){
	m_pGLWidget->m_display_facet_normals = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::showdt(bool _val)
{
	m_pGLWidget->m_display_input_dt = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::showvoronoi(bool _val)
{
	m_pGLWidget->m_display_input_voronoi = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::showmat(bool _val)
{
	m_pGLWidget->m_display_ma_edges= _val;
	m_pGLWidget->m_display_ma_faces= _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::showsphere(bool _val)
{
	m_pGLWidget->m_display_ma_spheres = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::showenvelop(bool _val)
{
	m_pGLWidget->m_display_ma_envelopes = _val;
	m_pGLWidget->updateGL();
}


void MedialAxisSimplification3D::showslabvertex(bool _val)
{
	m_pGLWidget->m_display_slab_vertex = _val;
	m_pGLWidget->updateGL();
}
void MedialAxisSimplification3D::showslabedge(bool _val)
{
	m_pGLWidget->m_display_slab_edge = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::showslabface(bool _val)
{
	m_pGLWidget->m_display_slab_faces= _val;
	m_pGLWidget->m_display_slab_edge = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::showslabenvelop(bool _val){
	m_pGLWidget->m_display_slab_envelop = _val;
	m_pGLWidget->updateGL();
}
void MedialAxisSimplification3D::showslabsphere(bool _val){
	m_pGLWidget->m_display_slab_sphere = _val;
	ui.actionSmooth_Render->setChecked(true);
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::simplifySlabByOne()
{
	if(m_pThreeDimensionalShape == NULL)
		return;
	SlabMesh *slabMesh = &(m_pThreeDimensionalShape->slab_mesh);
	int currentVertexCount = slabMesh->numVertices;
	if (currentVertexCount <= 100)
	{
		slabMesh->Simplify(5);	
		slabMesh->ComputeFacesNormal();
		slabMesh->ComputeVerticesNormal();
	}else
	{
		const int framesPerOrder = 15;
		float t = log10f(currentVertexCount);
		t = floorf(t);
		t = powf(10, t)*0.9f/framesPerOrder;
		slabMesh->Simplify(t);	
		slabMesh->ComputeFacesNormal();
		slabMesh->ComputeVerticesNormal();
	}
}

void MedialAxisSimplification3D::simplifySlab(){

	// ZY: This is for simplification...
	// before this there should be a process for weighting the edge.
	if(m_pThreeDimensionalShape == NULL || slab_initial == false)
		return;
	SlabMesh *slabMesh = &(m_pThreeDimensionalShape->slab_mesh);
	slabMesh->CleanIsolatedVertices();
	int threhold = 1;
	bool ok = FALSE;
	int simplifyNum = QInputDialog::getInt(this, tr( "Medial Mesh Simplify" ),
		tr( "Target number of vertice" ), min(10000, (int)(slabMesh->numVertices / 2)), 1, slabMesh->numVertices, 1, &ok);
	cout << simplifyNum << endl;
	if (ok)
	{
		
		threhold = simplifyNum;
	}
	else
	{
		return;
	}
	if(slabMesh == NULL)
		return;

	long start_time = clock();
	/*
	This is the core part of QMAT simplification...
	Zhiyang

	Here we make the following modification...
	- Get selected inside poles.
	- Remove the point to a better postion for further simplication.
	*/
	vector<vector<double> > selected_pole;
	cout << "Simplify_with_Selected_Pole" << endl;
	//int sss;
	//cin >> sss;
	slabMesh->Simplify_with_Selected_Pole(slabMesh->numVertices - simplifyNum, selected_pole);
	// This is the key method...
	//  To be improve
	long end_time = clock();

	m_pGLWidget->set3DSelected_Pole(selected_pole);
	m_pGLWidget->m_display_select_pole = true;
	cout << "szie selected_pole " << selected_pole.size() << endl;
	cout << "m_pGLWidget->m_display_select_pole  " << m_pGLWidget->m_display_select_pole << endl;
	slabMesh->Export_OBJ("C://Users//frank//Desktop//export//sim_MA");
// here!
	string res;
	std::stringstream ss;
	ss << end_time - start_time;
	ss >> res;

	slabMesh->ComputeFacesNormal();
	slabMesh->ComputeVerticesNormal();
	slabMesh->ComputeEdgesCone();
	slabMesh->ComputeFacesSimpleTriangles();
	m_pGLWidget->updateGL();

	QMessageBox::StandardButton reply;
	reply = QMessageBox::information(this, tr("Time"), tr(res.c_str()) + tr("ms"));
	m_isSimplified = true;
}

void MedialAxisSimplification3D::initialslab()
{
	if(m_pThreeDimensionalShape == NULL)
		return;
	SlabMesh *slabMesh = &(m_pThreeDimensionalShape->slab_mesh);
	if(slabMesh == NULL)
		return;
	if(slabMesh->numVertices == 0)
		return;
	if (ui.actionSet_k_value->isChecked() == false)
		return;

	if (slab_initial == false)
	{
		initialize();
		long ti = m_pThreeDimensionalShape->LoadSlabMesh();
		string res;
		std::stringstream ss;
		ss << ti;
		ss >> res;
		QMessageBox::StandardButton reply;
		reply = QMessageBox::information(this, tr("Initial Time"), tr(res.c_str()) + tr("ms"));
		slab_initial = true;
	}
}


void MedialAxisSimplification3D::showstabilityratio(bool _val)
{
	m_pGLWidget->m_display_stability_ratio = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::determineVertexType()
{
	if(m_pThreeDimensionalShape == NULL || slab_initial == false)
		return;
	SlabMesh *slabMesh = &(m_pThreeDimensionalShape->slab_mesh);
	if(slabMesh == NULL)
		return;

	slabMesh->DistinguishVertexType();
	slabMesh->initBoundaryCollapseQueue();
}

bool MedialAxisSimplification3D::importMA(QString name) {
	QString filename = name.append(".ma");
	QFileInfo file(filename);
	if(file.exists() == false){
		QMessageBox::StandardButton reply;
		//reply = QMessageBox::information(this, tr("Error"), tr("Related .ma file is missing."));
		//reply = QMessageBox::information(this, tr("Warning"), tr("Related .ma file is missing.\nGenerate .ma file?"));
		reply = QMessageBox::question(this, tr("Related .ma file is missing."), tr("Related .ma file is missing.\nGenerate .ma file?"), 
			QMessageBox::Yes | QMessageBox::No, QMessageBox::Yes);
		if (reply == QMessageBox::Yes)
		{
			// compute and export .ma file
			computeMA();
			//return true;
		}
		else
		{
			return false;
		}
		
	}
	
	//else
	{
		std::cout << "Loading file " << filename.toLocal8Bit().constData() << std::endl;
		bool suc = false;
		{
			m_pThreeDimensionalShape->input_nmm.meshname = filename.toLocal8Bit().constData();
			m_pThreeDimensionalShape->input_nmm.domain = m_pThreeDimensionalShape->input.domain;
			m_pThreeDimensionalShape->input_nmm.pmesh = &(m_pThreeDimensionalShape->input);
			m_pThreeDimensionalShape->slab_mesh.pmesh = &(m_pThreeDimensionalShape->input);
			m_pThreeDimensionalShape->slab_mesh.type = 1;
			m_pThreeDimensionalShape->slab_mesh.bound_weight = 1.0;
			loadconfigurefile();
			m_pThreeDimensionalShape->LoadInputNMM(filename.toLocal8Bit().constData());

			std::cout << "Done." << std::endl;
			suc = true;
		}
		if(suc)
		{
			// show the input mesh in the dialog.
			m_pGLWidget->set3DShape(m_pThreeDimensionalShape);
			//ui.actionShow_MAT->setChecked(true);
			//ui.actionSlab_Edge->setChecked(true);
			ui.actionSlab_face->setChecked(true);
			ui.actionHeadBar->setChecked(true);
		}
	}

	return true;
}

void MedialAxisSimplification3D::importMA(){
	if (m_pThreeDimensionalShape == NULL)
		return;
	if (m_pThreeDimensionalShape->input.empty())
		return;
	QString filename = QFileDialog::getOpenFileName(this, tr("Select an MA file to open"), NULL, tr("MA(*.ma)"));
	if (!filename.isEmpty())
	{
		QDir qd(filename);
		std::cout << "Loading file " << filename.toLocal8Bit().constData() << std::endl;
		std::cout << "checking here" << endl;
		int ip;
		std::cin >> ip;
		bool suc = false;
		{
			m_pThreeDimensionalShape->input_nmm.meshname = filename.toLocal8Bit().constData();
			m_pThreeDimensionalShape->input_nmm.domain = m_pThreeDimensionalShape->input.domain;
			m_pThreeDimensionalShape->input_nmm.pmesh = &(m_pThreeDimensionalShape->input);
			m_pThreeDimensionalShape->slab_mesh.pmesh = &(m_pThreeDimensionalShape->input);
			m_pThreeDimensionalShape->slab_mesh.type = 1;
			m_pThreeDimensionalShape->slab_mesh.bound_weight = 1.0;
			loadconfigurefile();
			m_pThreeDimensionalShape->LoadInputNMM(filename.toLocal8Bit().constData());
			//m_pThreeDimensionalShape->ComputeHausdorffDistance();
			
			std::cout << "Done." << std::endl;
			suc = true;
		}
		if(suc)
		{
			// show the input mesh in the dialog.
			m_pGLWidget->set3DShape(m_pThreeDimensionalShape);
			//ui.actionShow_MAT->setChecked(true);
			ui.actionSlab_Edge->setChecked(true);
			ui.actionSlab_face->setChecked(true);
			ui.actionHeadBar->setChecked(true);
		}
	}
}

void MedialAxisSimplification3D::exportMA(){

	// computeMA()
	//computeMA();
	if (m_pThreeDimensionalShape == NULL)
		return;
	m_pThreeDimensionalShape->input_nmm.Export(m_pThreeDimensionalShape->input_nmm.meshname);
	//m_pThreeDimensionalShape->slab_mesh.Export(m_pThreeDimensionalShape->input_nmm.meshname);
}

// Export simplified mesh
void MedialAxisSimplification3D::exportSimplifiedMA()
{
	if (m_isSimplified == false)
	{
		QMessageBox::information(this, tr("Haven't been simplified."), tr("Please press \"Simplify Medial Mesh\" button,\nbefore export simplified mesh."));
		return;
	}
	// Export simplified mesh
	m_pThreeDimensionalShape->slab_mesh.Export(m_pThreeDimensionalShape->input_nmm.meshname);
	QMessageBox::information(this, tr("Export simplified medial mesh."), tr("Export simplified medial mesh.\n Done."));
}
void MedialAxisSimplification3D::computeMA(){
	if (m_pThreeDimensionalShape == NULL)
		return;
	if (m_pThreeDimensionalShape->input.empty())
		return;
	string filename = m_pThreeDimensionalShape->input_nmm.meshname;
	// filename
	//filename = filename.substr(0, filename.find('.'))+".off";
	filename = m_offFileName.toLatin1();
	std::ifstream streampol(  filename );
	Polyhedron pol;
	streampol >> pol;

	// set the non manifold mesh
	Mesh_domain * pdom;
	pdom = new Mesh_domain(pol);
	m_pThreeDimensionalShape->input.domain = pdom;
	m_pThreeDimensionalShape->input_nmm.domain = pdom;
	m_pThreeDimensionalShape->input_nmm.pmesh = &(m_pThreeDimensionalShape->input);
	m_pThreeDimensionalShape->input.computedt();		// compute the delauney trianglation
	m_pThreeDimensionalShape->input.markpoles();		// mark the poles in the surface?
	m_pThreeDimensionalShape->ComputeInputNMM();
	std::cout << "Done." << std::endl;
	
}

void MedialAxisSimplification3D::loadconfigurefile()
{
	//ifstream infile("config.ini");
	//if (infile.fail()) 
	//{
	//	QMessageBox::StandardButton reply;
	//	reply = QMessageBox::information(this, tr("Warning"), tr("Configure file load fail, Load initial value."));
	//	return;
	//}

	//string line;
	//double confi_value[3];
	//int i = 0;
	//while (infile >> line) 
	//{
	//	if (line.empty())
	//		continue;

	//	int start_pos = 0, end_pos = line.size() - 1, pos;

	//	if ((pos = line.find('=')) == -1)
	//	{
	//		continue;
	//	}
	//	string na=line.substr(pos + 1, end_pos);
	//	confi_value[i++] = atof(na.c_str());
	//}
	//m_pThreeDimensionalShape->start_multiple = confi_value[0];
	//m_pThreeDimensionalShape->end_multiple = confi_value[1];
	////Math<double>::INVERSE_TOLERANCE = confi_value[2];
	//Math<double>::INVERSE_TOLERANCE = confi_value[2];
	//Math<float>::INVERSE_TOLERANCE = confi_value[2];

	//m_pThreeDimensionalShape->slab_mesh.start_multi = m_pThreeDimensionalShape->start_multiple;
	//m_pThreeDimensionalShape->slab_mesh.end_multi = m_pThreeDimensionalShape->end_multiple;

	//infile.close();
	//return ;
}

void MedialAxisSimplification3D::initialize()
{
	//if (ui.actionMethodOne->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.preserve_boundary_method = 1;
	//else if (ui.actionMethodTwo->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.preserve_boundary_method = 2;
	//else if (ui.actionMethodThree->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.preserve_boundary_method = 3;
	//else
		m_pThreeDimensionalShape->slab_mesh.preserve_boundary_method = 0;

	//if (ui.actionTypeOne->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.hyperbolic_weight_type = 1;
	//else if (ui.actionTypeTwo->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.hyperbolic_weight_type = 2;
	//else if (ui.actionTypeThree->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.hyperbolic_weight_type = 3;
	//else
		m_pThreeDimensionalShape->slab_mesh.hyperbolic_weight_type = 3;

	//if (ui.actionPreserve_Saved_Vertex->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.preserve_saved_vertices = true;
	//else
	//	m_pThreeDimensionalShape->slab_mesh.preserve_saved_vertices = false;

	//if (ui.actionComputer_Hausdorff->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.compute_hausdorff = true;
	//else
		m_pThreeDimensionalShape->slab_mesh.compute_hausdorff = false;

	//if (ui.actionClear_Error->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.clear_error = true;
	//else
	//	m_pThreeDimensionalShape->slab_mesh.clear_error = false;


	//if (ui.actionScale_One->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.boundary_compute_scale = 1;
	//else if (ui.actionScale_Two->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.boundary_compute_scale = 2;
	//else if (ui.action_Other_Scale->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.boundary_compute_scale = 3;
	//else
		m_pThreeDimensionalShape->slab_mesh.boundary_compute_scale = 0;

	//if (ui.actionPrevent_Inversion->isChecked())
	//	m_pThreeDimensionalShape->slab_mesh.prevent_inversion = true;
	//else
		m_pThreeDimensionalShape->slab_mesh.prevent_inversion = false;
}

void MedialAxisSimplification3D::zoomma(bool _val)
{
	m_pGLWidget->m_display_zoom_ma = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::setratiofactor()
{
	if(m_pThreeDimensionalShape == NULL || m_pThreeDimensionalShape->slab_mesh.numVertices == 0) {
		return;
	}


	double tk = 0.001;
	bool ok = FALSE;
	double factor = QInputDialog::getDouble(this,
		tr( "set ratio factor" ),
		tr( "Please set 'k' value: " ), k, 0, 1, 8, &ok);
	if (ok)
		tk = factor;
	else
		return;

	k = tk;

	QMessageBox::StandardButton reply;
	reply = QMessageBox::information(this, tr("Time"), tr("The setting will be valid in the next loading."));
}

void MedialAxisSimplification3D::saveEPS()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Save Screen as PNG"), NULL, tr("PNG(*.png)"));
	if (!filename.isEmpty()) 
	{
		string file_totle_add = filename.toLocal8Bit().constData();
		int po = file_totle_add.find_last_of('/');
		string addre = file_totle_add.substr(0, po + 1);

		if(m_pThreeDimensionalShape == NULL)
			return;

		SlabMesh *slabMesh = &(m_pThreeDimensionalShape->slab_mesh);

		//int number = slabMesh->numVertices / 50;
		//int offset = 1;
		//while(number / 50 > 0)
		//{
		//	offset++;
		//	number /= 50;
		//}

		int index = 0;

		//while (slabMesh->numVertices > 20)
		//{
			stringstream ss;
			ss << index;
			string file_na=ss.str();
			
			//if (file_na.size() < offset)
			//{
			//	string att = "";
			//	string add_0("0");
			//	int num0 = offset - file_na.size();
			//	while(num0 > 0)
			//	{
			//		file_na = add_0 + file_na;
			//		num0--;
			//	}
			//}

			string ind(".png");
			string ful_name =  addre + file_na + ind;
			/*
			const int size = (int)9e6;
			GLfloat *pFeedbackBuffer = new GLfloat[size];
			glFeedbackBuffer(size,GL_3D_COLOR,pFeedbackBuffer);
			glRenderMode(GL_FEEDBACK);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glPushMatrix();
			m_pGLWidget->updateGL();	
			glPopMatrix();
			glFlush();
			int NbValues = glRenderMode(GL_RENDER);
			CPsRenderer PsRenderer;
			PsRenderer.Run(ful_name.c_str(),pFeedbackBuffer,NbValues,true);
			delete [] pFeedbackBuffer;*/

			glRenderMode(GL_RENDER);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glPushMatrix();
			m_pGLWidget->updateGL();	
			glPopMatrix();
			glFlush();

			QPixmap p = QPixmap::grabWidget(m_pGLWidget);
			p.save(ful_name.c_str());

		//	index++;

		//	simplifySlabByOne();

		//}
		
	}
}

void MedialAxisSimplification3D::saveViewPoint()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Save current viewpoint"), NULL, tr("Viewpoint(*.viewpoint)"));
	if (!filename.isEmpty())
		m_pGLWidget->SaveViewpoint(filename.toLatin1());	
	m_pGLWidget->updateGL();	
}

void MedialAxisSimplification3D::importVP(QString viewpoint){
	QString filename = viewpoint.append(".viewpoint");
	QFileInfo file(filename);
	if(file.exists() == true){
		m_pGLWidget->LoadViewpoint(filename.toLatin1());
	}
}

void MedialAxisSimplification3D::loadViewPoint()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Select a viewpoint to open"), NULL, tr("Viewpoint(*.viewpoint)"));
	if(!filename.isEmpty())
		m_pGLWidget->LoadViewpoint(filename.toLatin1());
	m_pGLWidget->updateGL();	
}

void MedialAxisSimplification3D::setCheckReverseOrientation(bool _val)
{
	m_pGLWidget->m_display_reverse_orientation = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::setCheckHeadBar(bool _val)
{
	m_pGLWidget->m_display_headbar = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::setCheckColorBar(bool _val)
{
	m_pGLWidget->m_display_colorbar = _val;
	m_pGLWidget->updateGL();
}

void MedialAxisSimplification3D::setSmoothRendering(bool _val)
{
	m_pGLWidget->m_smoothshading = _val;
	m_pGLWidget->updateGL();
}

// combine each operate to the special function
void MedialAxisSimplification3D::createActions()
{
	connect(ui.actionOpen_Off, SIGNAL(triggered()), this, SLOT(openmeshfile()));
	connect(ui.actionOpen_Off_simple, SIGNAL(triggered()), this, SLOT(openmeshfile2()));
	connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
	connect(ui.actionShow_Face,SIGNAL(toggled(bool)),this,SLOT(showface(bool)));
	connect(ui.actionShow_Edge,SIGNAL(toggled(bool)),this,SLOT(showedge(bool)));
	connect(ui.actionShow_Vertex,SIGNAL(toggled(bool)),this,SLOT(showvertex(bool)));
	connect(ui.actionDT,SIGNAL(toggled(bool)),this,SLOT(showdt(bool)));
	connect(ui.actionVoronoi_Diagram,SIGNAL(toggled(bool)),this,SLOT(showvoronoi(bool)));
	connect(ui.actionShow_MAT,SIGNAL(toggled(bool)),this,SLOT(showmat(bool)));
	connect(ui.actionShow_Sphere,SIGNAL(toggled(bool)),this,SLOT(showsphere(bool)));
	connect(ui.actionShow_Envelop,SIGNAL(toggled(bool)),this,SLOT(showenvelop(bool)));

	connect(ui.actionShow_Facet_Normal,SIGNAL(toggled(bool)),this,SLOT(showfacetnormal(bool)));

	connect(ui.actionImport_Medial_Axis, SIGNAL(triggered()), this, SLOT(importMA()));
	connect(ui.actionExport_Medial_Axis, SIGNAL(triggered()), this, SLOT(exportSimplifiedMA()));
	connect(ui.actionCompute_Medial_Axis, SIGNAL(triggered()), this, SLOT(computeMA()));


	connect(ui.actionSlab_Vertex,SIGNAL(toggled(bool)), this, SLOT(showslabvertex(bool)));
	connect(ui.actionSlab_Edge,SIGNAL(toggled(bool)), this, SLOT(showslabedge(bool)));
	connect(ui.actionSlab_face,SIGNAL(toggled(bool)), this, SLOT(showslabface(bool)));

	connect(ui.actionSlab_Envelop,SIGNAL(toggled(bool)), this, SLOT(showslabenvelop(bool)));
	connect(ui.actionSlab_Sphere,SIGNAL(toggled(bool)), this, SLOT(showslabsphere(bool)));
	connect(ui.actionSimplify_Slab, SIGNAL(triggered()), this, SLOT(simplifySlab()));
	connect(ui.actionStability_Ratio, SIGNAL(toggled(bool)), this, SLOT(showmeansquareerror(bool)));

	connect(ui.actionShow_stability_ratio, SIGNAL(toggled(bool)), this, SLOT(showstabilityratio(bool)));

	connect(ui.actionInitial_Slab, SIGNAL(triggered()), this, SLOT(initialslab()));

	connect(ui.actionDetermine_Vertex_Type, SIGNAL(triggered()), this, SLOT(determineVertexType()));
	connect(ui.actionInitialize, SIGNAL(triggered()), this, SLOT(initialize()));
	connect(ui.actionZoom_MA, SIGNAL(toggled(bool)), this, SLOT(zoomma(bool)));

	connect(ui.actionSet_k_value, SIGNAL(triggered()), this, SLOT(setratiofactor()));
	
	connect(ui.actionSave_EPS,  SIGNAL(triggered()), this, SLOT(saveEPS()));
	connect(ui.actionSave_ViewPoint,  SIGNAL(triggered()), this, SLOT(saveViewPoint()));
	connect(ui.actionLoad_ViewPoint,  SIGNAL(triggered()), this, SLOT(loadViewPoint()));
	connect(ui.actionReverse_Orientation,  SIGNAL(toggled(bool)), this, SLOT(setCheckReverseOrientation(bool)));
	connect(ui.actionHeadBar,  SIGNAL(toggled(bool)), this, SLOT(setCheckHeadBar(bool)));
	connect(ui.actionColorBar,  SIGNAL(toggled(bool)), this, SLOT(setCheckColorBar(bool)));
	connect(ui.actionSmooth_Render,  SIGNAL(toggled(bool)), this, SLOT(setSmoothRendering(bool)));
}

