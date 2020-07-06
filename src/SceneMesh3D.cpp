#include "SceneMesh3D.h"
#include "ctime"

#define FILENAME "hand2"
#define CONTROL_POINTS 1
#define SAMPLE_POINTS 50
#define SECTION_SIZE 0.5
#define NUM_OF_FILES 18

#define COMPUTE_EIGENDECOMP 0
#define COMPUTE_COORDS 0
#define COMPUTE_GEODESIC 0

using namespace std;
using namespace vvr;
using namespace Eigen;

Mesh3DScene::Mesh3DScene()
{
	//! Load settings.
	vvr::Shape::DEF_LINE_WIDTH = 4;
	vvr::Shape::DEF_POINT_SIZE = 10;
	m_perspective_proj = true;
	m_bg_col = Colour("768E77");
	m_obj_col = Colour("454545");
	const string objDir = getBasePath() + "resources/obj/";
	const string objFile = objDir + FILENAME + ".obj";
	m_model_original = vvr::Mesh(objFile);
	m_model_newA = vvr::Mesh(m_model_original);
	m_model_newB = vvr::Mesh(m_model_original);
	m_model_newC = vvr::Mesh(m_model_original);
	m_model_newD = vvr::Mesh(m_model_original);
	m_model_new_draw = &m_model_newA;
	m_points_coords3D_draw = &m_points_coords3DA;
	m_points_difCoords3D_draw = &m_points_difCoords3DA;
	m_triangles_coords3D_draw = &m_triangles_coords3DA;
	m_triangles_difCoords3D_draw = &m_triangles_difCoords3DA;
	originalFiles = { "armadillo_low_low", "b66_L2", "bone", "bunny_low", "cube", "dolphin", "dragon_low_low",
		"flashlight", "flashlightNoCentered", "hand2", "icosahedron", "phone_v02", "polyhedron",
		"suzanne", "teapotMultiMesh", "unicorn_low", "unicorn_low_low", "vvrlab" };
	reset();
}

void Mesh3DScene::reset()
{
	Scene::reset();
	//! Define what will be vissible by default
	m_style_flag = 0;
	ratio_flag = 0;
	m_style_flag |= FLAG_SHOW_SOLID;
	m_style_flag |= FLAG_SHOW_WIRE;
	m_style_flag |= FLAG_SHOW_AXES;
	m_style_flag |= FLAG_SHOW_AABB;
	m_style_flag |= FLAG_SHOW_ORIGINAL_MODEL;
	m_style_flag |= FLAG_SHOW_TRIANGLES;
}

void Mesh3DScene::resize()
{
	//! By Making `first_pass` static and initializing it to true,
	//! we make sure that the if block will be executed only once.

	static bool first_pass = true;

	if (first_pass)
	{
		m_model_original.setBigSize(getSceneWidth() / 2);
		m_model_original.update();
		m_model_newA.setBigSize(getSceneWidth() / 2);
		m_model_newB.setBigSize(getSceneWidth() / 2);
		m_model_newC.setBigSize(getSceneWidth() / 2);
		m_model_newD.setBigSize(getSceneWidth() / 2);
		m_model_newA.update();
		m_model_newB.update();
		m_model_newC.update();
		m_model_newD.update();
		m_model = m_model_original;
		Tasks();
		first_pass = false;
	}
}

void Mesh3DScene::Tasks()
{
	vector<vvr::Triangle>& m_triangles = m_model.getTriangles();
	vector<vec>& m_vertices = m_model.getVertices();
	int verticesCount = m_vertices.size();
	int trianglesCount = m_triangles.size();
	cout << "verticesCount = " << verticesCount << endl;

	SparseMatrix<double> A(verticesCount, verticesCount);
	SparseMatrix<double> I(verticesCount, verticesCount);
	SparseMatrix<double> D(verticesCount, verticesCount);
	SparseMatrix<double> D_inverse(verticesCount, verticesCount);
	SparseMatrix<double> L(verticesCount, verticesCount);
	MatrixXd Coords(verticesCount, 3);
	MatrixXd DifCoords(verticesCount, 3);
	Task1(m_triangles, m_vertices, I, A, D, D_inverse, L, Coords, DifCoords);

	MatrixXd Ls = D - A;
	VectorXd eigenValues(verticesCount);
	MatrixXd eigenVectors(verticesCount, verticesCount);
	Task2(Ls, eigenValues, eigenVectors, verticesCount);
	MatrixXd Q = eigenVectors;
	MatrixXd CoordsNewA(verticesCount, 3);
	MatrixXd DifCoordsNewA(verticesCount, 3);
	MatrixXd CoordsNewB(verticesCount, 3);
	MatrixXd DifCoordsNewB(verticesCount, 3);
	MatrixXd CoordsNewC(verticesCount, 3);
	MatrixXd DifCoordsNewC(verticesCount, 3);
	MatrixXd CoordsNewD(verticesCount, 3);
	MatrixXd DifCoordsNewD(verticesCount, 3);
	MatrixXd LDense = L;

	Task3(Q, Coords, DifCoords, CoordsNewA, DifCoordsNewA, CoordsNewB, DifCoordsNewB, CoordsNewC, DifCoordsNewC, CoordsNewD, DifCoordsNewD, verticesCount, L);
	Task4(Coords, CoordsNewA, CoordsNewB, CoordsNewC, CoordsNewD, DifCoords, DifCoordsNewA, DifCoordsNewB, DifCoordsNewC, DifCoordsNewD, LDense);
	Task5A(verticesCount, eigenValues);
	Task5B();
	Task5C();
}

void Mesh3DScene::Task1(vector<vvr::Triangle>& m_triangles, vector<vec>& m_vertices, SparseMatrix<double>& I, SparseMatrix<double>& A, SparseMatrix<double>& D, SparseMatrix<double>& D_inverse, SparseMatrix<double>& L, MatrixXd& Coords, MatrixXd& DifCoords) {
	int verticesCount = m_vertices.size();
	int trianglesCount = m_triangles.size();
	SparseIdentity(I, verticesCount);
	//give values to A and D
	for (int i = 0; i < trianglesCount; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = j + 1; k < 3; k++) {
				if (A.coeffRef(m_triangles[i].v[j], m_triangles[i].v[k]) == 0) {
					A.coeffRef(m_triangles[i].v[j], m_triangles[i].v[k]) = 1;
					A.coeffRef(m_triangles[i].v[k], m_triangles[i].v[j]) = 1;
					D.coeffRef(m_triangles[i].v[j], m_triangles[i].v[j])++;
					D.coeffRef(m_triangles[i].v[k], m_triangles[i].v[k])++;
				}
			}
		}
	}
	//*/
	for (int i = 0; i < verticesCount; i++) {
		for (int j = i + 1; j < verticesCount; j++) {
			if (dist(m_vertices[i], m_vertices[j]) == 0) {
				A.coeffRef(i, j) = 1;
				A.coeffRef(j, i) = 1;
				D.coeffRef(i, i)++;
				D.coeffRef(j, j)++;
			}
		}
	}
	//*/

	SparseDiagonalInverse(D, D_inverse, verticesCount);
	L = I - D_inverse * A;
	for (int i = 0; i < verticesCount; i++) {
		Coords(i, 0) = (double)m_vertices[i].x;
		Coords(i, 1) = (double)m_vertices[i].y;
		Coords(i, 2) = (double)m_vertices[i].z;
	}

	DifCoords = L * Coords;

	//find max of dif coords on each axis
	double maxX = FindMax(DifCoords, verticesCount, 0);
	double maxY = FindMax(DifCoords, verticesCount, 1);
	double maxZ = FindMax(DifCoords, verticesCount, 2);
	MatrixXd newDifCoords(verticesCount, 3);
	GetDifCoordsInNormalDirection(m_triangles, m_vertices, DifCoords, newDifCoords);
	//compute the magnitude of the dif coords of all the vertices
	VectorXd DifCoordsMagnitude(verticesCount);
	double maxMag = 0;
	for (int i = 0; i < verticesCount; i++) {
		DifCoordsMagnitude(i) = sqrt(newDifCoords(i, 0) * newDifCoords(i, 0) + newDifCoords(i, 1) * newDifCoords(i, 1) + newDifCoords(i, 2) * newDifCoords(i, 2));
		if (DifCoordsMagnitude(i) > maxMag)
			maxMag = DifCoordsMagnitude(i);
	}
	//draw all vertices with the corresponding colour
	for (int i = 0; i < verticesCount; i++) {
		m_points3D.push_back(Point3D(m_vertices[i].x, m_vertices[i].y, m_vertices[i].z, vvr::Colour((DifCoordsMagnitude(i) / maxMag) * 255, (1 - (DifCoordsMagnitude(i) / maxMag)) * 255, 0)));
	}
	//draw all triangles with the corresponding color
	for (int i = 0; i < trianglesCount; i++) {
		m_triangles3D.push_back(Triangle3D(m_vertices[m_triangles[i].vi1].x, m_vertices[m_triangles[i].vi1].y, m_vertices[m_triangles[i].vi1].z,
			m_vertices[m_triangles[i].vi2].x, m_vertices[m_triangles[i].vi2].y, m_vertices[m_triangles[i].vi2].z,
			m_vertices[m_triangles[i].vi3].x, m_vertices[m_triangles[i].vi3].y, m_vertices[m_triangles[i].vi3].z, vvr::Colour::black));
		m_triangles3D[i].setColourPerVertex(vvr::Colour((DifCoordsMagnitude(m_triangles[i].vi1) / maxMag) * 255, (1 - (DifCoordsMagnitude(m_triangles[i].vi1) / maxMag)) * 255, 0),
			vvr::Colour((DifCoordsMagnitude(m_triangles[i].vi2) / maxMag) * 255, (1 - (DifCoordsMagnitude(m_triangles[i].vi2) / maxMag)) * 255, 0),
			vvr::Colour((DifCoordsMagnitude(m_triangles[i].vi3) / maxMag) * 255, (1 - (DifCoordsMagnitude(m_triangles[i].vi3) / maxMag)) * 255, 0));
	}
	cout << "Task1 done\n" << endl;
}

void Mesh3DScene::Task2(MatrixXd& Ls, VectorXd& eigenValues, MatrixXd& eigenVectors, int verticesCount) {
	const string eigenDir = getBasePath() + "resources/eigen/";
	const string eigenFile = eigenDir + FILENAME + ".txt";
	if (COMPUTE_EIGENDECOMP) {
		ComputeEigenDecomposition(Ls, eigenValues, eigenVectors);
		//SaveEigenToFile(eigenFile, eigenValues, eigenVectors, verticesCount);
	}
	else {
		ReadEigenFromFile(eigenFile, eigenValues, eigenVectors, verticesCount);
		cout << "fileIO complete" << endl;
	}
	cout << "Task2 done\n" << endl;
}

void Mesh3DScene::Task3(MatrixXd& Q, MatrixXd& Coords, MatrixXd& DifCoords, MatrixXd& CoordsNewA, MatrixXd& DifCoordsNewA, MatrixXd& CoordsNewB, MatrixXd& DifCoordsNewB, MatrixXd& CoordsNewC, MatrixXd& DifCoordsNewC, MatrixXd& CoordsNewD, MatrixXd& DifCoordsNewD, int verticesCount, SparseMatrix<double>& L) {
	double ratios[4] = { 0.3, 0.1, 0.05, 0.01 };
	const string coordsDir = getBasePath() + "resources/coords/";
	const string coordsFile = coordsDir + FILENAME + ".txt";
	if (COMPUTE_COORDS) {
		for (int i = 0; i < 4; i++) {
			MatrixXd Qcopy = Q;
			SparseMatrix<double> LCopy = L;
			if (i == 0)
				Task3Sub(Qcopy, Coords, DifCoords, CoordsNewA, DifCoordsNewA, m_model_newA.getVertices(), verticesCount, ratios[i], LCopy);
			if (i == 1)
				Task3Sub(Qcopy, Coords, DifCoords, CoordsNewB, DifCoordsNewB, m_model_newB.getVertices(), verticesCount, ratios[i], LCopy);
			if (i == 2)
				Task3Sub(Qcopy, Coords, DifCoords, CoordsNewC, DifCoordsNewC, m_model_newC.getVertices(), verticesCount, ratios[i], LCopy);
			if (i == 3)
				Task3Sub(Qcopy, Coords, DifCoords, CoordsNewD, DifCoordsNewD, m_model_newD.getVertices(), verticesCount, ratios[i], LCopy);
		}
		//SaveCoordsToFile(coordsFile, CoordsNewA, DifCoordsNewA, CoordsNewB, DifCoordsNewB, CoordsNewC, DifCoordsNewC, CoordsNewD, DifCoordsNewD);
	}
	else {
		ReadCoordsFromFile(coordsFile, CoordsNewA, DifCoordsNewA, CoordsNewB, DifCoordsNewB, CoordsNewC, DifCoordsNewC, CoordsNewD, DifCoordsNewD, verticesCount);
	}
	cout << "Task3 done\n" << endl;
}

void Mesh3DScene::Task3Sub(MatrixXd& Q, MatrixXd& Coords, MatrixXd& DifCoords, MatrixXd& CoordsNew, MatrixXd& DifCoordsNew, vector<vec>& m_vertices_new, int verticesCount, double ratio, SparseMatrix<double>& L) {
	//*
	MatrixXd Xtilda = Q.transpose() * Coords;
	for (int i = ratio * verticesCount; i < verticesCount; i++) {
		for (int j = 0; j < 3; j++) {
			Xtilda(i, j) = 0;
		}
	}
	DifCoordsNew = L * (Q * Xtilda);
	//*/
	/*
	MatrixXd DeltaTilda = Q.transpose() * DifCoords;
	for (int i = ratio * verticesCount; i < verticesCount; i++) {
		for (int j = 0; j < 3; j++) {
			DeltaTilda(i, j) = 0;
		}
	}
	DifCoordsNew = Q * DeltaTilda;
	//*/

	//CoordsNew = Q * Xtilda;

	int* options = new int[verticesCount];
	for (int i = 0; i < verticesCount; i++) {
		options[i] = i;
	}
	shuffle_arr(options, verticesCount);

	cout << "model reconstruction using "<< 100 * ratio <<"% of eigenvectors begins:" << endl;
	//int rank;
	int dimCount = verticesCount;
	//SparseQR <SparseMatrix<double>, COLAMDOrdering<int>> col_decomp;
	SimplicialLLT<SparseMatrix<double>> llt;
	while (true) {
		dimCount++;
		int count = dimCount - verticesCount - 1;
		L.conservativeResize(dimCount, verticesCount);
		L.coeffRef(dimCount - 1, options[count]) = 1;

		DifCoordsNew.conservativeResize(dimCount, NoChange);
		for (int i = 0; i < 3; i++) {
			DifCoordsNew(dimCount - 1, i) = Coords(options[count], i);
		}

		if (!L.isCompressed())
			L.makeCompressed();
		llt.compute(L.transpose() * L);
		double det = llt.determinant();

		if (abs(det) > 0 && count >= CONTROL_POINTS - 1) {
			cout << "det = " << det << "	count = " << count <<endl;
			break;
		}
	}

	CoordsNew = llt.solve(L.transpose() * DifCoordsNew);
	cout << "model reconstructed\n" << endl;
	CoordsToModel(m_vertices_new, CoordsNew, verticesCount);
}

void Mesh3DScene::Task4(const MatrixXd& Coords, const MatrixXd& CoordsNewA, const MatrixXd& CoordsNewB, const MatrixXd& CoordsNewC, const MatrixXd& CoordsNewD, const MatrixXd& DifCoords, const MatrixXd& DifCoordsNewA, const MatrixXd& DifCoordsNewB, const MatrixXd& DifCoordsNewC, const MatrixXd& DifCoordsNewD, const MatrixXd& L) {
	vector<vvr::Triangle>& m_triangles = m_model.getTriangles();
	Task4Sub(m_triangles, m_model_newA.getVertices(), m_points_coords3DA, m_triangles_coords3DA, m_points_difCoords3DA, m_triangles_difCoords3DA, Coords, CoordsNewA, DifCoords, DifCoordsNewA, L);
	Task4Sub(m_triangles, m_model_newB.getVertices(), m_points_coords3DB, m_triangles_coords3DB, m_points_difCoords3DB, m_triangles_difCoords3DB, Coords, CoordsNewB, DifCoords, DifCoordsNewB, L);
	Task4Sub(m_triangles, m_model_newC.getVertices(), m_points_coords3DC, m_triangles_coords3DC, m_points_difCoords3DC, m_triangles_difCoords3DC, Coords, CoordsNewC, DifCoords, DifCoordsNewC, L);
	Task4Sub(m_triangles, m_model_newD.getVertices(), m_points_coords3DD, m_triangles_coords3DD, m_points_difCoords3DD, m_triangles_difCoords3DD, Coords, CoordsNewD, DifCoords, DifCoordsNewD, L);
	cout << "Task4 done\n" << endl;
}

void Mesh3DScene::Task4Sub(vector<vvr::Triangle>& m_triangles, vector<vec>& m_vertices_new, vector<Point3D>& m_points_coords3D, vector<Triangle3D>& m_triangles_coords3D, vector<Point3D>& m_points_difCoords3D, vector<Triangle3D>& m_triangles_difCoords3D, const MatrixXd& Coords, const MatrixXd& CoordsNew, const MatrixXd& DifCoords, const MatrixXd& DifCoordsNew, const MatrixXd& L) {
	int verticesCount = m_vertices_new.size();
	int trianglesCount = m_triangles.size();
	
	VectorXd Dif(verticesCount);
	double maxDif = 0;
	for (int i = 0; i < verticesCount; i++) {
		vec diference = vec(CoordsNew(i, 0) - Coords(i, 0), CoordsNew(i, 1) - Coords(i, 1), CoordsNew(i, 2) - Coords(i, 2));
		Dif(i) = diference.Length();
		if (Dif(i) > maxDif) {
			maxDif = Dif(i);
		}
	}
	for (int i = 0; i < verticesCount; i++) {
		m_points_coords3D.push_back(Point3D(m_vertices_new[i].x, m_vertices_new[i].y, m_vertices_new[i].z, vvr::Colour((Dif(i) / maxDif) * 255, (1 - (Dif(i) / maxDif)) * 255, 0)));
	}
	for (int i = 0; i < trianglesCount; i++) {
		m_triangles_coords3D.push_back(Triangle3D(m_vertices_new[m_triangles[i].vi1].x, m_vertices_new[m_triangles[i].vi1].y, m_vertices_new[m_triangles[i].vi1].z,
			m_vertices_new[m_triangles[i].vi2].x, m_vertices_new[m_triangles[i].vi2].y, m_vertices_new[m_triangles[i].vi2].z,
			m_vertices_new[m_triangles[i].vi3].x, m_vertices_new[m_triangles[i].vi3].y, m_vertices_new[m_triangles[i].vi3].z, vvr::Colour::black));
		m_triangles_coords3D[i].setColourPerVertex(vvr::Colour((Dif(m_triangles[i].vi1) / maxDif) * 255, (1 - (Dif(m_triangles[i].vi1) / maxDif)) * 255, 0),
			vvr::Colour((Dif(m_triangles[i].vi2) / maxDif) * 255, (1 - (Dif(m_triangles[i].vi2) / maxDif)) * 255, 0),
			vvr::Colour((Dif(m_triangles[i].vi3) / maxDif) * 255, (1 - (Dif(m_triangles[i].vi3) / maxDif)) * 255, 0));
	}
	
	MatrixXd DifCoordsNewNew = L * CoordsNew;
	double totalMistake = 0;
	for (int i = 0; i < verticesCount; i++) {
		vec diff = vec(DifCoordsNewNew(i, 0) - DifCoordsNew(i, 0), DifCoordsNewNew(i, 1) - DifCoordsNew(i, 1), DifCoordsNewNew(i, 2) - DifCoordsNew(i, 2));
		totalMistake += diff.Length();
	}
	cout << "total mistake = " << totalMistake << endl;
	VectorXd DifInDifCoords(verticesCount);
	double maxDifInDifCoords = 0;
	for (int i = 0; i < verticesCount; i++) {
		vec diference = vec(DifCoordsNew(i, 0) - DifCoords(i, 0), DifCoordsNew(i, 1) - DifCoords(i, 1), DifCoordsNew(i, 2) - DifCoords(i, 2));
		DifInDifCoords(i) = diference.Length();
		if (DifInDifCoords(i) > maxDifInDifCoords) {
			maxDifInDifCoords = DifInDifCoords(i);
		}
	}
	for (int i = 0; i < verticesCount; i++) {
		m_points_difCoords3D.push_back(Point3D(m_vertices_new[i].x, m_vertices_new[i].y, m_vertices_new[i].z, vvr::Colour((DifInDifCoords(i) / maxDifInDifCoords) * 255, 0, (1 - (DifInDifCoords(i) / maxDifInDifCoords)) * 255)));
	}
	for (int i = 0; i < trianglesCount; i++) {
		m_triangles_difCoords3D.push_back(Triangle3D(m_vertices_new[m_triangles[i].vi1].x, m_vertices_new[m_triangles[i].vi1].y, m_vertices_new[m_triangles[i].vi1].z,
			m_vertices_new[m_triangles[i].vi2].x, m_vertices_new[m_triangles[i].vi2].y, m_vertices_new[m_triangles[i].vi2].z,
			m_vertices_new[m_triangles[i].vi3].x, m_vertices_new[m_triangles[i].vi3].y, m_vertices_new[m_triangles[i].vi3].z, vvr::Colour::black));
		m_triangles_difCoords3D[i].setColourPerVertex(vvr::Colour((DifInDifCoords(m_triangles[i].vi1) / maxDifInDifCoords) * 255, 0, (1 - (DifInDifCoords(m_triangles[i].vi1) / maxDifInDifCoords)) * 255),
			vvr::Colour((DifInDifCoords(m_triangles[i].vi2) / maxDifInDifCoords) * 255, 0, (1 - (DifInDifCoords(m_triangles[i].vi2) / maxDifInDifCoords)) * 255),
			vvr::Colour((DifInDifCoords(m_triangles[i].vi3) / maxDifInDifCoords) * 255, 0, (1 - (DifInDifCoords(m_triangles[i].vi3) / maxDifInDifCoords)) * 255));
	}
}

void Mesh3DScene::Task5A(int verts, VectorXd& eigenVal) {
	files = originalFiles;
	const string eigenDir = getBasePath() + "resources/eigen/";
	int verticesCount[NUM_OF_FILES + 1];
	VectorXd eigenValues[NUM_OF_FILES + 1];
	double maxEigenValue = 0;
	for (int i = 0; i < NUM_OF_FILES; i++) {
		string eigenFile = eigenDir + files[i] + ".txt";
		ifstream myfile;
		myfile.open(eigenFile);
		myfile >> verticesCount[i];
		eigenValues[i] = VectorXd(verticesCount[i]);
		for (int j = 0; j < verticesCount[i]; j++) {
			myfile >> eigenValues[i](j);
		}
		myfile.close();
		if (eigenValues[i](verticesCount[i] - 1) > maxEigenValue)
			maxEigenValue = eigenValues[i](verticesCount[i] - 1);
	}

	verticesCount[NUM_OF_FILES] = verts;
	eigenValues[NUM_OF_FILES] = eigenVal;
	if (eigenVal[verts - 1] > maxEigenValue)
		maxEigenValue = eigenVal[verts - 1];
	
	int sectionCount = ceil(maxEigenValue / SECTION_SIZE);
	double** eigenSections = new double* [NUM_OF_FILES + 1];
	for (int i = 0; i < NUM_OF_FILES + 1; i++) {
		eigenSections[i] = new double[sectionCount];
		for (int j = 0; j < sectionCount; j++) {
			eigenSections[i][j] = 0;
		}
	}
	double maxSection[NUM_OF_FILES + 1];
	for (int i = 0; i < NUM_OF_FILES + 1; i++) {
		maxSection[i] = 0;
	}
	for (int i = 0; i < NUM_OF_FILES + 1; i++) {
		for (int j = 0; j < verticesCount[i]; j++) {
			int section = ceil(eigenValues[i](j) / SECTION_SIZE);
			eigenSections[i][section] ++;
		}
		for (int j = 0; j < sectionCount; j++) {
			if (eigenSections[i][j] > maxSection[i])
				maxSection[i] = eigenSections[i][j];
		}
		for (int j = 0; j < sectionCount; j++) {
			eigenSections[i][j] /= maxSection[i];          //verticesCount[i]; 
		}
	}
	double errors[NUM_OF_FILES];
	int myIndex;
	for (int i = 0; i < NUM_OF_FILES; i++) {
		errors[i] = 0;
	}
	for (int i = 0; i < NUM_OF_FILES; i++) {
		for (int j = 0; j < sectionCount; j++) {
			errors[i] += (eigenSections[NUM_OF_FILES][j] - eigenSections[i][j]) * (eigenSections[NUM_OF_FILES][j] - eigenSections[i][j]);
		}
	}
	printResult(errors);

	for (int i = 0; i < NUM_OF_FILES + 1; i++) {
		delete[] eigenSections[i];
	}
	delete[] eigenSections;
}

void Mesh3DScene::Task5B() {
	files = originalFiles;
	const string objDir = getBasePath() + "resources/obj/";
	double*** distances = new double** [NUM_OF_FILES + 1];
	for (int i = 0; i < NUM_OF_FILES + 1; i++) {
		distances[i] = new double* [SAMPLE_POINTS];
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			distances[i][j] = new double[SAMPLE_POINTS];
		}
	}
	for (int i = 0; i < NUM_OF_FILES + 1; i++) {
		string objFile;
		if (i < NUM_OF_FILES)
			objFile = objDir + files[i] + ".obj";
		else
			objFile = objDir + FILENAME + ".obj";
		m_model_other = Mesh(objFile);
		vector<vec>& vertices = m_model_other.getVertices();
		vec* verticesSelected = new vec[SAMPLE_POINTS];
		verticesSelected[0] = vertices[0];       //vertices[rand() % vertices.size()];
		for (int j = 1; j < SAMPLE_POINTS; j++) {
			double biggestScore = -1;
			int biggestIndex = -1;
			for (int k = 0; k < vertices.size(); k++) {
				double minDist = 1e5;
				for (int l = 0; l < j; l++) {
					if (dist(vertices[k], verticesSelected[l]) < minDist)
						minDist = dist(vertices[k], verticesSelected[l]);
				}
				if (minDist > biggestScore) {
					biggestScore = minDist;
					biggestIndex = k;
				}
			}
			verticesSelected[j] = vertices[biggestIndex];
		}
		double maxDist = 0;
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			distances[i][j][j] = 0;
			for (int k = j + 1; k < SAMPLE_POINTS; k++) {
				distances[i][j][k] = dist(verticesSelected[j], verticesSelected[k]);
				distances[i][k][j] = distances[i][j][k];
				if (distances[i][j][k] > maxDist)
					maxDist = distances[i][j][k];
			}
		}
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			for (int k = 0; k < SAMPLE_POINTS; k++) {
				distances[i][j][k] = exp(-distances[i][j][k] * distances[i][j][k] / (2 * maxDist * maxDist));
			}
		}
		delete[] verticesSelected;
	}

	VectorXd eigenValues[NUM_OF_FILES + 1];
	for (int i = 0; i < NUM_OF_FILES + 1; i++) {
		MatrixXd A(SAMPLE_POINTS, SAMPLE_POINTS);
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			for (int k = 0; k < SAMPLE_POINTS; k++) {
				A(j, k) = distances[i][j][k];
			}
		}
		MatrixXd temp;
		ComputeEigenDecomposition(A, eigenValues[i], temp);
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			if (eigenValues[i](j) < 0)
				eigenValues[i](j) = -eigenValues[i](j);
		}
	}
	double errors[NUM_OF_FILES];
	for (int i = 0; i < NUM_OF_FILES; i++) {
		errors[i] = CalculateError(eigenValues[NUM_OF_FILES], eigenValues[i]);
	}
	printResult(errors);

	for (int i = 0; i < NUM_OF_FILES + 1; i++) {
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			delete[] distances[i][j];
		}
		delete[] distances[i];
	}
	delete[] distances;
}

void Mesh3DScene::Task5C() {
	files = originalFiles;
	const string objDir = getBasePath() + "resources/obj/";
	double*** distances = new double** [NUM_OF_FILES + 1];
	for (int i = 0; i < NUM_OF_FILES + 1; i++) {
		distances[i] = new double* [SAMPLE_POINTS];
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			distances[i][j] = new double[SAMPLE_POINTS];
		}
	}
	int start = 0;
	if (!COMPUTE_GEODESIC)
		start = NUM_OF_FILES;
	for (int i = start; i < NUM_OF_FILES + 1; i++) {
		string objFile;
		if (i < NUM_OF_FILES)
			objFile = objDir + files[i] + ".obj";
		else
			objFile = objDir + FILENAME + ".obj";
		m_model_other = Mesh(objFile);
		vector<vec>& vertices = m_model_other.getVertices();
		int* verticesSelected = new int[SAMPLE_POINTS];
		verticesSelected[0] = rand() % vertices.size();

		MatrixXd A = MatrixXd::Zero(vertices.size(), vertices.size());
		for (int l = 0; l < m_model_other.getTriangles().size(); l++) {
			for (int j = 0; j < 2; j++) {
				for (int k = j + 1; k < 3; k++) {
					if (A(m_model_other.getTriangles()[l].v[j], m_model_other.getTriangles()[l].v[k]) == 0) {
						A(m_model_other.getTriangles()[l].v[j], m_model_other.getTriangles()[l].v[k]) = 1;
						A(m_model_other.getTriangles()[l].v[k], m_model_other.getTriangles()[l].v[j]) = 1;
					}
				}
			}
		}
		for (int j = 0; j < vertices.size(); j++) {
			for (int k = j + 1; k < vertices.size(); k++) {
				if (dist(vertices[j], vertices[k]) == 0) {
					A(j, k) = 1;
					A(k, j) = 1;
				}
			}
		}

		MatrixXd distance = MatrixXd::Zero(vertices.size(), SAMPLE_POINTS);

		for (int j = 1; j < SAMPLE_POINTS; j++) {
			int biggestIndex = geodesicDist(vertices, verticesSelected, j, A, distance);
			verticesSelected[j] = biggestIndex;
		}
		double maxDist = 0;
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			distances[i][j][j] = 0;
			for (int k = j + 1; k < SAMPLE_POINTS; k++) {
				distances[i][j][k] = distance(verticesSelected[k], j);		//geodesic(vertices, m_model_other.getTriangles(), verticesSelected[j], verticesSelected[k], A);
				distances[i][k][j] = distances[i][j][k];
				if (distances[i][j][k] > maxDist)
					maxDist = distances[i][j][k];
			}
		}
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			for (int k = 0; k < SAMPLE_POINTS; k++) {
				distances[i][j][k] = exp(-distances[i][j][k] * distances[i][j][k] / (2 * maxDist * maxDist));
			}
		}
		delete[] verticesSelected;
	}
	VectorXd eigenValues[NUM_OF_FILES + 1];
	for (int i = start; i < NUM_OF_FILES + 1; i++) {
		MatrixXd A(SAMPLE_POINTS, SAMPLE_POINTS);
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			for (int k = 0; k < SAMPLE_POINTS; k++) {
				A(j, k) = distances[i][j][k];
			}
		}
		MatrixXd temp;
		ComputeEigenDecomposition(A, eigenValues[i], temp);
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			if (eigenValues[i](j) < 0)
				eigenValues[i](j) = -eigenValues[i](j);
		}
	}
	if (!COMPUTE_GEODESIC)
		ReadGeodesicFromFile(eigenValues);
	else 
		SaveGeodesicToFile(eigenValues);
	double errors[NUM_OF_FILES];
	for (int i = 0; i < NUM_OF_FILES; i++) {
		errors[i] = CalculateError(eigenValues[NUM_OF_FILES], eigenValues[i]);
	}
	printResult(errors);

	for (int i = 0; i < NUM_OF_FILES + 1; i++) {
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			delete[] distances[i][j];
		}
		delete[] distances[i];
	}
	delete[] distances;
}

void Mesh3DScene::GetDifCoordsInNormalDirection(vector<vvr::Triangle>& m_triangles, vector<vec>& m_vertices, MatrixXd& DifCoords, MatrixXd& newDifCoords){
	int verticesCount = m_vertices.size();
	int trianglesCount = m_triangles.size();
	for (int i = 0; i < verticesCount; i++) {
		vec normal = vec(0, 0, 0);
		for (int j = 0; j < trianglesCount; j++) {
			if (m_triangles[j].vi1 == i || m_triangles[j].vi2 == i || m_triangles[j].vi3 == i) {
				double a = (m_vertices[m_triangles[j].vi1] - m_vertices[m_triangles[j].vi2]).Length();
				double b = (m_vertices[m_triangles[j].vi1] - m_vertices[m_triangles[j].vi3]).Length();
				double c = (m_vertices[m_triangles[j].vi2] - m_vertices[m_triangles[j].vi3]).Length();
				double p = (a + b + c) / 2;
				double area = sqrt(p * (p - a) * (p - b) * (p - c));
				normal += m_triangles[j].getNormal() * area;
			}
		}
		normal = normal.Normalized();
		vec oldDifCoord = vec(DifCoords(i, 0), DifCoords(i, 1), DifCoords(i, 2));
		vec newDifCoord = (oldDifCoord.Dot(normal) * (normal / (normal.x * normal.x + normal.y * normal.y + normal.z * normal.z)));
		newDifCoords(i, 0) = newDifCoord.x;
		newDifCoords(i, 1) = newDifCoord.y;
		newDifCoords(i, 2) = newDifCoord.z;
	}
}

void Mesh3DScene::ComputeEigenDecomposition(MatrixXd& Ls, VectorXd& eigenValues, MatrixXd& eigenVectors) {
	SelfAdjointEigenSolver<MatrixXd> solver;
	solver.compute(Ls); 
	eigenValues = solver.eigenvalues().real();
	eigenVectors = solver.eigenvectors().real();
}

void Mesh3DScene::SaveEigenToFile(const string eigenFile, VectorXd& eigenValues, MatrixXd& eigenVectors, int verticesCount) {
	ofstream myfile;
	myfile.open(eigenFile);
	myfile << verticesCount << "\n\n";
	myfile << eigenValues << "\n\n";
	myfile << eigenVectors;
	myfile.close();
}

void Mesh3DScene::ReadEigenFromFile(const string eigenFile, VectorXd& eigenValues, MatrixXd& eigenVectors, int verticesCount) {
	ifstream myfile;
	myfile.open(eigenFile);
	double count;
	myfile >> count;
	for (int i = 0; i < verticesCount; i++) {
		myfile >> eigenValues(i);
	}
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			myfile >> eigenVectors(i, j);
		}
	}
	myfile.close();
}

void Mesh3DScene::SaveCoordsToFile(const string coordsFile, MatrixXd& CoordsNewA, MatrixXd& DifCoordsNewA, MatrixXd& CoordsNewB, MatrixXd& DifCoordsNewB, MatrixXd& CoordsNewC, MatrixXd& DifCoordsNewC, MatrixXd& CoordsNewD, MatrixXd& DifCoordsNewD) {
	ofstream myfile;
	myfile.open(coordsFile);
	myfile << CoordsNewA << "\n\n";
	myfile << DifCoordsNewA.rows() << "\n";
	myfile << DifCoordsNewA << "\n\n";
	myfile << CoordsNewB << "\n\n";
	myfile << DifCoordsNewB.rows() << "\n";
	myfile << DifCoordsNewB << "\n\n";
	myfile << CoordsNewC << "\n\n";
	myfile << DifCoordsNewC.rows() << "\n";
	myfile << DifCoordsNewC << "\n\n";
	myfile << CoordsNewD << "\n\n";
	myfile << DifCoordsNewD.rows() << "\n";
	myfile << DifCoordsNewD << "\n\n";
	myfile.close();
}

void Mesh3DScene::ReadCoordsFromFile(const string coordsFile, MatrixXd& CoordsNewA, MatrixXd& DifCoordsNewA, MatrixXd& CoordsNewB, MatrixXd& DifCoordsNewB, MatrixXd& CoordsNewC, MatrixXd& DifCoordsNewC, MatrixXd& CoordsNewD, MatrixXd& DifCoordsNewD, int verticesCount) {
	int count;
	ifstream myfile;
	myfile.open(coordsFile);
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < 3; j++) {
			myfile >> CoordsNewA(i, j);
		}
	}
	myfile >> count;
	DifCoordsNewA.conservativeResize(count, NoChange);
	for (int i = 0; i < count; i++) {
		for (int j = 0; j < 3; j++) {
			myfile >> DifCoordsNewA(i, j);
		}
	}
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < 3; j++) {
			myfile >> CoordsNewB(i, j);
		}
	}
	myfile >> count;
	DifCoordsNewB.conservativeResize(count, NoChange);
	for (int i = 0; i < count; i++) {
		for (int j = 0; j < 3; j++) {
			myfile >> DifCoordsNewB(i, j);
		}
	}
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < 3; j++) {
			myfile >> CoordsNewC(i, j);
		}
	}
	myfile >> count;
	DifCoordsNewC.conservativeResize(count, NoChange);
	for (int i = 0; i < count; i++) {
		for (int j = 0; j < 3; j++) {
			myfile >> DifCoordsNewC(i, j);
		}
	}
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < 3; j++) {
			myfile >> CoordsNewD(i, j);
		}
	}
	myfile >> count;
	DifCoordsNewD.conservativeResize(count, NoChange);
	for (int i = 0; i < count; i++) {
		for (int j = 0; j < 3; j++) {
			myfile >> DifCoordsNewD(i, j);
		}
	}
	myfile.close();
	CoordsToModel(m_model_newA.getVertices(), CoordsNewA, verticesCount);
	CoordsToModel(m_model_newB.getVertices(), CoordsNewB, verticesCount);
	CoordsToModel(m_model_newC.getVertices(), CoordsNewC, verticesCount);
	CoordsToModel(m_model_newD.getVertices(), CoordsNewD, verticesCount);
}

void Mesh3DScene::SaveGeodesicToFile(VectorXd eigenValues[]) {
	const string geodesicFile = getBasePath() + "resources/geodesic/file.txt";
	ofstream myfile;
	myfile.open(geodesicFile);
	for (int i = 0; i < NUM_OF_FILES; i++) {
		myfile << files[i] << "\n\n";
		myfile << eigenValues[i] << "\n\n";
	}
	myfile.close();
}

void Mesh3DScene::ReadGeodesicFromFile(VectorXd eigenValues[]) {
	const string geodesicFile = getBasePath() + "resources/geodesic/file.txt";
	string name;
	ifstream myfile;
	myfile.open(geodesicFile);
	for (int i = 0; i < NUM_OF_FILES; i++) {
		myfile >> name;
		eigenValues[i] = VectorXd(SAMPLE_POINTS);
		for (int j = 0; j < SAMPLE_POINTS; j++) {
			myfile >> eigenValues[i](j);
		}
	}
	myfile.close();
}

void Mesh3DScene::CoordsToModel(std::vector<vec>& m_vertices_new, MatrixXd& CoordsNew, int verticesCount) {
	for (int i = 0; i < verticesCount; i++) {
		m_vertices_new[i].x = CoordsNew(i, 0);
		m_vertices_new[i].y = CoordsNew(i, 1);
		m_vertices_new[i].z = CoordsNew(i, 2);
	}
}

double Mesh3DScene::CalculateError(VectorXd& l1, VectorXd& l2) {
	double error = 0;
	for (int i = 0; i < SAMPLE_POINTS; i++) {
		error += (((sqrt(l1(i)) - sqrt(l2(i))) * (sqrt(l1(i)) - sqrt(l2(i)))) / (sqrt(l1(i)) + sqrt(l2(i))));
	}
	return 0.5 * error;
}

void Mesh3DScene::printResult(double errors[]) {
	for (int i = 0; i < NUM_OF_FILES - 1; i++) {
		for (int j = 0; j < NUM_OF_FILES - 1 - i; j++) {
			if (errors[j] > errors[j + 1]) {
				double temp = errors[j];
				errors[j] = errors[j + 1];
				errors[j + 1] = temp;
				string tempfile = files[j];
				files[j] = files[j + 1];
				files[j + 1] = tempfile;
			}
		}
	}
	cout << "\nGiven the model \"" << FILENAME << "\" the most similar shapes and their errors in decreasing order are:\n" << endl;
	for (int i = 0; i < NUM_OF_FILES; i++) {
		cout << files[i] << "  " << errors[i] << endl;
	}
}

void Mesh3DScene::arrowEvent(ArrowDir dir, int modif)
{
	
}

void Mesh3DScene::keyEvent(unsigned char key, bool up, int modif)
{
	Scene::keyEvent(key, up, modif);
	key = tolower(key);

	switch (key)
	{
	case 's': m_style_flag ^= FLAG_SHOW_SOLID; break;
	case 'w': m_style_flag ^= FLAG_SHOW_WIRE; break;
	case 'n': m_style_flag ^= FLAG_SHOW_NORMALS; break;
	case 'a': m_style_flag ^= FLAG_SHOW_AXES; break;
	case 'b': m_style_flag ^= FLAG_SHOW_AABB; break;
	case 'd': 
		m_style_flag ^= FLAG_SHOW_DIFCOORDS;
		if (m_style_flag & FLAG_SHOW_DIFCOORDS)
			cout << "showing error in deferential coordinates" << endl;
		else
			cout << "showing error in cartesian coordinates" << endl;
		break;
	case 'o': 
		m_style_flag ^= FLAG_SHOW_ORIGINAL_MODEL;
		if (m_style_flag & FLAG_SHOW_ORIGINAL_MODEL)
			cout << "showing original model" << endl;
		else
			cout << "showing reconstructed model using eigenvectors" << endl;
		break;
	case 'p': m_style_flag ^= FLAG_SHOW_POINTS; break;
	case 't': m_style_flag ^= FLAG_SHOW_TRIANGLES; break;
	case 'v': 
		ratio_flag = (ratio_flag + 1) % 4;
		if (ratio_flag == 0) {
			cout << "showing model reconstructed using 30% of eigenvectors" << endl;
			m_model_new_draw = &m_model_newA;
			m_points_coords3D_draw = &m_points_coords3DA;
			m_points_difCoords3D_draw = &m_points_difCoords3DA;
			m_triangles_coords3D_draw = &m_triangles_coords3DA;
			m_triangles_difCoords3D_draw = &m_triangles_difCoords3DA;
		}
		else if (ratio_flag == 1) {
			cout << "showing model reconstructed using 10% of eigenvectors" << endl;
			m_model_new_draw = &m_model_newB;
			m_points_coords3D_draw = &m_points_coords3DB;
			m_points_difCoords3D_draw = &m_points_difCoords3DB;
			m_triangles_coords3D_draw = &m_triangles_coords3DB;
			m_triangles_difCoords3D_draw = &m_triangles_difCoords3DB;
		}
		else if (ratio_flag == 2) {
			cout << "showing model reconstructed using 5% of eigenvectors" << endl;
			m_model_new_draw = &m_model_newC;
			m_points_coords3D_draw = &m_points_coords3DC;
			m_points_difCoords3D_draw = &m_points_difCoords3DC;
			m_triangles_coords3D_draw = &m_triangles_coords3DC;
			m_triangles_difCoords3D_draw = &m_triangles_difCoords3DC;
		}
		else if (ratio_flag == 3) {
			cout << "showing model reconstructed using 1% of eigenvectors" << endl;
			m_model_new_draw = &m_model_newD;
			m_points_coords3D_draw = &m_points_coords3DD;
			m_points_difCoords3D_draw = &m_points_difCoords3DD;
			m_triangles_coords3D_draw = &m_triangles_coords3DD;
			m_triangles_difCoords3D_draw = &m_triangles_difCoords3DD;
		}
		break;
	}
}

void Mesh3DScene::draw()
{
	if (m_style_flag & FLAG_SHOW_ORIGINAL_MODEL) {
		if (m_style_flag & FLAG_SHOW_SOLID) m_model.draw(m_obj_col, SOLID);
		if (m_style_flag & FLAG_SHOW_WIRE) m_model.draw(Colour::black, WIRE);
		if (m_style_flag & FLAG_SHOW_NORMALS) m_model.draw(Colour::black, NORMALS);
		if (m_style_flag & FLAG_SHOW_AXES) m_model.draw(Colour::black, AXES);
		if (m_style_flag & FLAG_SHOW_POINTS) {
			for (int i = 0; i < m_points3D.size(); i++) {
				m_points3D[i].draw();
			}
		}
		if (m_style_flag & FLAG_SHOW_TRIANGLES) {
			for (int i = 0; i < m_triangles3D.size(); i++) {
				m_triangles3D[i].draw();
			}
		}
	}
	else {
		if (m_style_flag & FLAG_SHOW_SOLID) m_model_new_draw->draw(m_obj_col, SOLID);
		if (m_style_flag & FLAG_SHOW_WIRE) m_model_new_draw->draw(Colour::black, WIRE);
		if (m_style_flag & FLAG_SHOW_NORMALS) m_model_new_draw->draw(Colour::black, NORMALS);
		if (m_style_flag & FLAG_SHOW_AXES) m_model_new_draw->draw(Colour::black, AXES);
		if (m_style_flag & FLAG_SHOW_DIFCOORDS) {
			if (m_style_flag & FLAG_SHOW_POINTS) {
				for (int i = 0; i < m_points_difCoords3D_draw->size(); i++) {
					(*m_points_difCoords3D_draw)[i].draw();
				}
			}
			if (m_style_flag & FLAG_SHOW_TRIANGLES) {
				for (int i = 0; i < m_triangles_difCoords3D_draw->size(); i++) {
					(*m_triangles_difCoords3D_draw)[i].draw();
				}
			}
		}
		else {
			if (m_style_flag & FLAG_SHOW_POINTS) {
				for (int i = 0; i < m_points_coords3D_draw->size(); i++) {
					(*m_points_coords3D_draw)[i].draw();
				}
			}
			if (m_style_flag & FLAG_SHOW_TRIANGLES) {
				for (int i = 0; i < m_triangles_coords3D_draw->size(); i++) {
					(*m_triangles_coords3D_draw)[i].draw();
				}
			}
		}
	}
}
 
int main(int argc, char* argv[])
{
	try {
		return vvr::mainLoop(argc, argv, new Mesh3DScene);
	}
	catch (std::string exc) {
		cerr << exc << endl;
		return 1;
	}
	catch (...)
	{
		cerr << "Unknown exception" << endl;
		return 1;
	}
}

void SparseIdentity(SparseMatrix<double>& I, int n) {
	I.reserve(VectorXi::Constant(n, 1));
	for (int i = 0; i < n; i++) {
		I.insert(i, i) = 1;
	}
}

void SparseDiagonalInverse(SparseMatrix<double>& D, SparseMatrix<double>& D_inverse, int n) {
	for (int i = 0; i < n; i++) {
		D_inverse.coeffRef(i, i) = 1 / D.coeffRef(i, i);
	}
}

double FindMax(MatrixXd& M, int n, int index) {
	int max = 0;
	for (int i = 0; i < n; i++) {
		if (M(i, index) > max)
			max = M(i, index);
	}
	return max;
}

void shuffle_arr(int* arr, size_t n) {
	if (n > 1) {
		size_t i;
		srand(time(NULL));
		for (i = 0; i < n - 1; i++) {
			size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
			int t = arr[j];
			arr[j] = arr[i];
			arr[i] = t;
		}
	}
}

double dist(vec v1, vec v2) {
	return sqrt((v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y) + (v1.z - v2.z) * (v1.z - v2.z));
}

int geodesicDist(vector<vec>& vertices, int* verticesSelected, int count, MatrixXd& A, MatrixXd& distance) {
	VectorXd temp = distance.col(count - 1);
	computeDist(vertices, verticesSelected[count - 1], temp, A);
	distance.col(count - 1) = temp;
	double maxDist = -1;
	int index = -1;
	for (int i = 0; i < vertices.size(); i++) {
		double minDist = 1e5;
		for (int j = 0; j < count; j++) {
			if (distance(i, j) < minDist)
				minDist = distance(i, j);
		}
		if (minDist > maxDist) {
			maxDist = minDist;
			index = i;
		}
	}
	return index;
}

void computeDist(vector<vec>& vertices, int index1, VectorXd& distances, MatrixXd& A) {
	for (int i = 0; i < vertices.size(); i++) {
		distances(i) = 1e5;
	}
	distances(index1) = 0;
	VectorXd visited = VectorXd::Zero(vertices.size());
	int toBeVisited = vertices.size();
	while (toBeVisited > 0) {
		double smallest = 1e6;
		int next = -1;
		for (int i = 0; i < vertices.size(); i++) {
			if (distances(i) < smallest && visited(i) == 0) {
				smallest = distances(i);
				next = i;
			}
		}
		toBeVisited--;
		visited(next) = 1;
		for (int i = 0; i < vertices.size(); i++) {
			if (A(i, next) == 1) {
				if (visited(i) == 0 && (distances(next) + dist(vertices[next], vertices[i]) < distances(i))) {
					distances(i) = distances(next) + dist(vertices[next], vertices[i]);
				}
			}
		}
	}
}