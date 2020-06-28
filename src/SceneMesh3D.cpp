#include "SceneMesh3D.h"
#include "iostream"
#include "fstream"
#include "ctime"

//#define RATIO 0.05
#define FILENAME "bone"

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
	cout << "Task1 done" << endl;

	//Task2
	MatrixXd Ls = D - A;
	VectorXd eigenValues(verticesCount);
	MatrixXd eigenVectors(verticesCount, verticesCount);
	/*
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			cout << Ls(i, j)<<"  ";
		}
		cout << endl;
	}
	//*/
	Task2(Ls, eigenValues, eigenVectors, verticesCount);
	/*
	for (int i = 0; i < verticesCount; i++) {
		cout << eigenValues(i) << endl;
	}
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			cout << eigenVectors(i, j) << "   ";
		}
		cout << endl;
	}
	//*/
	//Task3
	MatrixXd Q = eigenVectors;
	MatrixXd CoordsNewA(verticesCount, 3);
	MatrixXd DifCoordsNewA(verticesCount, 3);
	MatrixXd CoordsNewB(verticesCount, 3);
	MatrixXd DifCoordsNewB(verticesCount, 3);
	MatrixXd CoordsNewC(verticesCount, 3);
	MatrixXd DifCoordsNewC(verticesCount, 3);
	MatrixXd CoordsNewD(verticesCount, 3);
	MatrixXd DifCoordsNewD(verticesCount, 3);
	Ls = D_inverse * Ls;
	Task3(Q, Coords, DifCoords, CoordsNewA, DifCoordsNewA, CoordsNewB, DifCoordsNewB, CoordsNewC, DifCoordsNewC, CoordsNewD, DifCoordsNewD, Ls, verticesCount);
	cout << "new Coords ready" << endl;
	Task4(Coords, CoordsNewA, CoordsNewB, CoordsNewC, CoordsNewD, DifCoords, DifCoordsNewA, DifCoordsNewB, DifCoordsNewC, DifCoordsNewD, L);
	cout << "done!!!" << endl;
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
	SparseDiagonalInverse(D, D_inverse, verticesCount);
	L = I - D_inverse * A;
	/*
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			cout << L.coeffRef(i, j) << "  ";
		}
		cout << endl;
	}
	*/
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
}

void Mesh3DScene::Task2(MatrixXd& Ls, VectorXd& eigenValues, MatrixXd& eigenVectors, int verticesCount) {
	const string eigenDir = getBasePath() + "resources/eigen/";
	const string eigenFile = eigenDir + FILENAME + ".txt";
	//ComputeEigenDecomposition(Ls, eigenValues, eigenVectors);
	//SaveEigenToFile(eigenFile, eigenValues, eigenVectors, verticesCount);
	ReadEigenFromFile(eigenFile, eigenValues, eigenVectors, verticesCount);
	cout << "fileIO complete" << endl;
}

void Mesh3DScene::Task3(MatrixXd& Q, MatrixXd& Coords, MatrixXd& DifCoords, MatrixXd& CoordsNewA, MatrixXd& DifCoordsNewA, MatrixXd& CoordsNewB, MatrixXd& DifCoordsNewB, MatrixXd& CoordsNewC, MatrixXd& DifCoordsNewC, MatrixXd& CoordsNewD, MatrixXd& DifCoordsNewD, MatrixXd& Ls, int verticesCount) {
	double ratios[4] = { 0.3, 0.1, 0.05, 0.01 };
	for (int i = 0; i < 4; i++) {
		MatrixXd Qcopy = Q;
		MatrixXd LsCopy = Ls;
		if (i == 0)
			Task3Sub(Qcopy, Coords, DifCoords, CoordsNewA, DifCoordsNewA, LsCopy, m_model_newA.getVertices(), verticesCount, ratios[i]);
		if (i == 1)
			Task3Sub(Qcopy, Coords, DifCoords, CoordsNewB, DifCoordsNewB, LsCopy, m_model_newB.getVertices(), verticesCount, ratios[i]);
		if (i == 2)
			Task3Sub(Qcopy, Coords, DifCoords, CoordsNewC, DifCoordsNewC, LsCopy, m_model_newC.getVertices(), verticesCount, ratios[i]);
		if (i == 3)
			Task3Sub(Qcopy, Coords, DifCoords, CoordsNewD, DifCoordsNewD, LsCopy, m_model_newD.getVertices(), verticesCount, ratios[i]);
	}
}

void Mesh3DScene::Task3Sub(MatrixXd& Q, MatrixXd& Coords, MatrixXd& DifCoords, MatrixXd& CoordsNew, MatrixXd& DifCoordsNew, MatrixXd& Ls, vector<vec>& m_vertices_new, int verticesCount, double ratio) {
	/*
	MatrixXd Xtilda = Q.transpose() * Coords;
	for (int i = RATIO * verticesCount; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			Q(j, i) = 0;
		}
		for (int j = 0; j < 3; j++) {
			Xtilda(i, j) = 0;
		}
	}
	//*/
	//*
	MatrixXd Deltatilda = Q.transpose() * DifCoords;
	for (int i = ratio * verticesCount; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			Q(j, i) = 0;
		}
		for (int j = 0; j < 3; j++) {
			Deltatilda(i, j) = 0;
		}
	}
	DifCoordsNew = Q * Deltatilda;
	//*/

	cout << "rank calculation begins" << endl;
	int rank;
	int dimCount = verticesCount;
	while (true) {
		dimCount++;
		int count = dimCount - verticesCount - 1;
		Ls.conservativeResize(dimCount, NoChange);
		for (int i = 0; i < verticesCount; i++) {
			Ls(dimCount - 1, i) = 0;
		}
		Ls(dimCount - 1, count) = 1;

		DifCoordsNew.conservativeResize(dimCount, NoChange);
		for (int i = 0; i < 3; i++) {
			DifCoordsNew(dimCount - 1, i) = Coords(count, i);
		}

		//*
		ColPivHouseholderQR<MatrixXd> col_decomp(Ls);
		rank = col_decomp.rank();
		cout << "rank = " << rank << endl;
		if (rank >= verticesCount) {
			cout << "hehehe" << count << endl;
			break;
		}
		//*/
		/*
		if (count >= 1) {
			//ColPivHouseholderQR<MatrixXd> col_decomp(Ls);
			//rank = col_decomp.rank();
			//cout << "rank = " << rank << endl;
			break;
		}
		//*/
	}

	cout << "new Dif Coords ready" << endl;
	//CoordsNew = ((Ls.transpose() * Ls).inverse()) * Ls.transpose() * DifCoordsNew;
	CoordsNew = Ls.colPivHouseholderQr().solve(DifCoordsNew);
	//CoordsNew = Q * Xtilda;

	for (int i = 0; i < verticesCount; i++) {
		m_vertices_new[i].x = CoordsNew(i, 0);
		m_vertices_new[i].y = CoordsNew(i, 1);
		m_vertices_new[i].z = CoordsNew(i, 2);
	}
}

void Mesh3DScene::Task4(const MatrixXd& Coords, const MatrixXd& CoordsNewA, const MatrixXd& CoordsNewB, const MatrixXd& CoordsNewC, const MatrixXd& CoordsNewD, const MatrixXd& DifCoords, const MatrixXd& DifCoordsNewA, const MatrixXd& DifCoordsNewB, const MatrixXd& DifCoordsNewC, const MatrixXd& DifCoordsNewD, const Eigen::MatrixXd& L) {
	vector<vvr::Triangle>& m_triangles = m_model.getTriangles();
	Task4Sub(m_triangles, m_model_newA.getVertices(), m_points_coords3DA, m_triangles_coords3DA, m_points_difCoords3DA, m_triangles_difCoords3DA, Coords, CoordsNewA, DifCoords, DifCoordsNewA, L);
	Task4Sub(m_triangles, m_model_newB.getVertices(), m_points_coords3DB, m_triangles_coords3DB, m_points_difCoords3DB, m_triangles_difCoords3DB, Coords, CoordsNewB, DifCoords, DifCoordsNewB, L);
	Task4Sub(m_triangles, m_model_newC.getVertices(), m_points_coords3DC, m_triangles_coords3DC, m_points_difCoords3DC, m_triangles_difCoords3DC, Coords, CoordsNewC, DifCoords, DifCoordsNewC, L);
	Task4Sub(m_triangles, m_model_newD.getVertices(), m_points_coords3DD, m_triangles_coords3DD, m_points_difCoords3DD, m_triangles_difCoords3DD, Coords, CoordsNewD, DifCoords, DifCoordsNewD, L);
}

void Mesh3DScene::Task4Sub(vector<vvr::Triangle>& m_triangles, vector<vec>& m_vertices_new, vector<Point3D>& m_points_coords3D, vector<Triangle3D>& m_triangles_coords3D, vector<Point3D>& m_points_difCoords3D, vector<Triangle3D>& m_triangles_difCoords3D, const MatrixXd& Coords, const MatrixXd& CoordsNew, const MatrixXd& DifCoords, const MatrixXd& DifCoordsNew, const Eigen::MatrixXd& L) {
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
		//cout << Dif(i) << endl;
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
	//
	MatrixXd DifCoordsNewNew = L * CoordsNew;
	double totalMistake = 0;
	for (int i = 0; i < verticesCount; i++) {
		vec diff = vec(DifCoordsNewNew(i, 0) - DifCoordsNew(i, 0), DifCoordsNewNew(i, 1) - DifCoordsNew(i, 1), DifCoordsNewNew(i, 2) - DifCoordsNew(i, 2));
		//cout << diff.Length() << endl;
		totalMistake += diff.Length();
	}
	cout << "total mistake = " << totalMistake << endl;
	//*
	VectorXd DifInDifCoords(verticesCount);
	double maxDifInDifCoords = 0;
	for (int i = 0; i < verticesCount; i++) {
		vec diference = vec(DifCoordsNew(i, 0) - DifCoords(i, 0), DifCoordsNew(i, 1) - DifCoords(i, 1), DifCoordsNew(i, 2) - DifCoords(i, 2));
		DifInDifCoords(i) = diference.Length();
		if (DifInDifCoords(i) > maxDifInDifCoords) {
			maxDifInDifCoords = DifInDifCoords(i);
		}
	}
	/*
	cout << "\n\n\n";
	for (int i = 0; i < verticesCount; i++) {
		cout << DifInDifCoords(i) << endl;
	}
	//*/
	//
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
	//*/
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
	//EigenSolver<MatrixXd> solver;
	SelfAdjointEigenSolver<MatrixXd> solver;
	solver.compute(Ls); 
	eigenValues = solver.eigenvalues().real();
	eigenVectors = solver.eigenvectors().real();
	/*
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			if (eigenVectors(i, j) < 1.0e-9) {
				eigenVectors(i, j) = 0;
			}
		}
	}
	//*/
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
		//cout << "in" << endl;
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
		/*
		if (ratio_flag == 0) {
			if (m_style_flag & FLAG_SHOW_SOLID) m_model_newA.draw(m_obj_col, SOLID);
			if (m_style_flag & FLAG_SHOW_WIRE) m_model_newA.draw(Colour::black, WIRE);
			if (m_style_flag & FLAG_SHOW_NORMALS) m_model_newA.draw(Colour::black, NORMALS);
			if (m_style_flag & FLAG_SHOW_AXES) m_model_newA.draw(Colour::black, AXES);
			if (m_style_flag & FLAG_SHOW_DIFCOORDS) {
				if (m_style_flag & FLAG_SHOW_POINTS) {
					for (int i = 0; i < m_points_difCoords3DA.size(); i++) {
						m_points_difCoords3DA[i].draw();
					}
				}
				if (m_style_flag & FLAG_SHOW_TRIANGLES) {
					for (int i = 0; i < m_triangles_difCoords3DA.size(); i++) {
						m_triangles_difCoords3DA[i].draw();
					}
				}
			}
			else {
				if (m_style_flag & FLAG_SHOW_POINTS) {
					for (int i = 0; i < m_points_coords3DA.size(); i++) {
						m_points_coords3DA[i].draw();
					}
				}
				if (m_style_flag & FLAG_SHOW_TRIANGLES) {
					for (int i = 0; i < m_triangles_coords3DA.size(); i++) {
						m_triangles_coords3DA[i].draw();
					}
				}
			}
		}
		else if (ratio_flag == 1) {
			if (m_style_flag & FLAG_SHOW_SOLID) m_model_newB.draw(m_obj_col, SOLID);
			if (m_style_flag & FLAG_SHOW_WIRE) m_model_newB.draw(Colour::black, WIRE);
			if (m_style_flag & FLAG_SHOW_NORMALS) m_model_newB.draw(Colour::black, NORMALS);
			if (m_style_flag & FLAG_SHOW_AXES) m_model_newB.draw(Colour::black, AXES);
			if (m_style_flag & FLAG_SHOW_DIFCOORDS) {
				if (m_style_flag & FLAG_SHOW_POINTS) {
					for (int i = 0; i < m_points_difCoords3DB.size(); i++) {
						m_points_difCoords3DB[i].draw();
					}
				}
				if (m_style_flag & FLAG_SHOW_TRIANGLES) {
					for (int i = 0; i < m_triangles_difCoords3DB.size(); i++) {
						m_triangles_difCoords3DB[i].draw();
					}
				}
			}
			else {
				if (m_style_flag & FLAG_SHOW_POINTS) {
					for (int i = 0; i < m_points_coords3DB.size(); i++) {
						m_points_coords3DB[i].draw();
					}
				}
				if (m_style_flag & FLAG_SHOW_TRIANGLES) {
					for (int i = 0; i < m_triangles_coords3DB.size(); i++) {
						m_triangles_coords3DB[i].draw();
					}
				}
			}
		}
		else if (ratio_flag == 2) {
			if (m_style_flag & FLAG_SHOW_SOLID) m_model_newC.draw(m_obj_col, SOLID);
			if (m_style_flag & FLAG_SHOW_WIRE) m_model_newC.draw(Colour::black, WIRE);
			if (m_style_flag & FLAG_SHOW_NORMALS) m_model_newC.draw(Colour::black, NORMALS);
			if (m_style_flag & FLAG_SHOW_AXES) m_model_newC.draw(Colour::black, AXES);
			if (m_style_flag & FLAG_SHOW_DIFCOORDS) {
				if (m_style_flag & FLAG_SHOW_POINTS) {
					for (int i = 0; i < m_points_difCoords3DC.size(); i++) {
						m_points_difCoords3DC[i].draw();
					}
				}
				if (m_style_flag & FLAG_SHOW_TRIANGLES) {
					for (int i = 0; i < m_triangles_difCoords3DC.size(); i++) {
						m_triangles_difCoords3DC[i].draw();
					}
				}
			}
			else {
				if (m_style_flag & FLAG_SHOW_POINTS) {
					for (int i = 0; i < m_points_coords3DC.size(); i++) {
						m_points_coords3DC[i].draw();
					}
				}
				if (m_style_flag & FLAG_SHOW_TRIANGLES) {
					for (int i = 0; i < m_triangles_coords3DC.size(); i++) {
						m_triangles_coords3DC[i].draw();
					}
				}
			}
		}
		else if (ratio_flag == 3) {
			if (m_style_flag & FLAG_SHOW_SOLID) m_model_newD.draw(m_obj_col, SOLID);
			if (m_style_flag & FLAG_SHOW_WIRE) m_model_newD.draw(Colour::black, WIRE);
			if (m_style_flag & FLAG_SHOW_NORMALS) m_model_newD.draw(Colour::black, NORMALS);
			if (m_style_flag & FLAG_SHOW_AXES) m_model_newD.draw(Colour::black, AXES);
			if (m_style_flag & FLAG_SHOW_DIFCOORDS) {
				if (m_style_flag & FLAG_SHOW_POINTS) {
					for (int i = 0; i < m_points_difCoords3DD.size(); i++) {
						m_points_difCoords3DD[i].draw();
					}
				}
				if (m_style_flag & FLAG_SHOW_TRIANGLES) {
					for (int i = 0; i < m_triangles_difCoords3DD.size(); i++) {
						m_triangles_difCoords3DD[i].draw();
					}
				}
			}
			else {
				if (m_style_flag & FLAG_SHOW_POINTS) {
					for (int i = 0; i < m_points_coords3DD.size(); i++) {
						m_points_coords3DD[i].draw();
					}
				}
				if (m_style_flag & FLAG_SHOW_TRIANGLES) {
					for (int i = 0; i < m_triangles_coords3DD.size(); i++) {
						m_triangles_coords3DD[i].draw();
					}
				}
			}
		}
		//*/
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