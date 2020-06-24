#include "SceneMesh3D.h"
#include "Eigen/Dense"
#include "iostream"
#include "fstream"
#include "Eigen/Sparse"
#include "Eigen/Eigenvalues"
#include "Eigen/SparseCholesky"
#include "ctime"

#define RATIO 0.1
#define FILENAME "dolphin"

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
	m_model_new = vvr::Mesh(m_model_original);
	reset();
}

void Mesh3DScene::reset()
{
	Scene::reset();

	//! Define what will be vissible by default
	m_style_flag = 0;
	m_style_flag |= FLAG_SHOW_SOLID;
	m_style_flag |= FLAG_SHOW_WIRE;
	m_style_flag |= FLAG_SHOW_AXES;
	m_style_flag |= FLAG_SHOW_AABB;
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
		m_model_new.setBigSize(getSceneWidth() / 2);
		m_model_new.update();
		m_model = m_model_original;
		Tasks();
		first_pass = false;
	}
}

void Mesh3DScene::Tasks()
{
	vector<vvr::Triangle>& m_triangles = m_model.getTriangles();
	vector<vec>& m_vertices = m_model.getVertices();
	//
	int verticesCount = m_vertices.size();
	int trianglesCount = m_triangles.size();
	cout << "verticesCount = " << verticesCount << endl;
	//
	for (int i = 0; i < verticesCount; i++) {
		//cout << m_vertices[i] << endl;
	}
	//Task1
	/*/
	MatrixXd A = MatrixXd::Zero(verticesCount, verticesCount);
	MatrixXd I = MatrixXd::Identity(verticesCount, verticesCount);
	MatrixXd D = MatrixXd::Zero(verticesCount, verticesCount);
	//
	for (int i = 0; i < trianglesCount; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = j + 1; k < 3; k++) {
				if (A(m_triangles[i].v[j], m_triangles[i].v[k]) == 0) {
					A(m_triangles[i].v[j], m_triangles[i].v[k]) = 1;
					A(m_triangles[i].v[k], m_triangles[i].v[j]) = 1;
					D(m_triangles[i].v[j], m_triangles[i].v[j])++;
					D(m_triangles[i].v[k], m_triangles[i].v[k])++;
				}
			}
		}
	}
	MatrixXd D_inverse = MatrixXd::Zero(verticesCount, verticesCount);
	for (int i = 0; i < verticesCount; i++) {
		D_inverse(i, i) = 1 / D(i, i);
	}
	MatrixXd L(verticesCount, verticesCount);
	//*/
	//*
	SparseMatrix<double> A(verticesCount, verticesCount);
	SparseMatrix<double> I(verticesCount, verticesCount);
	SparseMatrix<double> D(verticesCount, verticesCount);
	//give values to I
	I.reserve(VectorXi::Constant(verticesCount, 1));
	for (int i = 0; i < verticesCount; i++) {
		I.insert(i, i) = 1;
	}
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
	//invert D
	SparseMatrix<double> D_inverse(verticesCount, verticesCount);
	for (int i = 0; i < verticesCount; i++) {
		D_inverse.coeffRef(i, i) = 1 / D.coeffRef(i, i);
	}
	SparseMatrix<double> L(verticesCount, verticesCount);
	//*/
	L = I - D_inverse * A;
	//MatrixXd Ls = D - A;

	/*
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			cout << L.coeffRef(i, j) << "  ";
		}
		cout << endl;
	}
	*/

	MatrixXd Coords(verticesCount, 3);
	MatrixXd DifCoords(verticesCount, 3);
	for (int i = 0; i < verticesCount; i++) {
		Coords(i, 0) = (double)m_vertices[i].x;
		Coords(i, 1) = (double)m_vertices[i].y;
		Coords(i, 2) = (double)m_vertices[i].z;
	}

	cout << "print Coords" << endl;
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < 3; j++) {
			//cout << Coords(i, j) << "  ";
		}
		//cout << endl;
	}

	DifCoords = L * Coords;

	//find max of dif coords on each axis
	double maxX = 0;
	double maxY = 0;
	double maxZ = 0;
	for (int i = 0; i < verticesCount; i++) {
		if (DifCoords(i, 0) > maxX)
			maxX = DifCoords(i, 0);
		if (DifCoords(i, 1) > maxY)
			maxY = DifCoords(i, 1);
		if (DifCoords(i, 2) > maxZ)
			maxZ = DifCoords(i, 2);
	}
	//
	//*
	MatrixXd newDifCoords(verticesCount, 3);
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
		//cout << normal.x << "   " << normal.y << "   " << normal.z << endl;
		vec oldDifCoord = vec(DifCoords(i, 0), DifCoords(i, 1), DifCoords(i, 2));
		vec newDifCoord = (oldDifCoord.Dot(normal) * (normal/(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z)));
		newDifCoords(i, 0) = newDifCoord.x;
		newDifCoords(i, 1) = newDifCoord.y;
		newDifCoords(i, 2) = newDifCoord.z;
		//cout << newDifCoords(i, 0) << "   " << newDifCoords(i, 1) << "   " << newDifCoords(i, 2) << endl;
	}
	//*/
	//compute the magnitude of the dif coords of all the vertices
	VectorXd DifCoordsMagnitude(verticesCount);
	double maxMag = 0;
	for (int i = 0; i < verticesCount; i++) {
		DifCoordsMagnitude(i) = sqrt(newDifCoords(i, 0) * newDifCoords(i, 0) + newDifCoords(i, 1) * newDifCoords(i, 1) + newDifCoords(i, 2) * newDifCoords(i, 2));
		if (DifCoordsMagnitude(i) > maxMag)
			maxMag = DifCoordsMagnitude(i);
		//cout << DifCoordsMagnitude(i) << endl;
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
		/*
		m_triangles3D[i].setColourPerVertex(vvr::Colour(sqrt(DifCoordsMagnitude(m_triangles[i].vi1) / maxMag) * 255, sqrt(1 - (DifCoordsMagnitude(m_triangles[i].vi1) / maxMag)) * 255, 0),
			vvr::Colour(sqrt(DifCoordsMagnitude(m_triangles[i].vi2) / maxMag) * 255, sqrt(1 - (DifCoordsMagnitude(m_triangles[i].vi2) / maxMag)) * 255, 0),
			vvr::Colour(sqrt(DifCoordsMagnitude(m_triangles[i].vi3) / maxMag) * 255, sqrt(1 - (DifCoordsMagnitude(m_triangles[i].vi3) / maxMag)) * 255, 0));
		//*/
		m_triangles3D[i].setColourPerVertex(vvr::Colour((DifCoordsMagnitude(m_triangles[i].vi1) / maxMag) * 255, (1 - (DifCoordsMagnitude(m_triangles[i].vi1) / maxMag)) * 255, 0),
			vvr::Colour((DifCoordsMagnitude(m_triangles[i].vi2) / maxMag) * 255, (1 - (DifCoordsMagnitude(m_triangles[i].vi2) / maxMag)) * 255, 0),
			vvr::Colour((DifCoordsMagnitude(m_triangles[i].vi3) / maxMag) * 255, (1 - (DifCoordsMagnitude(m_triangles[i].vi3) / maxMag)) * 255, 0));
		//m_triangles3D[i].setColourPerVertex(vvr::Colour(abs(newDifCoords(m_triangles[i].vi1, 0) / maxX) * 255, abs(newDifCoords(m_triangles[i].vi1, 1) / maxY) * 255, abs(newDifCoords(m_triangles[i].vi1, 2) / maxZ) * 255),
		//vvr::Colour(abs(newDifCoords(m_triangles[i].vi2, 0) / maxX) * 255, abs(newDifCoords(m_triangles[i].vi2, 1) / maxY) * 255, abs(newDifCoords(m_triangles[i].vi2, 2) / maxZ) * 255),
		//vvr::Colour(abs(newDifCoords(m_triangles[i].vi3, 0) / maxX) * 255, abs(newDifCoords(m_triangles[i].vi3, 1) / maxY) * 255, abs(newDifCoords(m_triangles[i].vi3, 2) / maxZ) * 255));
	}
	//*/
	//
	//
	//
	//Task2
	//*
	cout << "Task1 done" << endl;
	
	MatrixXd Ls = D - A;
	/*
	for (int i = 0; i < verticesCount; i++) {
		//cout << L.coeffRef(0, i) << endl;
	}
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			cout << Ls(i, j)<<"  ";
		}
		cout << endl;
	}
	//*/
	//
	VectorXd eigenValues(verticesCount);
	MatrixXd eigenVectors(verticesCount, verticesCount);
	//
	/*
	//compute eigenvalues and eigenvectors
	//EigenSolver<MatrixXd> solver;
	SelfAdjointEigenSolver<MatrixXd> solver;
	solver.compute(Ls); //Ls
	//
	eigenValues = solver.eigenvalues().real();
	eigenVectors = solver.eigenvectors().real();
	//*/
	/*
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			if (eigenVectors(i, j) < 1.0e-9) {
				eigenVectors(i, j) = 0;
			}
		}
	}
	//*/
	//
	const string eigenDir = getBasePath() + "resources/eigen/";
	const string eigenFile = eigenDir + FILENAME + ".txt";
	//save eigenvalues and eigenvectors to file
	/*
	ofstream myfile;
	myfile.open(eigenFile);
	myfile << verticesCount << "\n\n";
	myfile << eigenValues << "\n\n";
	myfile << eigenVectors;
	myfile.close();
	//*/
	//
	//*
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
	//*/
	//
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
	//*
	MatrixXd Lamda = MatrixXd::Zero(verticesCount, verticesCount);
	for (int i = 0; i < verticesCount; i++) {
		Lamda(i, i) = eigenValues(i);
	}
	MatrixXd Q = eigenVectors;
	//*/
	/*
	for (int i = RATIO * verticesCount; i < verticesCount; i++) {
		Lamda(i, i) = 0;
		//cout << "cut" << endl;
		for (int j = 0; j < verticesCount; j++) {
			Q(j, i) = 0;
		}
	}
	//*/
	//*
	MatrixXd LsNew = (Q * Lamda) * Q.transpose(); //same as Q.inverse();
	for (int i = 0; i < verticesCount; i++) {
		//cout << LsNew(0, i) << endl;
	}

	MatrixXd temp = LsNew;
	//LsNew = D_inverse * temp;


	//*/
	/*
	cout << "\n\noriginal:" << endl;
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			//cout << LsNew(i, j) << "  ";
		}
		//cout << endl;
	}
	//*/
	//*
	MatrixXd LNew = D_inverse * LsNew;
	for (int i = 0; i < verticesCount; i++) {
		//cout << LNew(100, i) << endl;
	}

	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			//cout << LNew(i, j)<<"  ";
		}
		//cout << endl;
	}
	//*/
	/*
	//CompleteOrthogonalDecomposition<MatrixXd> co_decomp(LsNew);
	ColPivHouseholderQR<MatrixXd> col_decomp(LsNew);
	cout << "rank = " << col_decomp.rank() << endl;
	int rank;
	int dimCount = verticesCount;
	//
	int *options = new int[verticesCount];
	for (int i = 0; i < verticesCount; i++) {
		options[i] = i;
	}
	shuffle_arr(options, verticesCount);
	//
	while (true) {
		dimCount++;
		int count = dimCount - verticesCount - 1;
		LsNew.conservativeResize(dimCount, NoChange);
		for (int i = 0; i < verticesCount; i++) {
			LsNew(dimCount - 1, i) = 0;
		}
		LsNew(dimCount - 1, options[count]) = 1;
		//
		DifCoords.conservativeResize(dimCount, NoChange);
		for (int i = 0; i < 3; i++) {
			DifCoords(dimCount - 1, i) = Coords(options[count], i);
		}
		//
		//CompleteOrthogonalDecomposition<MatrixXd> co_decomp(LsNew);
		//*
		ColPivHouseholderQR<MatrixXd> col_decomp(LsNew);
		rank = col_decomp.rank();
		cout << "rank = " << rank << endl;
		if (rank >= verticesCount) {
			cout << "hehehe" << endl;
			break;
		}
		//*/
		//
		/*
		if (count >= (1 - RATIO) * verticesCount) {
			ColPivHouseholderQR<MatrixXd> col_decomp(LsNew);
			rank = col_decomp.rank();
			cout << "rank = " << rank << endl;
			//break;
		}
		//*/
	//}
	/*
	for (int i = 0; i < dimCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			//cout << LsNew(i, j) << "  ";
		}
		//cout << endl;
	}
	//*/

	//*
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
	MatrixXd CoordsNew(verticesCount, 3);
	//CoordsNew = ((LsNew.transpose() * LsNew).inverse()) * LsNew.transpose() * DifCoords;
	CoordsNew = Q * Xtilda;
	
	cout << "print new Coords" << endl;
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < 3; j++) {
			//cout << CoordsNew(i, j) << "  ";
		}
		//cout << endl;
	}
	
	vector<vec>& m_vertices_new = m_model_new.getVertices();
	for (int i = 0; i < verticesCount; i++) {
		m_vertices_new[i].x = CoordsNew(i, 0);
		m_vertices_new[i].y = CoordsNew(i, 1);
		m_vertices_new[i].z = CoordsNew(i, 2);
	}
	//*/


	cout << "done!!!" << endl;
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
	}
}

void Mesh3DScene::draw()
{
	/*
	if (m_style_flag & FLAG_SHOW_SOLID) m_model.draw(m_obj_col, SOLID);
	if (m_style_flag & FLAG_SHOW_WIRE) m_model.draw(Colour::black, WIRE);
	if (m_style_flag & FLAG_SHOW_NORMALS) m_model.draw(Colour::black, NORMALS);
	if (m_style_flag & FLAG_SHOW_AXES) m_model.draw(Colour::black, AXES);
	//*/
	//*
	if (m_style_flag & FLAG_SHOW_SOLID) m_model_new.draw(m_obj_col, SOLID);
	if (m_style_flag & FLAG_SHOW_WIRE) m_model_new.draw(Colour::black, WIRE);
	if (m_style_flag & FLAG_SHOW_NORMALS) m_model_new.draw(Colour::black, NORMALS);
	if (m_style_flag & FLAG_SHOW_AXES) m_model_new.draw(Colour::black, AXES);
	//*/

	//*
	for (int i = 0; i < m_points3D.size(); i++) {
		//m_points3D[i].draw();
	 }
	 //*/
	for (int i = 0; i < m_triangles3D.size(); i++) {
		//m_triangles3D[i].draw();
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

void Mesh3DScene::Task1(std::vector<vvr::Triangle>& m_triangles, std::vector<vec>& m_vertices){//, Eigen::MatrixXd& L, Eigen::MatrixXd& A, Eigen::MatrixXd& D){

}

