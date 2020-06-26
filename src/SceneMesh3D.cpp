#include "SceneMesh3D.h"
#include "iostream"
#include "fstream"
#include "ctime"

#define RATIO 0.1
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
	cout << "fileIO complete" << endl;
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
	MatrixXd Q = eigenVectors;

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
	for (int i = RATIO * verticesCount; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			Q(j, i) = 0;
		}
		for (int j = 0; j < 3; j++) {
			Deltatilda(i, j) = 0;
		}
	}
	MatrixXd DifCoordsNew(verticesCount, 3);
	DifCoordsNew = Q * Deltatilda;
	//*/

	Ls = D_inverse * Ls;
	//*
	cout << "rank calculation begins" << endl;
	//ColPivHouseholderQR<MatrixXd> col_decomp(Ls);
	//cout << "rank = " << col_decomp.rank() << endl;
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
		//
		DifCoordsNew.conservativeResize(dimCount, NoChange);
		for (int i = 0; i < 3; i++) {
			DifCoordsNew(dimCount - 1, i) = Coords(count, i);
		}
		//
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

	//*/
	cout << "new Dif Coords ready" << endl;
	//*
	MatrixXd CoordsNew(verticesCount, 3);
	//CoordsNew = ((Ls.transpose() * Ls).inverse()) * Ls.transpose() * DifCoordsNew;
	CoordsNew = Ls.colPivHouseholderQr().solve(DifCoordsNew);
	//CoordsNew = Q * Xtilda;
	
	vector<vec>& m_vertices_new = m_model_new.getVertices();
	for (int i = 0; i < verticesCount; i++) {
		m_vertices_new[i].x = CoordsNew(i, 0);
		m_vertices_new[i].y = CoordsNew(i, 1);
		m_vertices_new[i].z = CoordsNew(i, 2);
	}
	//*/
	cout << "new Coords ready" << endl;
	
	Task4(m_triangles, m_vertices_new, Coords, CoordsNew, DifCoords, DifCoordsNew, L);


	cout << "done!!!" << endl;
}

void Mesh3DScene::Task4(std::vector<vvr::Triangle>& m_triangles, std::vector<vec>& m_vertices_new, const MatrixXd& Coords, const MatrixXd& CoordsNew, const MatrixXd& DifCoords, const MatrixXd& DifCoordsNew, const Eigen::MatrixXd& L) {
	int verticesCount = m_vertices_new.size();
	int trianglesCount = m_triangles.size();
	//*
	VectorXd Dif(verticesCount);
	double maxDif = 0;
	for (int i = 0; i < verticesCount; i++) {
		vec diference = vec(CoordsNew(i, 0) - Coords(i, 0), CoordsNew(i, 1) - Coords(i, 1), CoordsNew(i, 2) - Coords(i, 2));
		Dif(i) = diference.Length();
		if (Dif(i) > maxDif) {
			maxDif = Dif(i);
		}
	}
	/*
	for (int i = 0; i < verticesCount; i++) {
		cout << Dif(i) << endl;
	}
	//*/
	//
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
	//*/
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
	case 'd': {
		m_style_flag ^= FLAG_SHOW_DIFCOORDS;
		if (m_style_flag & FLAG_SHOW_DIFCOORDS)
			cout << "showing error in deferential coordinates" << endl;
		else
			cout << "showing error in cartesian coordinates" << endl;
		break;
	}
	case 'o': {
		m_style_flag ^= FLAG_SHOW_ORIGINAL_MODEL;
		if (m_style_flag & FLAG_SHOW_ORIGINAL_MODEL)
			cout << "showing original model" << endl;
		else
			cout << "showing reconstructed model using eigenvectors" << endl;
		break;
	}
	case 'p': m_style_flag ^= FLAG_SHOW_POINTS; break;
	case 't': m_style_flag ^= FLAG_SHOW_TRIANGLES; break;
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
		if (m_style_flag & FLAG_SHOW_SOLID) m_model_new.draw(m_obj_col, SOLID);
		if (m_style_flag & FLAG_SHOW_WIRE) m_model_new.draw(Colour::black, WIRE);
		if (m_style_flag & FLAG_SHOW_NORMALS) m_model_new.draw(Colour::black, NORMALS);
		if (m_style_flag & FLAG_SHOW_AXES) m_model_new.draw(Colour::black, AXES);
		if (m_style_flag & FLAG_SHOW_DIFCOORDS) {
			if (m_style_flag & FLAG_SHOW_POINTS) {
				for (int i = 0; i < m_points3D.size(); i++) {
					m_points_difCoords3D[i].draw();
				}
			}
			if (m_style_flag & FLAG_SHOW_TRIANGLES) {
				for (int i = 0; i < m_triangles_difCoords3D.size(); i++) {
					m_triangles_difCoords3D[i].draw();
				}
			}
		}
		else {
			if (m_style_flag & FLAG_SHOW_POINTS) {
				for (int i = 0; i < m_points3D.size(); i++) {
					m_points_coords3D[i].draw();
				}
			}
			if (m_style_flag & FLAG_SHOW_TRIANGLES) {
				for (int i = 0; i < m_triangles_coords3D.size(); i++) {
					m_triangles_coords3D[i].draw();
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

void Mesh3DScene::Task1(std::vector<vvr::Triangle>& m_triangles, std::vector<vec>& m_vertices){

}

