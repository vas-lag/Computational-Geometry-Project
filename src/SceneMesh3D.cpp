#include "SceneMesh3D.h"
#include "Eigen/Dense"
#include "iostream"
#include "Eigen/Sparse"
#include "Eigen/Eigenvalues"
//#include "Eigen/LU"

#define RATIO 1.0

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
	const string objFile = objDir + "cube.obj";
	m_model_original = vvr::Mesh(objFile);
	reset();
}

void Mesh3DScene::reset()
{
	Scene::reset();

	//! Define plane
	m_plane_d = 0;
	m_plane = Plane(vec(0, 1, 1).Normalized(), m_plane_d);

	//! Define what will be vissible by default
	m_style_flag = 0;
	m_style_flag |= FLAG_SHOW_SOLID;
	m_style_flag |= FLAG_SHOW_WIRE;
	m_style_flag |= FLAG_SHOW_AXES;
	m_style_flag |= FLAG_SHOW_AABB;
	//m_style_flag |= FLAG_SHOW_PLANE;
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
	MatrixXd Coords(verticesCount, 3);
	MatrixXd DifCoords(verticesCount, 3);
	for (int i = 0; i < verticesCount; i++) {
		Coords(i, 0) = (double)m_vertices[i].x;
		Coords(i, 1) = (double)m_vertices[i].y;
		Coords(i, 2) = (double)m_vertices[i].z;
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
		DifCoords(i, 0) = newDifCoord.x;
		DifCoords(i, 1) = newDifCoord.y;
		DifCoords(i, 2) = newDifCoord.z;
		//cout << DifCoords(i, 0) << "   " << DifCoords(i, 1) << "   " << DifCoords(i, 2) << endl;
	}
	//*/
	//compute the magnitude of the dif coords of all the vertices
	VectorXd DifCoordsMagnitude(verticesCount);
	double maxMag = 0;
	for (int i = 0; i < verticesCount; i++) {
		DifCoordsMagnitude(i) = sqrt(DifCoords(i, 0) * DifCoords(i, 0) + DifCoords(i, 1) * DifCoords(i, 1) + DifCoords(i, 2) * DifCoords(i, 2));
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
		//m_triangles3D[i].setColourPerVertex(vvr::Colour(abs(DifCoords(m_triangles[i].vi1, 0) / maxX) * 255, abs(DifCoords(m_triangles[i].vi1, 1) / maxY) * 255, abs(DifCoords(m_triangles[i].vi1, 2) / maxZ) * 255),
		//vvr::Colour(abs(DifCoords(m_triangles[i].vi2, 0) / maxX) * 255, abs(DifCoords(m_triangles[i].vi2, 1) / maxY) * 255, abs(DifCoords(m_triangles[i].vi2, 2) / maxZ) * 255),
		//vvr::Colour(abs(DifCoords(m_triangles[i].vi3, 0) / maxX) * 255, abs(DifCoords(m_triangles[i].vi3, 1) / maxY) * 255, abs(DifCoords(m_triangles[i].vi3, 2) / maxZ) * 255));
	}
	//*/
	//
	//
	//
	//Task2
	//*
	cout << "Task1 done" << endl;
	MatrixXd Ls = D - A;
	for (int i = 0; i < verticesCount; i++) {
		//cout << L.coeffRef(0, i) << endl;
	}
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			//cout << Ls(i, j)<<"  ";
		}
		//cout << endl;
	}
	EigenSolver<MatrixXd> solver;
	solver.compute(Ls); //Ls
	//
	VectorXd eigenValues(verticesCount);
	MatrixXd eigenVectors(verticesCount, verticesCount);
	//
	eigenValues = solver.eigenvalues().real();
	eigenVectors = solver.eigenvectors().real();
	//
	//std::sort(eigenValues.data(), eigenValues.data() + eigenValues.size());
	cout << "sort" << endl;
	//*
	//sort eigenvalues in increasing order 
	//and their eigenvectors accordingly
	for (int i = 0; i < verticesCount; i++) {
		for (int j = i + 1; j < verticesCount; j++) {
			if (eigenValues(i) > eigenValues(j)) {
				double temp = eigenValues(i);
				eigenValues(i) = eigenValues(j);
				eigenValues(j) = temp;
				for (int k = 0; k < verticesCount; k++) {
					temp = eigenVectors(k, i);
					eigenVectors(k, i) = eigenVectors(k, j);
					eigenVectors(k, j) = temp;
				}
			}
		}
	}
	//
	for (int i = 0; i < verticesCount; i++) {
		//cout << eigenValues(i) << endl;
	}
	MatrixXd Lamda = MatrixXd::Zero(verticesCount, verticesCount);
	for (int i = 0; i < verticesCount; i++) {
		Lamda(i, i) = eigenValues(i);
	}
	MatrixXd Q = eigenVectors;
	//*
	for (int i = RATIO * verticesCount; i < verticesCount; i++) {
		Lamda(i, i) = 0;
		cout << "cut" << endl;
		for (int j = 0; j < verticesCount; j++) {
			Q(j, i) = 0;
		}
	}
	//*/
	MatrixXd LsNew = (Q * Lamda) * Q.inverse(); //same as eigenVectors.inverse();
	for (int i = 0; i < verticesCount; i++) {
		//cout << LsNew(0, i) << endl;
	}
	//*/
	
	//*
	MatrixXd LNew = D_inverse * LsNew;
	for (int i = 0; i < verticesCount; i++) {
		//cout << LNew(100, i) << endl;
	}

	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < verticesCount; j++) {
			cout << LNew(i, j)<<"  ";
		}
		cout << endl;
	}

	FullPivLU<MatrixXd> lu_decomp(LNew);
	cout << "rank = " << lu_decomp.rank() << endl;

	cout << "done!!!" << endl;
	//*/
}

void Mesh3DScene::arrowEvent(ArrowDir dir, int modif)
{
	math::vec n = m_plane.normal;
	if (dir == UP) m_plane_d += 1;
	if (dir == DOWN) m_plane_d -= 1;
	else if (dir == LEFT) n = math::float3x3::RotateY(DegToRad(1)).Transform(n);
	else if (dir == RIGHT) n = math::float3x3::RotateY(DegToRad(-1)).Transform(n);
	m_plane = Plane(n.Normalized(), m_plane_d);
 
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
	case 'p': m_style_flag ^= FLAG_SHOW_PLANE; break;
	case 'b': m_style_flag ^= FLAG_SHOW_AABB; break;
	}
}

void Mesh3DScene::draw()
{
	//! Draw plane
	if (m_style_flag & FLAG_SHOW_PLANE) {
		vvr::Colour colPlane(0x41, 0x14, 0xB3);
		float u = 20, v = 20;
		math::vec p0(m_plane.Point(-u, -v, math::vec(0, 0, 0)));
		math::vec p1(m_plane.Point(-u, v, math::vec(0, 0, 0)));
		math::vec p2(m_plane.Point(u, -v, math::vec(0, 0, 0)));
		math::vec p3(m_plane.Point(u, v, math::vec(0, 0, 0)));
		math2vvr(math::Triangle(p0, p1, p2), colPlane).draw();
		math2vvr(math::Triangle(p2, p1, p3), colPlane).draw();
	}

	if (m_style_flag & FLAG_SHOW_SOLID) m_model.draw(m_obj_col, SOLID);
	if (m_style_flag & FLAG_SHOW_WIRE) m_model.draw(Colour::black, WIRE);
	if (m_style_flag & FLAG_SHOW_NORMALS) m_model.draw(Colour::black, NORMALS);
	if (m_style_flag & FLAG_SHOW_AXES) m_model.draw(Colour::black, AXES);

	//*
	for (int i = 0; i < m_points3D.size(); i++) {
		//m_points3D[i].draw();
	 }
	 //*/
	for (int i = 0; i < m_triangles3D.size(); i++) {
		m_triangles3D[i].draw();
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

