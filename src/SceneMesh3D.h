#include <VVRScene/canvas.h>
#include <VVRScene/mesh.h>
#include <VVRScene/settings.h>
#include <VVRScene/utils.h>
#include <MathGeoLib.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <set>
#include "symmetriceigensolver3x3.h"
#include "canvas.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Eigenvalues"
#include "Eigen/SparseCholesky"

#define FLAG_SHOW_AXES       1
#define FLAG_SHOW_WIRE       2
#define FLAG_SHOW_SOLID      4
#define FLAG_SHOW_NORMALS    8
#define FLAG_SHOW_PLANE     16
#define FLAG_SHOW_AABB      32
#define FLAG_SHOW_DIFCOORDS 64
#define FLAG_SHOW_ORIGINAL_MODEL 128
#define FLAG_SHOW_POINTS 256
#define FLAG_SHOW_TRIANGLES 512

class Mesh3DScene : public vvr::Scene
{
public:
	Mesh3DScene();
	const char* getName() const { return "3D Scene"; }
	void keyEvent(unsigned char key, bool up, int modif) override;
	void arrowEvent(vvr::ArrowDir dir, int modif) override;

private:
	void draw() override;
	void reset() override;
	void resize() override;
	void Tasks();
	void Task1(std::vector<vvr::Triangle>& m_triangles, std::vector<vec>& m_vertices, Eigen::SparseMatrix<double>& I, Eigen::SparseMatrix<double>& A, Eigen::SparseMatrix<double>& D, Eigen::SparseMatrix<double>& D_inverse, Eigen::SparseMatrix<double>& L, Eigen::MatrixXd& Coords, Eigen::MatrixXd& DifCoords);
	void Task2(Eigen::MatrixXd& Ls, Eigen::VectorXd& eigenValues, Eigen::MatrixXd& eigenVectors, int verticesCount);
	void Task3(Eigen::MatrixXd& Q, Eigen::MatrixXd& Coords, Eigen::MatrixXd& DifCoords, Eigen::MatrixXd& CoordsNewA, Eigen::MatrixXd& DifCoordsNewA, Eigen::MatrixXd& CoordsNewB, Eigen::MatrixXd& DifCoordsNewB, Eigen::MatrixXd& CoordsNewC, Eigen::MatrixXd& DifCoordsNewC, Eigen::MatrixXd& CoordsNewD, Eigen::MatrixXd& DifCoordsNewD, Eigen::MatrixXd& Ls, int verticesCount);
	void Task3Sub(Eigen::MatrixXd& Q, Eigen::MatrixXd& Coords, Eigen::MatrixXd& DifCoords, Eigen::MatrixXd& Ls, Eigen::MatrixXd& CoordsNew, Eigen::MatrixXd& DifCoordsNew, std::vector<vec>& m_vertices_new, int verticesCount, double ratio);
	void Task4(const Eigen::MatrixXd& Coords, const Eigen::MatrixXd& CoordsNewA, const Eigen::MatrixXd& CoordsNewB, const Eigen::MatrixXd& CoordsNewC, const Eigen::MatrixXd& CoordsNewD, const Eigen::MatrixXd& DifCoords, const Eigen::MatrixXd& DifCoordsNewA, const Eigen::MatrixXd& DifCoordsNewB, const Eigen::MatrixXd& DifCoordsNewC, const Eigen::MatrixXd& DifCoordsNewD, const Eigen::MatrixXd& L);
	void Task4Sub(std::vector<vvr::Triangle>& m_triangles, std::vector<vec>& m_vertices_new, std::vector<vvr::Point3D>& m_points_coords_3D, std::vector<vvr::Triangle3D>& m_triangles_coords3D, std::vector<vvr::Point3D>& m_points_difCoords_3D, std::vector<vvr::Triangle3D>& m_triangles_difCoords3D, const Eigen::MatrixXd& Coords, const Eigen::MatrixXd& CoordsNew, const Eigen::MatrixXd& DifCoords, const Eigen::MatrixXd& DifCoordsNew, const Eigen::MatrixXd& L);
	void GetDifCoordsInNormalDirection(std::vector<vvr::Triangle>& m_triangles, std::vector<vec>& m_vertices, Eigen::MatrixXd& DifCoords, Eigen::MatrixXd& newDifCoords);
	void ComputeEigenDecomposition(Eigen::MatrixXd& Ls, Eigen::VectorXd& eigenValues, Eigen::MatrixXd& eigenVectors);
	void SaveEigenToFile(const std::string eigenFile, Eigen::VectorXd& eigenValues, Eigen::MatrixXd& eigenVectors, int verticesCount);
	void ReadEigenFromFile(const std::string eigenFile, Eigen::VectorXd& eigenValues, Eigen::MatrixXd& eigenVectors, int verticesCount);

private:
	int m_style_flag, ratio_flag;
	float m_plane_d;
	vvr::Canvas2D m_canvas;
	vvr::Colour m_obj_col;
	vvr::Mesh m_model_original, m_model, m_model_newA, m_model_newB, m_model_newC, m_model_newD;
	vvr::Mesh *m_model_new_draw;
	vvr::Box3D m_aabb;
	math::vec m_center_mass;
	math::vec m_pca_cen;
	math::vec m_pca_dir;
	math::Plane m_plane;
	std::vector<int> m_intersections;
	std::vector<vvr::Point3D> m_points3D, m_points_coords3DA, m_points_difCoords3DA, m_points_coords3DB, m_points_difCoords3DB, m_points_coords3DC, m_points_difCoords3DC, m_points_coords3DD, m_points_difCoords3DD;
	std::vector<vvr::Point3D> *m_points_coords3D_draw, *m_points_difCoords3D_draw;
	std::vector<vvr::Triangle3D> m_triangles3D, m_triangles_coords3DA, m_triangles_difCoords3DA, m_triangles_coords3DB, m_triangles_difCoords3DB, m_triangles_coords3DC, m_triangles_difCoords3DC, m_triangles_coords3DD, m_triangles_difCoords3DD;
	std::vector<vvr::Triangle3D> *m_triangles_coords3D_draw, *m_triangles_difCoords3D_draw;
};

void SparseIdentity(Eigen::SparseMatrix<double>& I, int n);
void SparseDiagonalInverse(Eigen::SparseMatrix<double>& D, Eigen::SparseMatrix<double>& D_inverse, int n);
double FindMax(Eigen::MatrixXd& M, int n, int index);
