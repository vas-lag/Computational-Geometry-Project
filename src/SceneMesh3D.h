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
	void Task1(std::vector<vvr::Triangle>& m_triangles, std::vector<vec>& m_vertices);
	void Task4(std::vector<vvr::Triangle>& m_triangles, std::vector<vec>& m_vertices_new, const Eigen::MatrixXd& Coords, const Eigen::MatrixXd& CoordsNew, const Eigen::MatrixXd& DifCoords, const Eigen::MatrixXd& DifCoordsNew, const Eigen::MatrixXd& L);

private:
	int m_style_flag;
	float m_plane_d;
	vvr::Canvas2D m_canvas;
	vvr::Colour m_obj_col;
	vvr::Mesh m_model_original, m_model, m_model_new;
	vvr::Box3D m_aabb;
	math::vec m_center_mass;
	math::vec m_pca_cen;
	math::vec m_pca_dir;
	math::Plane m_plane;
	std::vector<int> m_intersections;
	std::vector<vvr::Point3D> m_points3D, m_points_coords3D, m_points_difCoords3D;
	std::vector<vvr::Triangle3D> m_triangles3D, m_triangles_coords3D, m_triangles_difCoords3D;
};

