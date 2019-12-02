#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/octree.h>
#include <igl/knn.h>
#include <igl/writeOBJ.h>
#include <igl/edges.h>
#include <iostream>
#include <ostream>
#include <queue>

#include "window.cpp"

using namespace Eigen;
using namespace std;

void addColorEdge(igl::opengl::glfw::Viewer &viewer, Window &w, RowVector3d color = RowVector3d(1, 0, 0), bool left = true, bool right = true)
{
    std::cout << "Color("
              << color(0) << ","
              << color(1) << ","
              << color(2) << ") ";
    //w.print();
    std::cout << std::endl;

    MatrixXi edgeVertices(1, 2);
    edgeVertices(0, 0) = w.get_v0id();
    edgeVertices(0, 1) = w.get_v1id();

    Vector3d normalized = (w.get_v1() - w.get_v0()) / (w.get_v1() - w.get_v0()).norm();

    MatrixXd point1(1, 3), point2(1, 3);
    point1.row(0) = w.get_v0() + w.get_b0() * normalized;
    point2.row(0) = w.get_v0() + w.get_b1() * normalized;

    if (left)
    {
        viewer.data(0).add_points(point1, color);
    }

    if (right)
    {

        viewer.data(0).add_points(point2, color);
    }

    viewer.data(0).add_edges(point1, point2, color);
}