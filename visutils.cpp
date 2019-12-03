#include "visutils.h"

void addColorEdge(igl::opengl::glfw::Viewer &viewer, Window &w, RowVector3d color, bool left, bool right)
{
    cout << "Color("
         << color(0) << ","
         << color(1) << ","
         << color(2) << ") ";
    w.print();
    cout << endl;

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