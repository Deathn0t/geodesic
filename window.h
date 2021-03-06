#ifndef Window_H
#define Window_H
#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;
using namespace std;

const double EPS = 1e-10;

class Window
{
private:
    double b0, b1;
    double d0, d1;
    Vector2d *s;
    double sigma; // geodesic distance from s to the source vs
    int dir;
    int edge_id;
    Vector3d v0;
    int v0id;
    Vector3d v1;
    int v1id;

    // direction of s fron the directed edge e

public:
    Window();

    Window(double f_b0, double f_b1, double f_d0, double f_d1, double f_sigma, int f_dir, int f_edge_id, Vector3d &f_v0, Vector3d &f_v1, int f_v0id, int f_v1id);

    Window(Window &cw);

    int get_edge_id();

    int get_dir();

    double get_sigma();

    double get_b0();

    double get_b1();

    void set_b0(double new_b0);

    void set_b1(double new_b1);

    Vector3d get_v0();

    Vector3d get_v1();

    int get_v0id();

    int get_v1id();

    double get_d0();

    double get_d1();

    void set_d0(double f_d0);

    void set_d1(double f_d1);

    Vector2d get_s();

    double norm_edge();

    double min_geodist();

    double min_geodist(Vector2d inter);

    double max_geodist();

    double max_geodist(Vector2d inter);

    double geodist_b0();

    double geodist_b1();

    double geodist(double x);

    void print();
};

#endif // Included_NameModel_H