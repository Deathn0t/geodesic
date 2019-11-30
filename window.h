#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;
using namespace std;

class Window
{
private:
    double b0, b1;
    double d0, d1;
    double sigma; // geodesic distance from s to the source vs
    int dir;
    int edge_id;
    Vector3d v0;
    int v0id;
    Vector3d v1;
    int v1id;

    // direction of s fron the directed edge e

public:
    Window() {}

    Window(double f_b0, double f_b1, double f_d0, double f_d1, double f_sigma, int f_dir, int f_edge_id, Vector3d &f_v0, Vector3d &f_v1, int f_v0id, int f_v1id)
    {
        b0 = f_b0;
        b1 = f_b1;
        d0 = f_d0;
        d1 = f_d1;
        sigma = f_sigma;
        dir = f_dir;
        edge_id = f_edge_id;
        v0 = f_v0;
        v1 = f_v1;
        v0id = f_v0id;
        v1id = f_v1id;
    }

    int get_edge_id()
    {
        return edge_id;
    }

    int get_dir()
    {
        return dir;
    }

    double get_sigma()
    {
        return sigma;
    }

    double get_b0()
    {
        return b0;
    }

    double get_b1()
    {
        return b1;
    }

    Vector3d get_v0()
    {
        return v0;
    }

    Vector3d get_v1()
    {
        return v1;
    }

    int get_v0id()
    {
        return v0id;
    }

    int get_v1id()
    {
        return v1id;
    }

    double get_d0()
    {
        return d0;
    }

    double get_d1()
    {
        return d1;
    }

    void print() {
        std::cout
            << "W(b0=" << b0
            << ", b1=" << b1
            << ", v0=[" << v0(0) << "," << v0(1) << "," << v0(2) << "]"
            << ", v1=[" << v1(0) << "," << v1(1) << "," << v1(2) << "]"
            << ")";
    }
};