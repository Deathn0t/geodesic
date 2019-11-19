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
    Vector3d v1;
        // direction of s fron the directed edge e

public:
    Window() {}

    Window(double f_b0, double f_b1, double f_d0, double f_d1, double f_sigma, int f_dir, int f_edge_id,Vector3d &f_v0,Vector3d &f_v1)
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
};