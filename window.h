#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;
using namespace std;

class Window
{
private:
    Vector3d b0, b1;
    double d0, d1;
    double sigma; // geodesic distance from s to the source vs
    int dir; 
    int edge_id;    // direction of s fron the directed edge e

public:
    Window() {}

    Window(Vector3d b0, Vector3d b1, double d0, double d1, double sigma, int dir, int edge_id)
    {
        b0 = b0;
        b1 = b1;
        d0 = d0;
        d1 = d1;
        sigma = sigma;
        dir = dir;
        edge_id = edge_id;
    }
};