#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;

class Window
{
private:
    Vector3d b0, b1;
    double d0, d1;
    double sigma; // geodesic distance from s to the source vs
    int dir;      // direction of s fron the directed edge e

public:
    Window() {}

    Window(Vector3d b0, Vector3d b1, double d0, double d1, double sigma, int dir)
    {
        this->b0 = b0;
        this->b1 = b1;
        this->d0 = d0;
        this->d1 = d1;
        this->sigma = sigma;
        this->dir = dir;
    }
};