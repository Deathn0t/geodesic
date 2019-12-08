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

#include "window.h"

using namespace Eigen;
using namespace std;

// const double EPS = 1e-7;

Vector2d intersect(Vector2d u, Vector2d v);

Vector2d compute_line(Vector2d u, Vector2d v);

/**
 * Check if c is in range [a,b]
**/
bool point_in_range(Vector2d c, Vector2d a, Vector2d b);
