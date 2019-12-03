#pragma once

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

void addColorEdge(igl::opengl::glfw::Viewer &viewer, Window &w, RowVector3d color = RowVector3d(1, 0, 0), bool left = true, bool right = true);

void colorWindow(igl::opengl::glfw::Viewer &viewer, Window &w, RowVector3d color = RowVector3d(1, 0, 0), bool left = true, bool right = true);