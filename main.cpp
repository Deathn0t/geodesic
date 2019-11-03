#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/octree.h>
#include <igl/knn.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>

using namespace Eigen; // to use the classes provided by Eigen library
MatrixXd V1;           // matrix storing vertex coordinates of the input curve
MatrixXi F1;
MatrixXd V2; // matrix storing vertex coordinates of the input curve
MatrixXi F2;

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;

  if (key == '1')
  {
    viewer.data(0).clear(); // Clear should be called before drawing the mesh
    viewer.data(0).set_mesh(V1, F1);
    viewer.data(0).set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3)); // update the mesh (both coordinates and faces)
  }

  return false;
}

void set_meshes(igl::opengl::glfw::Viewer &viewer)
{
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().set_mesh(V1, F1);
  viewer.append_mesh();
  viewer.data().set_mesh(V2, F2);
  viewer.data(0).set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3));
  viewer.data(1).set_colors(Eigen::RowVector3d(0.8, 0.3, 0.3));
}

void set_pc(igl::opengl::glfw::Viewer &viewer)
{
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().add_points(V1, Eigen::RowVector3d(0.3, 0.8, 0.3));
}

void ex1()
{
  igl::readOFF("../data/star.off", V1, F1);
  igl::readOFF("../data/star_rotated.off", V2, F2);
  igl::opengl::glfw::Viewer viewer;
  set_meshes(viewer);
  viewer.launch();
}

int main(int argc, char *argv[])
{
  ex1();
}
