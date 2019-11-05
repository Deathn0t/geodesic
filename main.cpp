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

using namespace Eigen; // to use the classes provided by Eigen library
MatrixXd V1;           // matrix storing vertex coordinates of the input curve
MatrixXi F1;

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
  viewer.data(0).set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3));
}

void set_pc(igl::opengl::glfw::Viewer &viewer)
{
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().add_points(V1, Eigen::RowVector3d(0.3, 0.8, 0.3));
}

/**
 * Initialize a queue Q with a window for each edge adjacent to source vs.
 * input:
 *    id_vs: the source vertex id.
 *    V: vertexes of the mesh.
 *    E: edges of the mesh.
 * output:
 *  - Q: the queue to initialize.
 */
void init_Q_igl(int id_vs, MatrixXd &V, MatrixXi E, std::queue<Window> &Q)
{
  // TODO
  Q = std::queue<Window>();
  Window cur_w; // current window
  for (int e = 0; e < E.rows(); e++)
  {
    // D(E(e, 0), 0) = D(E(e, 0), 0) + 1;
    // D(E(e, 1), 0) = D(E(e, 1), 0) + 1;
    // std::cout << E(e, 0) << ", " << E(e, 1) << std::endl;
    if (id_vs == E(e, 0))
    {
      cur_w = Window();
    }
    else if (id_vs == E(e, 1))
    {
    }
  }
}

// void init_Q_He(queue<Window> &Q, Vector3d &vs)
// {
//   // TODO
// }

/**
 * Compute exact geodesics from a source vertex to very over vertices of the mesh.
 * input:
 *    F: faces of the mesh where to compute geodesics pathes
 *    V: vertexes of the mesh where to compute geodesics pathes
 *    id_vs: the source vertex id
 * output:
 *    ...: something with all windows computed to apply backtracing?
 */
void exact_geodesics_igl(MatrixXi &F, MatrixXd &V, int id_vs)
{
  // initialize the queue Q with a window for each edge adjacent
  // to source: source_v_t
  std::queue<Window> Q; // Q.push(..); Q.front(); Q.pop();
  Window cur_w;         // current window during iterations
  MatrixXi E;
  igl::edges(F, E); //! computed once here for efficiency

  init_Q_igl(id_vs, V, E, Q); // init Q with libigl DS

  while (!Q.empty())
  {
    // select and remove a Window from Q
    cur_w = Q.front();
    Q.pop();

    // propagate selected window
    // TODO

    // update queue with new windows
    // TODO
  }
}

/**
 * Retrieve the geodesics path from a source vertex to an end vertex, applying
 * backtracing on the result of exact_geodesics.
 * input:
 *    vs: the source/start vertex for the geodesic path.
 *    ve: the end/arrival vertex for the geodescic path.
 * output:
 *    ...: a kind of a list probably which represents the path?
 */
void retrieve_path(Vector3d &vs, Vector3d &ve)
{
  // TODO
}

void example_1()
{
  igl::readOFF("../data/star.off", V1, F1);
  int id_vs = 0;

  exact_geodesics_igl(F1, V1, id_vs);

  igl::opengl::glfw::Viewer viewer;
  set_meshes(viewer);
  viewer.launch();
}

int main(int argc, char *argv[])
{
  example_1();
}
