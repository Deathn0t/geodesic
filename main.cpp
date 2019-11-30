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

#include "HalfedgeBuilder.cpp"

#include "window.h"

using namespace Eigen;
using namespace std;
// to use the classes provided by Eigen library
MatrixXd V1; // matrix storing vertex coordinates of the input curve
MatrixXi F1;

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;

  if (key == '1')
  {

    viewer.data(0).clear(); // Clear should be called before drawing the mesh
    viewer.data(0).set_mesh(V1, F1);
    viewer.data(0).set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3));
    // update the mesh (both coordinates and faces)
  }

  return false;
}

void set_meshes(igl::opengl::glfw::Viewer &viewer, MatrixXd &V, MatrixXi &F)
{
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data(0).set_mesh(V, F);
  viewer.append_mesh();
  // viewer.data(0).show_faces = false;
  viewer.data(0).show_lines = false;
}

void addColorEdge(igl::opengl::glfw::Viewer &viewer, MatrixXd &V, Window &w, RowVector3d color = RowVector3d(1, 0, 0), bool left = true, bool right = true)
{
  std::cout << "Color("
            << color(0) << ","
            << color(1) << ","
            << color(2) << ") ";
  w.print();
  std::cout << std::endl;

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

void set_pc(igl::opengl::glfw::Viewer &viewer)
{
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().add_points(V1, Eigen::RowVector3d(0.3, 0.8, 0.3));
}

void add_window_Q(std::map<int, list<Window *> *> &e2w, std::queue<Window *> &Q, Vector3d &vs, Vector3d &v0, Vector3d &v1, int edge_id, int v0id, int v1id)
{
  double b0 = 0.;
  double b1 = (v1 - v0).norm();

  double d0 = (vs - v0).norm();
  double d1 = (vs - v1).norm();

  double sigma = 0;
  int dir = 0;

  Window *newWindow;
  newWindow = new Window(b0, b1, d0, d1, sigma, dir, edge_id, v0, v1, v0id, v1id);
  Q.push(newWindow);

  e2w[edge_id]->push_front(newWindow);
}

void init_edgeid2windows(std::map<int, list<Window *> *> &e2w, HalfedgeDS &he)
{
  for (int e_id = 0; e_id < he.sizeOfHalfedges(); e_id++)
  {
    if (e2w.find(e_id) == e2w.end())
    {
      e2w[e_id] = new list<Window *>();
      e2w[he.getOpposite(e_id)] = e2w[e_id];
    }
  }
}
/**
 * Initialize a queue Q with a window for each edge adjacent to source vs.
 * input:
 *    id_vs: the source vertex id.
 *    V: vertexes of the mesh.
 *    E: edges of the mesh.
 * output:
 *  - Q: the queue to initialize.
 *
 *
 */
void init_Q(HalfedgeDS &he, int id_vs, MatrixXd &V, std::queue<Window *> &Q, std::map<int, list<Window *> *> &e2w)
{

  Vector3d vs, v_init, v_b0, v_b1;
  int incidentEdge, nextEdge;
  int v0id, v1id;

  vs = V.row(id_vs);
  incidentEdge = he.getEdge(id_vs);
  nextEdge = he.getNext(incidentEdge);
  v_init = V.row(he.getTarget(nextEdge));
  v0id = he.getTarget(nextEdge);

  nextEdge = he.getNext(nextEdge);
  v1id = he.getTarget(nextEdge);
  v_b1 = V.row(he.getTarget(nextEdge));

  add_window_Q(e2w, Q, vs, v_init, v_b1, nextEdge, v0id, v1id);

  while (v_init != v_b1)
  {
    v_b0 = v_b1;
    v0id = v1id;
    nextEdge = he.getNext(he.getOpposite(he.getNext(nextEdge)));
    v_b1 = V.row(he.getTarget(nextEdge));
    v1id = he.getTarget(nextEdge);

    add_window_Q(e2w, Q, vs, v_b0, v_b1, nextEdge, v0id, v1id);
  }
  // std::cout<<Q.size()<<std::endl;
}

Vector2d intersect(Vector2d u, Vector2d v)
{
  Vector2d inter;
  // x coord. of the intersection between l1 and l2
  inter(0) = (u(1) - v(1)) / (v(0) - u(0));
  // y coord. of the intersection between l1 and l2
  inter(1) = v(0) * inter(0) + v(1);
  return inter;
}

void propagate_window(igl::opengl::glfw::Viewer &viewer, MatrixXd &V, HalfedgeDS &he, Window *p_w, std::queue<Window *> &Q, std::map<int, list<Window *> *> &e2w)
{
  Window w = *p_w;
  // addColorEdge(viewer, V, w);

  // we take the 3 points of our face in 3d representation

  int edgeid_p0p2, edgeid_p2p1, edgeid_p1p0;
  edgeid_p0p2 = he.getNext(he.getOpposite(w.get_edge_id()));
  edgeid_p2p1 = w.get_edge_id();
  edgeid_p1p0 = he.getOpposite(w.get_edge_id());

  int p0id, p1id, p2id;
  p0id = he.getTarget(edgeid_p1p0);
  p1id = he.getTarget(edgeid_p2p1);
  p2id = he.getTarget(edgeid_p0p2);

  Vector3d p03d, p13d, p23d;
  p03d = w.get_v0();
  p13d = w.get_v1();
  p23d = V.row(he.getTarget(edgeid_p0p2)).transpose();

  // our x-axis is the (p0,p1) ligne then we compute the angle theta between (p0,p1) and (p0,p2) to get p2 in our 2D representation
  Vector3d u = p03d - p13d;
  Vector3d v = p03d - p23d;
  double theta = acos(u.dot(v) / (u.norm() * v.norm()));

  Vector2d p02d, p12d, p22d;
  p02d = Vector2d(0, 0);
  p12d = Vector2d((p13d - p03d).norm(), 0);
  p22d = Vector2d(cos(theta) * (p23d - p03d).norm(), sin(theta) * (p23d - p03d).norm());

  // Computation vs (pseudo source) for new windows, intersection of 2 circles
  // the center of these two circles have the same y coord. for their center, thouse it simplify the resolution
  double x0, x1, d0, d1, c, delta, x, y_vs1, y_vs2;
  x0 = w.get_b0();
  x1 = w.get_b1();
  d0 = w.get_d0();
  d1 = w.get_d1();
  x = (d1 * d1 - d0 * d0 - x1 * x1 + x0 * x0) / (2 * (x0 - x1));
  c = x1 * x1 + x * x - 2 * x1 * x - d1 * d1;
  delta = sqrt(-4 * c);
  y_vs1 = delta / 2;
  y_vs2 = -delta / 2;

  // we have the two possibles sources, we choosed the one with a positive y
  Vector2d vs;
  if (y_vs1 > 0.)
  {
    vs = Vector2d(x, y_vs1);
  }
  else
  {
    vs = Vector2d(x, y_vs2);
  }

  // solve linear system to find both lines
  // a*b0 + b = 0
  // a*vsx + b = vsy
  Vector2d l0, l1, lb; // two sides of the pencil of lights
  MatrixXd lA = Matrix2d(2, 2);
  lA(0, 0) = w.get_b0();

  lA(0, 1) = 1.;
  lA(1, 0) = vs(0);
  lA(1., 1.) = 1.;
  lb(0) = 0;
  lb(1) = vs(1);
  l0 = lA.colPivHouseholderQr().solve(lb); // first line

  lA(0, 0) = w.get_b1();
  l1 = lA.colPivHouseholderQr().solve(lb); // second line

  // filter different cases
  if (w.get_b0() < 1e-10 or w.get_b1() < 1e-10)
  {
    if (w.get_b0() < 1e-10)
    {
      Window pw = Window(0., p22d.norm(), 0., p22d.norm(), w.get_sigma() + w.get_d0(), 0., edgeid_p0p2, p03d, p23d, p0id, p2id); // red window on left side (p0,p2)
      addColorEdge(viewer, V, pw, RowVector3d(0, 1, 0), false);

      Vector2d lp2p1; // line (p1,p2)
      lA(0, 0) = p12d(0);
      lA(1, 0) = p22d(0);
      lb(1) = p22d(1);
      lp2p1 = lA.colPivHouseholderQr().solve(lb);

      Vector2d int_l0_lp2p1 = intersect(lp2p1, l0);

      pw = Window(
          0.,
          (p22d - int_l0_lp2p1).norm(),
          p22d.norm(), int_l0_lp2p1.norm(),
          w.get_sigma() + w.get_d0(), 0.,
          edgeid_p2p1, p23d, p13d, p2id, p1id); // red window on left side (p2,p1)
      addColorEdge(viewer, V, pw, RowVector3d(0, 0, 1), false);

      Vector2d int_l1_lp2p1 = intersect(lp2p1, l1);
      pw = Window(
          (p22d - int_l0_lp2p1).norm(),
          (p22d - int_l1_lp2p1).norm(),
          w.get_d0() + int_l0_lp2p1.norm(),
          w.get_d1() + int_l1_lp2p1.norm(),
          w.get_sigma(),
          0.,
          edgeid_p2p1, p23d, p13d, p2id, p1id); // yellow window on side (p2,p1)
      addColorEdge(viewer, V, pw, RowVector3d(0, 1, 1), false);
    }
    else
    {                                                                                                                            // w.get_b1() < 1e-10 is true
      Window pw = Window(0., p22d.norm(), 0., p22d.norm(), w.get_sigma() + w.get_d0(), 0., edgeid_p0p2, p03d, p23d, p0id, p2id); // red window on left side (p0,p2)
      addColorEdge(viewer, V, pw);

      Vector2d lp2p1; // line (p1,p2)
      lA(0, 0) = p12d(0);
      lA(1, 0) = p22d(0);
      lb(1) = p22d(1);
      lp2p1 = lA.colPivHouseholderQr().solve(lb);

      Vector2d int_l0_lp2p1 = intersect(lp2p1, l0);

      pw = Window(
          0.,
          (p22d - int_l0_lp2p1).norm(),
          p22d.norm(), int_l0_lp2p1.norm(),
          w.get_sigma() + w.get_d0(), 0.,
          edgeid_p2p1, p13d, p23d, p1id, p2id); // red window on left side (p2,p1)
      addColorEdge(viewer, V, pw);

      Vector2d int_l1_lp2p1 = intersect(lp2p1, l1);
      pw = Window(
          (p22d - int_l0_lp2p1).norm(),
          (p22d - int_l1_lp2p1).norm(),
          w.get_d0() + int_l0_lp2p1.norm(),
          w.get_d1() + int_l1_lp2p1.norm(),
          w.get_sigma(),
          0.,
          edgeid_p2p1, p13d, p23d, p1id, p2id); // yellow window on side (p2,p1)
    }
  }
  else
  {

    Vector2d lp0p2 = Vector2d(p22d(1) / p22d(0), 0.); // line (p0,p2)

    Vector2d lp1p2; // line (p1,p2)
    lA(0, 0) = p12d(0);
    lA(1, 0) = p22d(0);
    lb(1) = p22d(1);
    lp1p2 = lA.colPivHouseholderQr().solve(lb);

    Vector2d int_l0_lp0p2;
    // x coord. of the intersection between l0 and (p0,p2)
    int_l0_lp0p2(0) = (lp0p2(1) - l0(1)) / (lp0p2(0) - l0(0));
    // y coord. of the intersection between l0 and (p0,p2)
    int_l0_lp0p2(1) = l0(0) * int_l0_lp0p2(0) + l0(1);

    Vector2d int_l1_lp0p2;
    // x coord. of the intersection between l1 and (p0,p2)
    int_l1_lp0p2(0) = (lp0p2(1) - l1(1)) / (lp0p2(0) - l1(0));
    // y coord. of the intersection between l1 and (p0,p2)
    int_l1_lp0p2(1) = l1(0) * int_l1_lp0p2(0) + l1(1);

    std::cout << "int. l0 and lp0p2: " << int_l0_lp0p2 << std::endl;
    std::cout << "int. l1 and lp0p2: " << int_l1_lp0p2 << std::endl;
  }
  //http://math.15873.pagesperso-orange.fr/IntCercl.html
  /* HOW TO PROPAGATE?

    1. Transpose old b0,b1 to new bo,b1 new plane (plane of edge),
    2. Take vector of d0, d1 from old face and intersect with new edge (vector = x=infiny, y =0)
    3. From this vectors, find intersections corresponding to new b0,b1.
    4. From this new b0,b1 find distance d0' , d1'.
    5. Add d0', d1' to d0  and d1.
    6. Compute new source using new distances and new b0,b1.
    7. This give the new window : new b0 to b1), new distances to source d0,d1, edge_id

    */

  // DO PLANE TRANSLATIONr
}

/**
 * Compute exact geodesics from a source vertex to very over vertices of the mesh.
 * input:
 *    F: faces of the mesh where to compute geodesics pathes
 *    V: vertexes of the mesh where to compute geodesics pathes
 *    id_vs: the source vertex id
 * output:
 *    ...: something with all windows computed to apply backtracing?
 */
void exact_geodesics(igl::opengl::glfw::Viewer &viewer, HalfedgeDS &he, MatrixXd &V, MatrixXi &F, int id_vs)
{

  // initialize the queue Q with a window for each edge adjacent
  // to source: source_v_t
  std::queue<Window *> Q;              // Q.push(..); Q.front(); Q.pop();
  std::map<int, list<Window *> *> e2w; // map edge id to windows

  init_edgeid2windows(e2w, he);
  std::cout << "init edgeid2windows(): OK" << std::endl;
  init_Q(he, id_vs, V, Q, e2w);
  std::cout << "init init_Q(): OK" << std::endl;

  Window *cur_w; // current window during iterations

  while (!Q.empty())
  {
    // select and remove a Window from Q
    cur_w = Q.front();
    Q.pop();

    // propagate selected window
    propagate_window(viewer, V, he, cur_w, Q, e2w);
    // return;
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
  igl::readOFF("../data/sphere.off", V1, F1);
  int vs = 0;
  HalfedgeBuilder *builder = new HalfedgeBuilder();
  HalfedgeDS he = (builder->createMeshWithFaces(V1.rows(), F1));

  igl::opengl::glfw::Viewer viewer;

  exact_geodesics(viewer, he, V1, F1, vs);
  viewer.data().add_points(V1.row(vs), Eigen::RowVector3d(1, 0, 0));
  set_meshes(viewer, V1, F1);
  viewer.launch();
}

// void print_l(list<int> l)
// {
//   for (auto v : l)
//     std::cout << v << ", ";
//   std::cout << std::endl;
// }

int main(int argc, char *argv[])
{
  example_1();
  // std::map<int, list<int> *> l;
  // l[0] = new list<int>();
  // l[1] = l[0];
  // l[0]->push_front(1);
  // std::cout << "l0: ";
  // print_l(*l[0]);
  // std::cout << "l1: ";
  // print_l(*l[1]);
}
