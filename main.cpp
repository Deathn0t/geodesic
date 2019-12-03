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
#include <tuple>

#include "HalfedgeBuilder.cpp"
#include "geoutils.h"
#include "visutils.h"
#include "window.h"

using namespace Eigen;
using namespace std;
// to use the classes provided by Eigen library
MatrixXd V1; // matrix storing vertex coordinates of the input curve
MatrixXi F1;
igl::opengl::glfw::Viewer VIEWER;
// const double EPS = 1e-7;

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
  viewer.data(0).show_faces = true;
  viewer.data(0).show_lines = false;
}

void set_pc(igl::opengl::glfw::Viewer &viewer)
{
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().add_points(V1, Eigen::RowVector3d(0.3, 0.8, 0.3));
}

Vector2d compute_pseudo_source(Window &w)
{
  // Computation vs (pseudo source) for new windows, intersection of 2 circles
  // the center of these two circles have the same y coord. for their center, thouse it simplify the resolution
  double x0, x1, d0, d1, c, delta, x, y_s1, y_s2;
  x0 = w.get_b0();
  x1 = w.get_b1();
  d0 = w.get_d0();
  d1 = w.get_d1();
  x = (d1 * d1 - d0 * d0 - x1 * x1 + x0 * x0) / (2 * (x0 - x1));
  c = x1 * x1 + x * x - 2 * x1 * x - d1 * d1;
  delta = sqrt(-4 * c);
  y_s1 = delta / 2;
  y_s2 = -delta / 2;

  // HERE VS SHOUDL BE CALLED S AND S SHOULD BE STORE IN WINDOW TO ALLOW FOR WINDOWS INTERSECTION.
  // we have the two possibles sources, we choosed the one with a positive y
  Vector2d s;
  if (y_s1 > 0.)
  {
    s = Vector2d(x, y_s1);
  }
  else
  {
    s = Vector2d(x, y_s2);
  }
  return s;
}

void add_window_Q(std::map<int, list<Window *> *> &e2w, std::queue<Window *> &Q, Vector3d &vs, Vector3d &v0, Vector3d &v1, int edge_id, int v0id, int v1id)
{
  double b0 = 0.;
  double b1 = (v1 - v0).norm();

  double d0 = (vs - v0).norm();
  double d1 = (vs - v1).norm();

  double sigma = 0;
  int dir = 0;

  double x0, x1, c, delta, x, y_s1, y_s2;

  x = (d1 * d1 - d0 * d0 - x1 * x1 + x0 * x0) / (2 * (x0 - x1));
  c = x1 * x1 + x * x - 2 * x1 * x - d1 * d1;
  delta = sqrt(-4 * c);
  y_s1 = delta / 2;
  y_s2 = -delta / 2;

  Vector2d s;
  if (y_s1 > 0.)
  {
    s = Vector2d(x, y_s1);
  }
  else
  {
    s = Vector2d(x, y_s2);
  }

  Window *newWindow;
  newWindow = new Window(b0, b1, d0, d1, s, sigma, dir, edge_id, v0, v1, v0id, v1id);
  addColorEdge(VIEWER, *newWindow, RowVector3d(0, 1, 0));
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
}

std::tuple<Vector2d, Vector2d, double>  computeIntersection(Window &leftWindow, Window &rightWindow)
{

  double intervalMin, intervalMax;
  intervalMin = rightWindow.get_b0();
  intervalMax = leftWindow.get_b1();

  double alpha, beta, gamma, A, B, C;
  Vector2d s_lw, s_rw;
  s_lw = compute_pseudo_source(leftWindow);
  s_rw = compute_pseudo_source(rightWindow);
  alpha = s_rw[0] - s_lw[0];
  beta = rightWindow.get_sigma() - leftWindow.get_sigma();
  gamma = s_rw.norm() * s_rw.norm() - s_lw.norm() * s_lw.norm() - beta * beta;
  A = alpha * alpha - beta * beta;
  B = gamma * alpha + 2 * s_rw[0] * beta * beta;
  C = (1 / 4.0) * gamma * gamma - s_rw.norm() * s_rw.norm() * beta * beta;

  double px1, px2, px;
  double delta = B * B - 4 * A * C;

  if (delta > 0.0)
  {
    px1 = (-B - sqrt(delta)) / (2 * A);
    px2 = (-B + sqrt(delta)) / (2 * A);

    if (intervalMin <= px1 && px1 <= intervalMax)
    {
      px = px1;
    }
    else
    {
      px = px2;
    }
  }
  else if (delta == 0.0)
  {
    px = -B / (2 * A);
  }
  else
  {
    std::cout << "IMPOSSIBLE BUT NO SOLUTION" << std::endl;
  }

  return std::make_tuple(s_lw,s_rw,px);
}

void push_window(Window &w, std::queue<Window *> &Q, std::map<int, list<Window *> *> &e2w)
{
  list<Window *> lw = *e2w[w.get_edge_id()];
  bool add_in_lw = true;

  // TODO: compute intersections with other window of same edges and replace if necessary
  std::cout << "For edge id " << w.get_edge_id() << ", the current windows are:" << std::endl;

  Window *curr_w;

  // CHECK IF INTERSECTION:
  for (Window *curr_w : lw) {
    std::cout << "One intersection:" << std::endl;
    
    curr_w->print();

    // 0 => new window, 1=> one window in list
    int leftWindow;

    double px;
    Vector2d s_lw,s_rw;

    if (w.get_b1() > curr_w->get_b0() && w.get_b0() < curr_w->get_b0())
    {
      /*           /\    /\
           //     /  \  /  \
           //    / w  \/  c \
           //   /     /\     \
               /     /  \     \  */
      leftWindow = 0;

      auto intersection_tuple = computeIntersection(w, *curr_w);
      s_lw = std::get<0>(intersection_tuple);
      s_rw = std::get<1>(intersection_tuple);
      px = std::get<2>(intersection_tuple);
      Vector2d px2d = Vector2d(px,0);

      w.set_b1(px);
      w.set_d1((s_lw - px2d).norm());

      curr_w->set_b0(px);
      curr_w->set_d0((s_rw - px2d).norm());
      
          
    }
    else if (curr_w->get_b1() > w.get_b0() && curr_w->get_b0() < w.get_b0())
    {
      /*           /\    /\
           //     /  \  /  \
           //    / c  \/  w \
           //   /     /\     \
               /     /  \     \  */
      leftWindow = 1;
      auto intersection_tuple = computeIntersection(*curr_w, w);
      s_lw = std::get<0>(intersection_tuple);
      s_rw = std::get<1>(intersection_tuple);
      px = std::get<2>(intersection_tuple);
      Vector2d px2d = Vector2d(px,0);

      curr_w->set_b1(px);
      curr_w->set_d1((s_lw - px2d).norm());

      w.set_b0(px);
      w.set_d0((s_rw-px2d).norm());

    }
    else if (w.get_b0() > curr_w->get_b0() && w.get_b1() < curr_w->get_b1())
    {
      // curr_w totally englobes w
      /*           /\            /\
           //     /c \         |/c \
           //    / /\ \        / \  \
           //   / /  \ \      /|  \  \
               / / w  \ \    / | w \  \  */
          std::cout<<"curr_w totally englobes w"<<std::endl;
    }
    else if (curr_w->get_b0() > w.get_b0() && curr_w->get_b1() < w.get_b1())
    {
      // w totally englobes curr_w
      /*           /\            /\
           //     /w \         |/w \
           //    / /\ \        / \  \
           //   / /  \ \      /|  \  \
               / / c  \ \    / | c \  \  */
          std::cout<<"w totally englobes curr_w"<<std::endl;

    }
    else if (curr_w->get_b0() > w.get_b1() || w.get_b0() > curr_w->get_b0())
    {
      std::cout << "NO INTERSECTION BETWEEN WINDOWS" << std::endl;
    }
  }

  // COMPARE DISTANCE AND DECIDE WHETHER THE WINDOW SHOULD BE ADDED

  add_in_lw = add_in_lw && w.get_d0() > 0. && w.get_d1() > 0.;
  add_in_lw = add_in_lw && abs(w.get_b0()-w.get_b1()) > EPS;

  if (add_in_lw)
  {
    lw.push_front(&w);
    Q.push(&w);
  }
}

void propagate_window(MatrixXd &V, HalfedgeDS &he, Window *p_w, std::queue<Window *> &Q, std::map<int, list<Window *> *> &e2w)
{
  Window w = *p_w;
  addColorEdge(VIEWER, w, RowVector3d(1, 0, 1));

  // we take the 3 points of our face in 3d representation

  int edgeid_p0p2, edgeid_p2p1, edgeid_p1p0;
  edgeid_p1p0 = he.getOpposite(w.get_edge_id());
  edgeid_p0p2 = he.getNext(edgeid_p1p0);
  edgeid_p2p1 = he.getNext(edgeid_p0p2);

  int p0id, p1id, p2id;
  p0id = he.getTarget(edgeid_p1p0);
  p1id = he.getTarget(edgeid_p2p1);
  p2id = he.getTarget(edgeid_p0p2);

  Vector3d p03d, p13d, p23d;
  p03d = w.get_v0();
  p13d = w.get_v1();
  p23d = V.row(p2id).transpose();

  // our x-axis is the (p0,p1) ligne then we compute the angle theta between (p0,p1) and (p0,p2) to get p2 in our 2D representation
  Vector3d u = p03d - p13d;
  Vector3d v = p03d - p23d;
  double theta = acos(u.dot(v) / (u.norm() * v.norm()));

  Vector2d p02d, p12d, p22d;
  p02d = Vector2d(0, 0);
  p12d = Vector2d((p13d - p03d).norm(), 0);
  p22d = Vector2d(cos(-theta) * (p23d - p03d).norm(), sin(-theta) * (p23d - p03d).norm());

  // Computation vs (pseudo source) for new windows, intersection of 2 circles
  // the center of these two circles have the same y coord. for their center, thouse it simplify the resolution
  Vector2d s = compute_pseudo_source(w);

  // solve linear system to find both lines
  // a*b0 + b = 0
  // a*vsx + b = vsy
  Vector2d l0, l1, lb; // two sides of the pencil of lights
  MatrixXd lA = Matrix2d(2, 2);
  lA(0, 0) = w.get_b0();

  lA(0, 1) = 1.;
  lA(1, 0) = s(0);
  lA(1., 1.) = 1.;
  lb(0) = 0;
  lb(1) = s(1);
  l0 = lA.colPivHouseholderQr().solve(lb); // first line

  lA(0, 0) = w.get_b1();
  l1 = lA.colPivHouseholderQr().solve(lb); // second line

  // line (p0,p2)  and line (p1,p2)
  Vector2d lp0p2, lp2p1;
  lp0p2 = Vector2d(p22d(1) / p22d(0), 0.);
  lA(0, 0) = p12d(0);
  lA(1, 0) = p22d(0);
  lb(1) = p22d(1);
  lp2p1 = lA.colPivHouseholderQr().solve(lb);

  // intersections
  Vector2d int_l0_lp0p2, int_l0_lp2p1, int_l1_lp0p2, int_l1_lp2p1;
  int_l0_lp0p2 = intersect(lp0p2, l0);
  int_l0_lp2p1 = intersect(lp2p1, l0);
  int_l1_lp0p2 = intersect(lp0p2, l1);
  int_l1_lp2p1 = intersect(lp2p1, l1);

  Window *pw;
  // filter different cases
  if (w.get_b0() < EPS || ((w.get_v1() - w.get_v0()).norm() - w.get_b1()) < EPS)
  {

    // the window w correspond to the whole edge
    if (w.get_b0() < EPS && ((w.get_v1() - w.get_v0()).norm() - w.get_b1()) < EPS) // case I - 1
    {

      // case where the whole face is inside the pencile of light, we create 2 windows
      if (!point_in_range(int_l0_lp2p1, p22d, p12d) &&
          !point_in_range(int_l1_lp0p2, p02d, p22d)) // case I - 1
      {
        pw = new Window(
            0,
            p22d.norm(),
            w.get_d0(),
            (s - p22d).norm(),
            s,
            w.get_sigma(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

        pw = new Window(
            0.,
            (p22d - p12d).norm(),
            (s - p22d).norm(),
            w.get_d1(),
            s,
            w.get_sigma(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
      }
      else if (point_in_range(int_l0_lp2p1, p22d, p12d)) // case I - 3
      {
        pw = new Window(
            (int_l0_lp2p1 - p22d).norm(),
            (p22d - p12d).norm(),
            (s - int_l0_lp2p1).norm(),
            w.get_d1(),
            s,
            w.get_sigma(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

        // p0 becomes the new pseudo source
        pw = new Window(
            0.,
            p22d.norm(),
            0.,
            p22d.norm(),
            p02d,
            w.get_sigma() + w.get_d0(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));

        pw = new Window(
            0.,
            (int_l0_lp2p1 - p22d).norm(),
            p22d.norm(),
            int_l0_lp2p1.norm(),
            p02d,
            w.get_sigma() + w.get_d0(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));
      }
      else if (point_in_range(int_l1_lp0p2, p02d, p22d)) // case I - 2
      {
        pw = new Window(
            0.,
            int_l1_lp0p2.norm(),
            w.get_d0(),
            (s - int_l1_lp0p2).norm(),
            s,
            w.get_sigma(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

        // p1 becomes the new pseudo source
        pw = new Window(
            int_l1_lp0p2.norm(),
            p22d.norm(),
            (p12d - int_l1_lp0p2).norm(),
            (p12d - p22d).norm(),
            p12d,
            w.get_sigma() + w.get_d1(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        // addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));

        pw = new Window(
            0.,
            (p12d - p22d).norm(),
            (p12d - p22d).norm(),
            0.,
            p12d,
            w.get_sigma() + w.get_d1(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        // addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));
      }
      else //! this case should not happen
      {
        std::cout << "error I" << std::endl;
      }
    }
    else if (w.get_b0() < EPS) // case II
    {
      if (!point_in_range(int_l0_lp2p1, p22d, p12d) && point_in_range(int_l1_lp2p1, p22d, p12d)) // case II - 1
      {
        pw = new Window(
            0,
            p22d.norm(),
            w.get_d0(),
            (s - p22d).norm(),
            s,
            w.get_sigma(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

        pw = new Window(
            0,
            (p22d - int_l1_lp2p1).norm(),
            (s - p22d).norm(),
            (s - int_l1_lp2p1).norm(),
            s,
            w.get_sigma(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
      }
      else if (point_in_range(int_l0_lp2p1, p22d, p12d) && point_in_range(int_l1_lp2p1, p22d, p12d)) // case II - 2
      {
        pw = new Window(
            (p22d - int_l0_lp2p1).norm(),
            (p22d - int_l1_lp2p1).norm(),
            (s - int_l0_lp2p1).norm(),
            (s - int_l1_lp2p1).norm(),
            s,
            w.get_sigma(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

        // p0 becomes the new pseudo source
        pw = new Window(
            0.,
            p22d.norm(),
            0.,
            p22d.norm(),
            p02d,
            w.get_sigma() + w.get_d0(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));

        pw = new Window(
            0.,
            (int_l0_lp2p1 - p22d).norm(),
            p22d.norm(),
            int_l0_lp2p1.norm(),
            p02d,
            w.get_sigma() + w.get_d0(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));
      }
      else if (point_in_range(int_l1_lp0p2, p02d, p22d)) // case II - 3
      {
        pw = new Window(
            0,
            int_l1_lp0p2.norm(),
            w.get_d0(),
            (s - int_l1_lp0p2).norm(),
            s,
            w.get_sigma(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
      }
      else //! this case should not happen
      {
        std::cout << "error case II" << std::endl;
      }
    }
    else if (((w.get_v1() - w.get_v0()).norm() - w.get_b1()) < EPS) // case II SYM
    {
      if (!point_in_range(int_l1_lp0p2, p02d, p22d) && point_in_range(int_l0_lp0p2, p02d, p22d)) // case II SYM - 1
      {
        pw = new Window(
            0,
            (p22d - p12d).norm(),
            (s - p22d).norm(),
            w.get_d1(),
            s,
            w.get_sigma(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

        pw = new Window(
            int_l0_lp0p2.norm(),
            p22d.norm(),
            (s - int_l0_lp0p2).norm(),
            (s - p22d).norm(),
            s,
            w.get_sigma(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
      }
      else if (point_in_range(int_l0_lp0p2, p02d, p22d) && point_in_range(int_l1_lp0p2, p02d, p22d)) // case II SYM - 2
      {
        pw = new Window(
            int_l0_lp0p2.norm(),
            int_l1_lp0p2.norm(),
            (s - int_l0_lp0p2).norm(),
            (s - int_l1_lp0p2).norm(),
            s,
            w.get_sigma(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

        // p1 becomes the new pseudo source
        pw = new Window(
            int_l1_lp0p2.norm(),
            p22d.norm(),
            (p12d - int_l1_lp0p2).norm(),
            (p12d - p22d).norm(),
            p12d,
            w.get_sigma() + w.get_d1(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));

        pw = new Window(
            0.,
            (p12d - p22d).norm(),
            (p12d - p22d).norm(),
            0.,
            p12d,
            w.get_sigma() + w.get_d1(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));
      }
      else if (point_in_range(int_l0_lp2p1, p22d, p12d)) // case II SYM - 3
      {
        pw = new Window(
            int_l0_lp2p1.norm(),
            (p22d - p12d).norm(),
            (s - int_l0_lp2p1).norm(),
            w.get_d1(),
            s,
            w.get_sigma(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
      }
      else //! this case should not happen
      {
        std::cout << "error case II SYM" << std::endl;
      }
    }
    else //! this case should not happen
    {
      std::cout << "error I - II - II SYM" << std::endl;
    }
  }
  else //case IV
  {
    if (point_in_range(int_l0_lp0p2, p02d, p22d) && point_in_range(int_l1_lp2p1, p22d, p12d)) // case IV - 1
    {
      pw = new Window(
          int_l0_lp0p2.norm(),
          p22d.norm(),
          (s - int_l0_lp0p2).norm(),
          (s - p22d).norm(),
          s,
          w.get_sigma(),
          0., edgeid_p0p2, p03d, p23d, p0id, p2id);
      push_window(*pw, Q, e2w);
      addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

      pw = new Window(
          0,
          (p22d - int_l1_lp2p1).norm(),
          (s - p22d).norm(),
          (s - int_l1_lp2p1).norm(),
          s,
          w.get_sigma(),
          0., edgeid_p2p1, p23d, p13d, p2id, p1id);
      push_window(*pw, Q, e2w);
      addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
    }
    else if (point_in_range(int_l0_lp0p2, p02d, p22d) && point_in_range(int_l1_lp0p2, p02d, p22d)) // case IV - 2
    {
      pw = new Window(
          int_l0_lp0p2.norm(),
          int_l1_lp0p2.norm(),
          (s - int_l0_lp0p2).norm(),
          (s - int_l1_lp0p2).norm(),
          s,
          w.get_sigma(),
          0., edgeid_p0p2, p03d, p23d, p0id, p2id);
      push_window(*pw, Q, e2w);
      addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
    }
    else if (point_in_range(int_l0_lp2p1, p22d, p12d) && point_in_range(int_l1_lp2p1, p22d, p12d)) // case IV - 3
    {
      pw = new Window(
          (p22d - int_l0_lp2p1).norm(),
          (p22d - int_l1_lp2p1).norm(),
          (s - int_l0_lp2p1).norm(),
          (s - int_l1_lp2p1).norm(),
          s,
          w.get_sigma(),
          0., edgeid_p2p1, p23d, p13d, p2id, p1id);
      push_window(*pw, Q, e2w);
      addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
    }
    else //! this case should not happen
    {
      std::cout << "error IV" << std::endl;
    }
  }
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
void exact_geodesics(HalfedgeDS &he, MatrixXd &V, MatrixXi &F, int id_vs)
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

  int it = 0;

  while (!Q.empty())
  {
    it++;

    // select and remove a Window from Q
    cur_w = Q.front();
    Q.pop();

    // propagate selected window
    propagate_window(V, he, cur_w, Q, e2w);

    if (it > 1000)
    {
      std::cout << "break after " << it << "iterations" << std::endl;
      return;
    }
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

void example_1(char *file)
{
  cout << "file: " << file << endl;
  igl::readOFF(file, V1, F1);
  int vs = 0;
  HalfedgeBuilder *builder = new HalfedgeBuilder();
  HalfedgeDS he = (builder->createMeshWithFaces(V1.rows(), F1));

  igl::opengl::glfw::Viewer &viewer = VIEWER;

  exact_geodesics(he, V1, F1, vs);
  viewer.data().add_points(V1.row(vs), Eigen::RowVector3d(1, 0, 0));
  set_meshes(viewer, V1, F1);
  viewer.launch();
}

void example_2(char *file)
{
  cout << "file: " << file << endl;
  igl::readOFF(file, V1, F1);
  int vs = 1195;
  HalfedgeBuilder *builder = new HalfedgeBuilder();
  HalfedgeDS he = (builder->createMeshWithFaces(V1.rows(), F1));

  igl::opengl::glfw::Viewer &viewer = VIEWER;

  exact_geodesics(he, V1, F1, vs);

  viewer.data()
      .add_points(V1.row(vs), Eigen::RowVector3d(1, 0, 0));
  set_meshes(viewer, V1, F1);
  viewer.launch();
}

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    char *file = "../data/../data/gargoyle_tri.off";
    example_2(file);
  }
  else
  {
    example_1(argv[1]);
  }
}
