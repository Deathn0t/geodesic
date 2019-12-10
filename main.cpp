#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/octree.h>
#include <igl/knn.h>
#include <igl/writeOBJ.h>
#include <igl/edges.h>
#include <igl/exact_geodesic.h>
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
int ID_VS, ID_VT;
HalfedgeDS *HE;

struct GreaterThanByDist
{
  bool operator()(Window *w1, Window *w2) const
  {
    // return w1->get_sigma() > w2->get_sigma();
    return w1->min_geodist() > w2->min_geodist();
  }
};

void remove_from_queue(priority_queue<Window *, vector<Window *>, GreaterThanByDist> &Q, Window *pw)
{
  std::queue<Window *> q;
  Window *curr_pw;
  for (int i = 0; i < Q.size(); i++)
  {
    curr_pw = Q.top();
    Q.pop();

    if (curr_pw == pw)
    {
      break;
    }
    else
    {
      q.push(curr_pw);
    }
  }
  for (int i = 0; i < q.size(); i++)
  {
    curr_pw = q.front();
    q.pop();
    Q.push(curr_pw);
  }
}

void add_window_Q(std::map<int, list<Window *> *> &e2w, priority_queue<Window *, vector<Window *>, GreaterThanByDist> &Q, Vector3d &vs, Vector3d &v0, Vector3d &v1, int edge_id, int v0id, int v1id)
{
  double b0 = 0.;
  double b1 = (v1 - v0).norm();

  double d0 = (vs - v0).norm();
  double d1 = (vs - v1).norm();

  double sigma = 0;
  int dir = 0;

  Window *newWindow;
  newWindow = new Window(b0, b1, d0, d1, sigma, dir, edge_id, v0, v1, v0id, v1id);
  addColorEdge(VIEWER, *newWindow, RowVector3d(0, 1, 0));
  Q.push(newWindow);

  e2w[edge_id]->push_front(newWindow);
}

void init_edgeid2windows(std::map<int, list<Window *> *> &e2w, HalfedgeDS &he)
{
  for (int e_id = 0; e_id < he.sizeOfHalfedges(); e_id++)
  {
    e2w[e_id] = new list<Window *>();
  }
}
/**
 * Initialize a queue Q with a window for each edge adjacent to source vs.
 */
void init_Q(HalfedgeDS &he, int id_vs, priority_queue<Window *, vector<Window *>, GreaterThanByDist> &Q, std::map<int, list<Window *> *> &e2w)
{

  Vector3d vs, v_init, v_b0, v_b1;
  int incidentEdge, nextEdge;
  int v0id, v1id;

  vs = V1.row(id_vs);
  incidentEdge = he.getEdge(id_vs);
  nextEdge = he.getNext(incidentEdge);
  v_init = V1.row(he.getTarget(nextEdge));
  v0id = he.getTarget(nextEdge);

  nextEdge = he.getNext(nextEdge);
  v1id = he.getTarget(nextEdge);
  v_b1 = V1.row(he.getTarget(nextEdge));

  add_window_Q(e2w, Q, vs, v_init, v_b1, nextEdge, v0id, v1id);

  while (v_init != v_b1)
  {
    v_b0 = v_b1;
    v0id = v1id;
    nextEdge = he.getNext(he.getOpposite(he.getNext(nextEdge)));
    v_b1 = V1.row(he.getTarget(nextEdge));
    v1id = he.getTarget(nextEdge);

    add_window_Q(e2w, Q, vs, v_b0, v_b1, nextEdge, v0id, v1id);
  }
}

Vector2d range_intersect(Window &w0, Window &w1)
{
  double a0, a1, b0, b1, mini, maxi;
  a0 = w0.get_b0();
  a1 = w1.get_b0();
  b0 = w0.get_b1();
  b1 = w1.get_b1();
  if (a1 > a0)
  {
    mini = a1;
  }
  else
  {
    mini = a0;
  }
  if (b1 > b0)
  {
    maxi = b0;
  }
  else
  {
    maxi = b1;
  }
  return Vector2d(mini, maxi); //((a1 < b0 && a0 < b1) || (a0 < b1 && a1 < b0));
}

std::tuple<Vector2d, Vector2d, double> point_equidist(Window &leftWindow, Window &rightWindow, Vector2d inter)
{

  double alpha, beta, gamma, A, B, C;
  Vector2d s_lw, s_rw;
  alpha = s_rw[0] - s_lw[0];
  beta = rightWindow.get_sigma() - leftWindow.get_sigma();
  gamma = (s_lw.norm() * s_lw.norm()) - (s_rw.norm() * s_rw.norm()) - (beta * beta);
  A = (alpha * alpha) - (beta * beta);
  B = (gamma * alpha) + (2 * s_rw[0] * beta * beta);
  C = (1.0 / 4.0) * (gamma * gamma) - (s_rw.norm() * s_rw.norm() * beta * beta);

  double px1, px2, px;
  double delta = (B * B) - (4 * A * C);

  // cout << "delta: " << delta << endl;
  if (abs(delta) <= EPS)
  {
    px = -B / (2 * A);
  }
  else if (delta > 0)
  {
    px1 = (-B - sqrt(delta)) / (2 * A);
    px2 = (-B + sqrt(delta)) / (2 * A);
    if (inter(0) <= px1 && px1 <= inter(1))
    {
      px = px1;
    }
    else
    {
      px = px2;
    }
  }
  else
  {
    cout << "IMPOSSIBLE BUT NO SOLUTION" << endl;
    cout << "A=" << A << endl;
    cout << "B=" << B << endl;
    cout << "C=" << C << endl;
    cout << "leftWindow: " << endl;
    leftWindow.print();
    cout << endl;
    cout << "rightWindow: " << endl;
    rightWindow.print();
    cout << endl;
    // exit(0);
  }

  px = max(px, inter(0));
  px = min(px, inter(1));

  return std::make_tuple(s_lw, s_rw, px);
}

bool intersect(Window &w0, Window &w1)
{
  double a0, a1, b0, b1;
  a0 = w0.get_b0();
  a1 = w1.get_b0();
  b0 = w0.get_b1();
  b1 = w1.get_b1();
  return ((a1 < b0 && a0 < b1) || (a0 < b1 && a1 < b0));
}

void push_window(Window &w, priority_queue<Window *, vector<Window *>, GreaterThanByDist> &Q, map<int, list<Window *> *> &e2w)
{
  cout << "PUSH: " << endl;
  w.print();
  cout << endl;
  list<Window *> &lw = *e2w[w.get_edge_id()];
  list<Window *> copy_lw = *e2w[w.get_edge_id()];
  bool add_in_Q = true;
  bool add_in_lw = true;

  // cout
  // << "For edge id " << w.get_edge_id() << " evaluating conflicts:" << endl;

  Window *curr_w;
  // CHECK IF INTERSECTION:
  cout << "start loop" << endl;
  for (Window *curr_w : copy_lw)
  {
    // 0 => new window, 1=> one window in list
    Vector2d inter = range_intersect(w, *curr_w);
    if (intersect(w, *curr_w) && abs(inter(0) - inter(1)) > EPS)
    {
      cout << " if 1" << endl;

      double min_dist_w_s = w.min_geodist(inter);
      double max_dist_w_s = w.max_geodist(inter);
      double min_dist_curr_w_s = curr_w->min_geodist(inter);
      double max_dist_curr_w_s = curr_w->max_geodist(inter);

      if (max_dist_w_s <= min_dist_curr_w_s)
      {
        cout << "sub if 1" << endl;
        // w is always better than curr_w, replace curr_w par w
        if (w.get_b1() > curr_w->get_b0() && w.get_b0() < curr_w->get_b0() + EPS)
        {
          cout << "#1" << endl;
          cout << "w: " << endl;
          w.print();
          cout << endl;
          cout << "curr_w: " << endl;
          curr_w->print();
          cout << endl;
          cout << "inter: " << inter << endl;
          // w is the left window
          curr_w->set_d0((curr_w->get_s() - Vector2d(w.get_b1(), 0)).norm());
          curr_w->set_b0(w.get_b1());
          cout << "#2" << endl;
        }
        else if (curr_w->get_b1() > w.get_b0() && curr_w->get_b0() < w.get_b0() + EPS)
        {
          // curr_w is the left window
          curr_w->set_d1((curr_w->get_s() - Vector2d(w.get_b0(), 0)).norm());
          curr_w->set_b1(w.get_b0());
        }
      }
      else if (max_dist_curr_w_s <= min_dist_w_s)
      {
        cout << "sub if 2" << endl;
        // curr_w is always better than w
        if (w.get_b1() > curr_w->get_b0() && w.get_b0() < curr_w->get_b0() + EPS)
        {
          // w is the left window
          cout << "### 1" << endl;
          cout << "w: " << endl;
          w.print();
          cout << endl;
          cout << "curr_w: " << endl;
          curr_w->print();
          cout << endl;
          w.set_d1((w.get_s() - Vector2d(curr_w->get_b0(), 0)).norm());
          w.set_b1(curr_w->get_b0());
          cout << "### 2" << endl;
        }
        else if (curr_w->get_b1() > w.get_b0() && curr_w->get_b0() < w.get_b0() + EPS)
        {
          // curr_w is the left window
          cout << "### 3" << endl;
          cout << "w: " << endl;
          w.print();
          cout << endl;
          cout << "curr_w: " << endl;
          curr_w->print();
          cout << endl;
          w.set_d0((w.get_s() - Vector2d(curr_w->get_b1(), 0)).norm());
          w.set_b0(curr_w->get_b1());
          cout << "### 4" << endl;
        }
      }
      else
      {
        cout << "sub else" << endl;
        auto intersection_tuple = point_equidist(w, *curr_w, inter);
        double px;
        Vector2d s_lw, s_rw, px2d;
        s_lw = std::get<0>(intersection_tuple);
        s_rw = std::get<1>(intersection_tuple);
        px = std::get<2>(intersection_tuple);
        cout << "w: " << endl;
        w.print();
        cout << endl;
        cout << "curr_w: " << endl;
        curr_w->print();
        cout << endl;
        cout << "min_dist_w_s: " << min_dist_w_s << endl;
        cout << "max_dist_w_s: : " << max_dist_w_s << endl;
        cout << "min_dist_curr_w_s: " << min_dist_curr_w_s << endl;
        cout << "max_dist_curr_w_s: " << max_dist_curr_w_s << endl;
        cout << "inter: " << endl
             << inter << endl;
        cout << "px: " << px << endl;

        if (w.get_b1() > curr_w->get_b0() && w.get_b0() < curr_w->get_b0() + EPS)
        {
          // w is the left window
          cout << "if" << endl;
          w.set_d1((w.get_s() - px2d).norm());
          w.set_b1(px);
          cout << " ##1 " << endl;

          curr_w->set_d0((curr_w->get_s() - px2d).norm());
          curr_w->set_b0(px);
          cout << " ##2 " << endl;
        }
        else if (curr_w->get_b1() > w.get_b0() && curr_w->get_b0() < w.get_b0() + EPS)
        {
          // curr_w is the left window
          cout << "else if" << endl;
          w.set_d0((w.get_s() - px2d).norm());
          w.set_b0(px);

          curr_w->set_d1((curr_w->get_s() - px2d).norm());
          curr_w->set_b1(px);
        }
        else
        {
          cout << "dernier cas" << endl;
        }
      }
    }
    else
    {
      cout << "else 1" << endl;
    }
    if ((abs(curr_w->get_b1() - curr_w->get_b0()) <= EPS))
    {
      lw.remove(curr_w);
      // remove_from_queue(Q, curr_w);
    }
    else if (curr_w->get_d0() <= EPS || curr_w->get_d1() <= EPS)
    {
      // remove_from_queue(Q, curr_w);
    }
  }

  // w.print();
  // cout << endl;
  // cout << "w.get_d0=" << w.get_d0() << endl;
  // cout << "w.get_d1=" << w.get_d1() << endl;
  // COMPARE DISTANCE AND DECIDE WHETHER THE WINDOW SHOULD BE ADDED

  add_in_Q = add_in_Q && 0. < w.get_d0() && 0. < w.get_d1();
  add_in_Q = add_in_Q && abs(w.get_b1() - w.get_b0()) > EPS;
  add_in_lw = add_in_lw && (w.get_b1() - w.get_b0()) > EPS;

  if (add_in_Q)
  {
    Q.push(&w);
  }

  if (add_in_lw)
  {
    lw.push_front(&w);
  }
}

void propagate_window(HalfedgeDS &he, Window *p_w, priority_queue<Window *, vector<Window *>, GreaterThanByDist> &Q, map<int, list<Window *> *> &e2w)
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
  p23d = V1.row(p2id).transpose();

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
  Vector2d s = w.get_s();

  // solve linear system to find both lines
  // a*b0 + b = 0
  // a*vsx + b = vsy
  Vector2d l0, l1, lb; // two sides of the pencil of lights

  l0 = compute_line(Vector2d(w.get_b0(), 0), s);
  l1 = compute_line(Vector2d(w.get_b1(), 0), s);

  // line (p0,p2)  and line (p1,p2)
  Vector2d lp0p2, lp2p1;
  lp0p2 = Vector2d(p22d(1) / p22d(0), 0.);
  lp2p1 = compute_line(p22d, p12d);

  // intersections
  Vector2d int_l0_lp0p2, int_l0_lp2p1, int_l1_lp0p2, int_l1_lp2p1;
  int_l0_lp0p2 = intersect(lp0p2, l0);
  int_l0_lp2p1 = intersect(lp2p1, l0);
  int_l1_lp0p2 = intersect(lp0p2, l1);
  int_l1_lp2p1 = intersect(lp2p1, l1);

  if (isnan(s(1)))
  {
    // cout << endl
    //      << endl;
    // cout << "IF -> ";
    // w.print();
    // cout << "s: " << endl;
    // cout << s << endl;
    // cout
    //     << "l0: " << endl
    //     << l0 << endl;
    // cout << "l1: " << endl
    //      << l1 << endl;

    // cout << "int_l0_lp0p2: " << endl
    //      << int_l0_lp0p2 << endl;
    // cout << "int_l1_lp0p2: " << endl
    //      << int_l1_lp0p2 << endl;
    // cout << "int_l0_lp2p1: " << endl
    //      << int_l0_lp2p1 << endl;
    // cout << "int_l1_lp2p1: " << endl
    //      << int_l1_lp2p1 << endl;
    e2w[w.get_edge_id()]->remove(p_w);
    return;
  }
  else
  {
    cout << endl
         << endl;
    cout << "ELSE -> ";
    w.print();
    cout << endl;
    cout << "p02d: " << endl
         << p02d << endl;
    cout << "p12d: " << endl
         << p12d << endl;
    cout << "p22d: " << endl
         << p22d << endl;
    cout << "s: " << endl;
    cout << s << endl;
    cout
        << "l0: " << endl
        << l0 << endl;
    cout << "l1: " << endl
         << l1 << endl;

    cout << "int_l0_lp0p2: " << endl
         << int_l0_lp0p2 << endl;
    cout << "int_l1_lp0p2: " << endl
         << int_l1_lp0p2 << endl;
    cout << "int_l0_lp2p1: " << endl
         << int_l0_lp2p1 << endl;
    cout << "int_l1_lp2p1: " << endl
         << int_l1_lp2p1 << endl;
  }

  Window *pw;
  // filter different cases
  if (w.get_b0() <= EPS || ((w.get_v1() - w.get_v0()).norm() - w.get_b1()) <= EPS)
  {

    // the window w correspond to the whole edge
    if (w.get_b0() <= EPS && ((w.get_v1() - w.get_v0()).norm() - w.get_b1()) <= EPS) // case I - 1
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
            w.get_sigma(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        // addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

        pw = new Window(
            0.,
            (p22d - p12d).norm(),
            (s - p22d).norm(),
            w.get_d1(),
            w.get_sigma(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        // addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
      }
      else if (point_in_range(int_l0_lp2p1, p22d, p12d)) // case I - 3
      {
        pw = new Window(
            (int_l0_lp2p1 - p22d).norm(),
            (p22d - p12d).norm(),
            (s - int_l0_lp2p1).norm(),
            w.get_d1(),
            w.get_sigma(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        // addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

        // p0 becomes the new pseudo source
        pw = new Window(
            0.,
            p22d.norm(),
            0.,
            p22d.norm(),
            w.get_sigma() + w.get_d0(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        // addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));

        pw = new Window(
            0.,
            (int_l0_lp2p1 - p22d).norm(),
            p22d.norm(),
            int_l0_lp2p1.norm(),
            w.get_sigma() + w.get_d0(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        // addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));
      }
      else if (point_in_range(int_l1_lp0p2, p02d, p22d)) // case I - 2
      {
        pw = new Window(
            0.,
            int_l1_lp0p2.norm(),
            w.get_d0(),
            (s - int_l1_lp0p2).norm(),
            w.get_sigma(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        // addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

        // p1 becomes the new pseudo source
        pw = new Window(
            int_l1_lp0p2.norm(),
            p22d.norm(),
            (p12d - int_l1_lp0p2).norm(),
            (p12d - p22d).norm(),
            w.get_sigma() + w.get_d1(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        // addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));

        pw = new Window(
            0.,
            (p12d - p22d).norm(),
            (p12d - p22d).norm(),
            0.,
            w.get_sigma() + w.get_d1(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        // addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));
      }
      else //! this case should not happen
      {
        std::cout << "error case I" << std::endl;
        exit(0);
      }
    }
    else if (w.get_b0() <= EPS) // case II
    {
      if (!point_in_range(int_l0_lp2p1, p22d, p12d) && point_in_range(int_l1_lp2p1, p22d, p12d)) // case II - 1
      {
        pw = new Window(
            0,
            p22d.norm(),
            w.get_d0(),
            (s - p22d).norm(),
            w.get_sigma(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

        pw = new Window(
            0,
            (p22d - int_l1_lp2p1).norm(),
            (s - p22d).norm(),
            (s - int_l1_lp2p1).norm(),
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
            w.get_sigma() + w.get_d0(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));

        pw = new Window(
            0.,
            (int_l0_lp2p1 - p22d).norm(),
            p22d.norm(),
            int_l0_lp2p1.norm(),
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
            w.get_sigma(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
      }
      else //! this case should not happen
      {
        cout << "error case II" << endl;
        w.print();
        cout << "s: " << endl;
        cout << s << endl;
        cout << "l0: " << endl
             << l0 << endl;
        cout << "l1: " << endl
             << l1 << endl;

        cout << "p02d: " << endl
             << p02d << endl;
        cout << "p12d: " << endl
             << p12d << endl;
        cout << "p22d: " << endl
             << p22d << endl;

        cout << "int_l0_lp0p2: " << endl
             << int_l0_lp0p2 << endl;
        cout << "int_l1_lp0p2: " << endl
             << int_l1_lp0p2 << endl;
        cout << "int_l0_lp2p1: " << endl
             << int_l0_lp2p1 << endl;
        cout << "int_l1_lp2p1: " << endl
             << int_l1_lp2p1 << endl;

        cout << "bool 1: " << point_in_range(int_l0_lp2p1, p22d, p12d) << endl;
        cout << "bool 2: " << point_in_range(int_l1_lp2p1, p22d, p12d) << endl;
        exit(0);
      }
    }
    else if (((w.get_v1() - w.get_v0()).norm() - w.get_b1()) <= EPS) // case II SYM
    {
      cout << "else if case II SYM" << endl;
      if (!point_in_range(int_l1_lp0p2, p02d, p22d) && point_in_range(int_l0_lp0p2, p02d, p22d)) // case II SYM - 1
      {
        cout << "else if case II SYM - 1" << endl;
        pw = new Window(
            0,
            (p22d - p12d).norm(),
            (s - p22d).norm(),
            w.get_d1(),
            w.get_sigma(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

        pw = new Window(
            int_l0_lp0p2.norm(),
            p22d.norm(),
            (s - int_l0_lp0p2).norm(),
            (s - p22d).norm(),
            w.get_sigma(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
      }
      else if (point_in_range(int_l0_lp0p2, p02d, p22d) && point_in_range(int_l1_lp0p2, p02d, p22d)) // case II SYM - 2
      {
        cout << "else if case II SYM - 2" << endl;
        pw = new Window(
            int_l0_lp0p2.norm(),
            int_l1_lp0p2.norm(),
            (s - int_l0_lp0p2).norm(),
            (s - int_l1_lp0p2).norm(),
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
            w.get_sigma() + w.get_d1(),
            0., edgeid_p0p2, p03d, p23d, p0id, p2id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));

        pw = new Window(
            0.,
            (p12d - p22d).norm(),
            (p12d - p22d).norm(),
            0.,
            w.get_sigma() + w.get_d1(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(1, 0, 0));
      }
      else if (point_in_range(int_l0_lp2p1, p22d, p12d)) // case II SYM - 3
      {
        cout << "else if case II SYM - 3" << endl;
        pw = new Window(
            (p22d - int_l0_lp2p1).norm(),
            (p22d - p12d).norm(),
            (s - int_l0_lp2p1).norm(),
            w.get_d1(),
            w.get_sigma(),
            0., edgeid_p2p1, p23d, p13d, p2id, p1id);
        push_window(*pw, Q, e2w);
        addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
      }
      else //! this case should not happen
      {
        std::cout << "error case II SYM" << std::endl;
        w.print();
        cout << "s: " << endl;
        cout << s << endl;
        cout << "l0: " << endl
             << l0 << endl;
        cout << "l1: " << endl
             << l1 << endl;

        cout << "p02d: " << endl
             << p02d << endl;
        cout << "p12d: " << endl
             << p12d << endl;
        cout << "p22d: " << endl
             << p22d << endl;

        cout << "int_l0_lp0p2: " << endl
             << int_l0_lp0p2 << endl;
        cout << "int_l1_lp0p2: " << endl
             << int_l1_lp0p2 << endl;
        cout << "int_l0_lp2p1: " << endl
             << int_l0_lp2p1 << endl;
        cout << "int_l1_lp2p1: " << endl
             << int_l1_lp2p1 << endl;

        cout << "bool 1: " << point_in_range(int_l0_lp2p1, p22d, p12d) << endl;
        cout << "bool 2: " << point_in_range(int_l1_lp2p1, p22d, p12d) << endl;
        exit(0);
      }
    }
    else //! this case should not happen
    {
      std::cout << "error I - II - II SYM" << std::endl;
      exit(0);
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
          w.get_sigma(),
          0., edgeid_p0p2, p03d, p23d, p0id, p2id);
      push_window(*pw, Q, e2w);
      addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));

      pw = new Window(
          0,
          (p22d - int_l1_lp2p1).norm(),
          (s - p22d).norm(),
          (s - int_l1_lp2p1).norm(),
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
          w.get_sigma(),
          0., edgeid_p2p1, p23d, p13d, p2id, p1id);
      push_window(*pw, Q, e2w);
      addColorEdge(VIEWER, *pw, RowVector3d(0, 0, 1));
    }
    else //! this case should not happen
    {
      std::cout << "error IV" << std::endl;
      exit(0);
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
map<int, list<Window *> *> *exact_geodesics(HalfedgeDS &he, int id_vs)
{

  // initialize the queue Q with a window for each edge adjacent
  // to source: source_v_t
  priority_queue<Window *, vector<Window *>, GreaterThanByDist> Q;
  map<int, list<Window *> *> *e2w = new map<int, list<Window *> *>(); // map edge id to windows

  init_edgeid2windows(*e2w, he);
  cout << "init edgeid2windows(): OK" << endl;
  init_Q(he, id_vs, Q, *e2w);
  cout << "init init_Q(): OK" << endl;

  Window *cur_w; // current window during iterations

  int it = 0;

  while (!Q.empty())
  {
    it++;

    // select and remove a Window from Q
    cur_w = Q.top();
    Q.pop();

    // propagate selected window
    propagate_window(he, cur_w, Q, *e2w);

    if (it > 100000)
    {
      cout << "break after " << it << "iterations" << endl;
      break;
    }
    else if (it % 10000 == 0)
    {
      cout << "it: " << it << endl;
    }
  }

  return e2w;
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

double compute_geodist(int v, HalfedgeDS &he, map<int, list<Window *> *> &e2w)
{
  // TO BE COMPLETED
  Window *w = NULL;
  double best_geodist = std::numeric_limits<double>::max(); // INF

  int e = he.getEdge(v);
  int pEdge = e;
  int acc = 0;

  int in_e = e;
  int out_e = he.getOpposite(in_e);

  Window *cw;
  while (pEdge != e or acc == 0)
  {

    for (list<Window *>::iterator it = (*e2w[in_e]).begin(); it != (*e2w[in_e]).end(); ++it)
    {
      cw = *it;

      if (cw->get_v0() == V1.row(v).transpose())
      { // v is b0
        if (cw->get_b0() <= EPS)
        {
          if (w == NULL)
          {
            w = cw;
            best_geodist = w->geodist_b0();
          }
          else
          {
            if (cw->geodist_b0() < best_geodist)
            {
              w = cw;
              best_geodist = w->geodist_b0();
            }
          }
        }
      }
      else
      { // v is b1
        if (abs(cw->get_b1() - cw->norm_edge()) <= EPS)
        {
          if (w == NULL)
          {
            w = cw;
            best_geodist = w->geodist_b1();
          }
          else
          {
            if (cw->geodist_b1() < best_geodist)
            {
              w = cw;
              best_geodist = w->geodist_b1();
            }
          }
        }
      }
    }

    for (list<Window *>::iterator it = (*e2w[out_e]).begin(); it != (*e2w[out_e]).end(); ++it)
    {
      cw = *it;

      if (cw->get_v0() == V1.row(v).transpose())
      { // v is b0
        if (cw->get_b0() <= EPS)
        {
          if (w == NULL)
          {
            w = cw;
            best_geodist = w->geodist_b0();
          }
          else
          {
            if (cw->geodist_b0() < best_geodist)
            {
              w = cw;
              best_geodist = w->geodist_b0();
            }
          }
        }
      }
      else
      { // v is b1
        if (abs(cw->get_b1() - cw->norm_edge()) <= EPS)
        {
          if (w == NULL)
          {
            w = cw;
            best_geodist = w->geodist_b1();
          }
          else
          {
            if (cw->geodist_b1() < best_geodist)
            {
              w = cw;
              best_geodist = w->geodist_b1();
            }
          }
        }
      }
    }

    pEdge = he.getOpposite(he.getNext(pEdge));
    in_e = pEdge;
    out_e = he.getOpposite(in_e);
    acc++;
  }

  return best_geodist;
}

void set_meshes(igl::opengl::glfw::Viewer &viewer, MatrixXd &V, MatrixXi &F);

void show_source_target(igl::opengl::glfw::Viewer &viewer)
{
  viewer.data()
      .add_points(V1.row(ID_VS), RowVector3d(0, 1, 0));
  viewer.data()
      .add_points(V1.row(ID_VT), RowVector3d(0, 0.5, 0));
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;

  if (key == '1')
  {
    viewer.data().clear();
    map<int, list<Window *> *> *e2w = exact_geodesics(*HE, ID_VS);

    double d = compute_geodist(ID_VT, *HE, *e2w);
    cout << "PERSO(1): geodesic distance for target vertex " << ID_VT << " -> " << d << endl;

    set_meshes(viewer, V1, F1);
    for (map<int, list<Window *> *>::iterator it = e2w->begin(); it != e2w->end(); ++it)
    {
      for (Window *pw_i : *(it->second))
      {
        colorWindow(VIEWER, *pw_i, RowVector3d(1, 0, 0), false, false);
      }
    }
    show_source_target(viewer);
  }
  else if (key == '2')
  {
    viewer.data().clear();
    map<int, list<Window *> *> *e2w = exact_geodesics(*HE, ID_VS);

    Eigen::VectorXd d(V1.rows(), 1);
    for (int vt_i = 0; vt_i < V1.rows(); vt_i++)
    {
      d(vt_i, 0) = compute_geodist(vt_i, *HE, *e2w);
      cout << "V " << vt_i << ", geodist= " << d(vt_i, 0) << endl;
    }
    cout << "PERSO(2): geodesic distance for target vertex " << ID_VT << " -> " << d(ID_VT, 0) << endl;
    const double strip_size = 0.05;
    // The function should be 1 on each integer coordinate
    d = (d / strip_size * igl::PI).array().sin().abs().eval();
    // Compute per-vertex colors
    Eigen::MatrixXd C;
    igl::colormap(igl::COLOR_MAP_TYPE_INFERNO, d, false, C);
    // set_meshes(viewer, V1, F1);
    viewer.data().set_mesh(V1, F1);
    viewer.data().set_colors(C);
    show_source_target(viewer);
    // update the mesh (both coordinates and faces)
  }
  else if (key == '3')
  {
    cout << "Execute LIB IGL" << endl;
    // FROM LIB IGL: https://github.com/libigl/libigl/blob/master/tutorial/206_GeodesicDistance/main.cpp
    viewer.data().clear();
    Eigen::VectorXi VS, FS, VT, FT;
    // The selected vertex is the source
    VS.resize(1);
    VS << ID_VS;
    // All vertices are the targets
    VT.setLinSpaced(V1.rows(), 0, V1.rows() - 1);
    Eigen::VectorXd d;
    igl::exact_geodesic(V1, F1, VS, FS, VT, FT, d);
    cout << "LIB IGL: geodesic distance for target vertex " << ID_VT << " -> " << d(ID_VT) << endl;
    const double strip_size = 0.05;
    // The function should be 1 on each integer coordinate
    d = (d / strip_size * igl::PI).array().sin().abs().eval();
    // Compute per-vertex colors
    Eigen::MatrixXd C;
    igl::colormap(igl::COLOR_MAP_TYPE_INFERNO, d, false, C);
    // Plot the mesh
    viewer.data().set_mesh(V1, F1);
    viewer.data().set_colors(C);
    show_source_target(viewer);
  }

  return false;
}

void set_meshes(igl::opengl::glfw::Viewer &viewer, MatrixXd &V, MatrixXi &F)
{
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().set_mesh(V, F);
  // viewer.data().show_faces = true;
  viewer.data().show_lines = false;
}

void set_pc(igl::opengl::glfw::Viewer &viewer)
{
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().add_points(V1, Eigen::RowVector3d(0.3, 0.8, 0.3));
}

int main(int argc, char *argv[])
{
  char *file;
  if (argc < 2)
  {
    file = "../data/../data/gargoyle_tri.off";
  }
  else
  {
    file = argv[1];
  }

  ID_VS = 0;
  ID_VT = 5; // sphere

  //ID_VS = 0; ID_VT = 1902; // gargoyle

  // ID_VS = 1195;
  // ID_VT = 1902; // gargoyle paper

  igl::readOFF(file, V1, F1);
  igl::opengl::glfw::Viewer &viewer = VIEWER;

  HalfedgeBuilder *builder = new HalfedgeBuilder();
  HalfedgeDS he = builder->createMeshWithFaces(V1.rows(), F1);
  HE = &he;

  set_meshes(viewer, V1, F1);
  show_source_target(viewer);
  viewer.launch();
}