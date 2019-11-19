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

void set_meshes(igl::opengl::glfw::Viewer &viewer, MatrixXd &V, MatrixXi &F)
{
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().set_mesh(V, F);
  viewer.append_mesh();
  viewer.data(0).set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3));
}

void set_pc(igl::opengl::glfw::Viewer &viewer)
{
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().add_points(V1, Eigen::RowVector3d(0.3, 0.8, 0.3));
}





float distance(Vector3d &v1, Vector3d &v2)
{
    return sqrt(pow(v2(2) - v1(2), 2) +  pow(v2(1) - v1(1), 2) +  pow(v2(0) - v1(0), 2));
}

void add_window_Q(std::queue<Window> &Q, Vector3d &vs, Vector3d &v0, Vector3d &v1,int edge_id)
{
  double b0 = 0.;
  double b1 = distance(v0, v1);

  double d0 = distance(vs, v0);
  double d1 = distance(vs, v1);

  double sigma = 0;
  int dir = 0;
  Q.push(Window(b0, b1, d0, d1, sigma, dir, edge_id,v0,v1));
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
void init_Q(HalfedgeDS &he, int id_vs, MatrixXd &V, std::queue<Window> &Q)
{
  // TODO
  Q = std::queue<Window>();

  Vector3d vs, v_init, v_b0, v_b1;
  int incidentEdge, nextEdge;

  vs = V.row(id_vs);
  incidentEdge = he.getEdge(id_vs);
	nextEdge = he.getNext(incidentEdge);
	v_init = V.row(he.getTarget(nextEdge));

  nextEdge = he.getNext(nextEdge);
	v_b1 = V.row(he.getTarget(nextEdge));
  add_window_Q(Q, vs, v_init, v_b1, nextEdge);

	while (v_init != v_b1)
	{
        v_b0 = v_b1;

        nextEdge = he.getNext(he.getOpposite(he.getNext(nextEdge)));
				v_b1 = V.row(he.getTarget(nextEdge));

        add_window_Q(Q, vs, v_b0, v_b1, nextEdge);
  }
  // std::cout<<Q.size()<<std::endl;
}

void propagate_window(MatrixXd &V,HalfedgeDS &he,Window &w)
{

   // On calcule nos angles de rotation entre notre ancien repere et le nouveau

    int nextEdge = he.getNext(he.getOpposite(w.get_edge_id()));
    Vector3d u = (w.get_v0() - w.get_v1()).transpose();
    Vector3d v = w.get_v0() - V.row(he.getTarget(nextEdge)).transpose();
    std::cout<<"u "<<u<<std::endl;
    std::cout<<"v "<<v<<std::endl;
    double theta = acos(u.dot(v)/(u.norm()*v.norm()));
    std::cout<<"theta "<<theta<<std::endl;
    
    // Calcul des nouveaux points b0, b1 (old window) apres rotation.

    MatrixXd R(2,2);
    R(0,0) = cos(theta); R(0,1) = -sin(theta);
    R(1,0) = sin(theta); R(1,1) = cos(theta);
    std::cout<<"R "<<R<<std::endl;

    Vector2d new_v0 = R*Vector2d(w.get_b0(),0);
    Vector2d new_v1 = R*Vector2d(w.get_b1(),0);
    
    std::cout<<"new_vo "<<new_v0<<std::endl;
    std::cout<<"new_v1 "<<new_v1<<std::endl;
    
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
 * Compute exact geodesics fron a source vertex to very over vertices of the mesh.
 * input:
 *    F: faces of the mesh where to compute geodesics pathes
 *    V: vertexes of the mesh where to compute geodesics pathes
 *    id_vs: the source vertex id
 * output:
 *    ...: something with all windows computed to apply backtracing?
 */
void exact_geodesics(HalfedgeDS &he,MatrixXd &V, MatrixXi &F, int id_vs)
{

  // initialize the queue Q with a window for each edge adjacent
  // to source: source_v_t
  std::queue<Window> Q; // Q.push(..); Q.front(); Q.pop();
  init_Q(he,id_vs, V, Q);

  Window cur_w;         // current window during iterations

  while (!Q.empty())
  {
    // select and remove a Window from Q
    cur_w = Q.front();
    Q.pop();
    
    propagate_window(V,he,cur_w);
    return;
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

}

void example_1()
{
  igl::readOFF("../data/star.off", V1, F1);
  int vs = 0;
  HalfedgeBuilder *builder = new HalfedgeBuilder();
  HalfedgeDS he = (builder->createMeshWithFaces(V1.rows(), F1));

  exact_geodesics(he, V1, F1, vs);

  igl::opengl::glfw::Viewer viewer;
  set_meshes(viewer,V1,F1);
  viewer.launch();
}

int main(int argc, char *argv[])
{
  example_1();
}
