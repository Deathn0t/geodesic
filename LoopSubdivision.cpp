/* ========================================================================= *
 *                                                                           *
 *                       Luca Castelli Aleardi                       		 *
 *           Copyright (c) 2019, Ecole Polytechnique                		 *
 *           Department of Computer Science                  				 *
 *                          All rights reserved.                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of the course material developed for		             *
 *   INF574 Digital Representation and Analysis of Shapes (2019/20)			 *
 * ========================================================================= */
#include <igl/opengl/glfw/Viewer.h>

#ifndef HALFEDGE_DS_HEADER
#define HALFEDGE_DS_HEADER
#include "HalfedgeDS.cpp"
#endif

using namespace Eigen;
using namespace std;

/**
 * @author Luca Castelli Aleardi (2019)
 */
class LoopSubdivision
{

public:
    /**
	 * Initialize the data structures
	 **/
    LoopSubdivision(MatrixXd &V_original, MatrixXi &F_original, HalfedgeDS &mesh)
    {
        he = &mesh;
        V = &V_original;
        F = &F_original;                   // NOT NEEDED if using the half-edge data structure
        int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
        int n = V_original.rows();         // number of vertices in the original mesh

        // TO BE COMPLETED (initialize arrays V1 and F1)
        nVertices = n + e;
        nFaces = 4 * F->rows();
        V1 = new MatrixXd(nVertices, 3);
        F1 = new MatrixXi(nFaces, 3);
    }

    /**
	 * Perform the subdivision of the mesh (just perform one round subdivision). <b>
	 * As result, the resulting subdivided mesh is stored in arrays 'V1' and 'F1'
	 **/
    void subdivide()
    {
        std::cout << "Performing one round subdivision" << endl;
        int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
        int n = he->sizeOfVertices();      // number of vertices in the original mesh
        int F = he->sizeOfFaces();         // number of vertices in the original mesh

        // TO BE COMPLETED
        // first step: perform a copy of the original points
        for (int i = 0; i < n; i++)
        {
            V1->row(i) = V->row(i);
        }

        // second step: compute new midpoint vertices and assign a number, between 0..e-1, to all halfedges
        int h, count;
        std::map<int, int> m;
        count = 0;
        for (int f = 0; f < F; f++)
        {
            h = he->getEdgeInFace(f);
            for (int i = 0; i < 3; i++)
            {
                if (m.find(h) == m.end())
                {
                    // not found
                    m[h] = n + count;
                    m[he->getOpposite(h)] = m[h];
                    count++;
                    V1->row(m[h]) = computeEdgePoint(h).row(0);
                }
                h = he->getNext(h);
            }
        }
        // third step: set the face/vertex incidence relations
        for (int f = 0; f < F; f++)
        {
            h = he->getEdgeInFace(f);
            for (int i = 0; i < 3; i++)
            {
                F1->row(f * 4 + i)(0) = m[h];
                F1->row(f * 4 + i)(1) = he->getTarget(h);
                F1->row(f * 4 + i)(2) = m[he->getNext(h)];
                h = he->getNext(h);
            }
            F1->row(f * 4 + 3)(0) = m[he->getPrev(h)];
            F1->row(f * 4 + 3)(1) = m[h];
            F1->row(f * 4 + 3)(2) = m[he->getNext(h)];
        }

        for (int v = 0; v < n; v++)
        {
            V1->row(v) = updateOriginalPoint(v).row(0);
        }
    }

    /**
	 * Return the number of half-edges
	 **/
    MatrixXd getVertexCoordinates()
    {
        return *V1;
    }

    /**
	 * Return the number of faces
	 **/
    MatrixXi getFaces()
    {
        return *F1;
    }

    /**
	 * Print the combinatorial information of the subdivided mesh <b>
	 * verbosity=0: print only the number of vertices and faces <b>
	 * verbosity=1: print all incidence relations
	 **/
    void print(int verbosity)
    {
        cout << "\tn=" << nVertices << ", f=" << nFaces << endl;

        if (verbosity > 0) // print all vertex coordinates and face/vertex incidence relations
        {
            for (int i = 0; i < nVertices; i++)
            {
                cout << "v" << i << ": " << V1->row(i) << endl;
            }

            std::cout << "new faces: " << nFaces << endl;
            for (int i = 0; i < nFaces; i++)
            {
                cout << "f" << i << ": " << F1->row(i) << endl;
            }
        }
    }

private:
    /**
	 * Compute the midpoint of the given half-edge 'h=(u,v)'
	 */
    MatrixXd computeEdgePoint(int h)
    {
        // TO BE COMPLETED
        MatrixXd m(1, 3);
        int v0, v1, v2, v4;
        v0 = he->getTarget(h);
        v1 = he->getTarget(he->getNext(h));
        v2 = he->getTarget(he->getOpposite(h));
        v4 = he->getTarget(he->getNext(he->getOpposite(h)));
        m.row(0) = (3. / 8.) * V->row(v0);
        m.row(0) += (3. / 8.) * V->row(v2);
        m.row(0) += (1. / 8.) * V->row(v1);
        m.row(0) += (1. / 8.) * V->row(v4);
        return m;
    }

    /**
	 * Given a vertex 'v' of the original mesh, compute and return its new coordinates
	 */
    MatrixXd updateOriginalPoint(int v)
    {
        // TO BE COMPLETED
        int h = he->getEdge(v);
        int d = vertexDegree(v);
        double a;
        MatrixXd vp = MatrixXd::Zero(1, 3);
        if (d == 3)
        {
            a = 3. / 16.;
        }
        else
        {
            a = 3. / (8. * d);
        }

        for (int j = 0; j < d; j++)
        {
            vp.row(0) = vp.row(0) + V->row(he->getTarget(he->getOpposite(h)));
            h = he->getOpposite(he->getNext(h));
        }
        vp.row(0) = a * vp.row(0);
        vp.row(0) = vp.row(0) + (1. - a * d) * V->row(v);
        return vp;
    }

    int vertexDegree(int v)
    {
        // TO BE COMPLETED
        int result = 0;
        int e = he->getEdge(v);

        int pEdge = he->getOpposite(he->getNext(e));

        while (pEdge != e)
        {
            pEdge = he->getOpposite(he->getNext(pEdge));
            result++;
        }
        return result + 1;
    }

    /** Half-edge representation of the original input mesh */
    HalfedgeDS *he;
    MatrixXd *V; // vertex coordinates of the original input mesh

    /** faces/vertex incidence relations in the original mesh */
    MatrixXi *F; // REMARK: not needed if using the half-edge data structure

    int nVertices, nFaces; // number of vertices, faces in the new subdivided mesh
    MatrixXd *V1;          // vertex coordinates of the new subdivided mesh
    MatrixXi *F1;          // faces of the new subdivided mesh
};
