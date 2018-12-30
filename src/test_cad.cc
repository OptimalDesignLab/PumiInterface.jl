#include <iostream>
#include <fstream>
#include <cassert>
#include <string.h>
#include <apf.h>
#include <gmi.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include "mpi.h"
//#include "pumiInterface_config.h"
#define HAVE_SIMMETRIX // avoid CMake dependency
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif

/* The purpose of this program is to test the accuracy of the derivatives
 * computes by the CAD model by comparing against finite difference.
 * This is done by defining a scalar-valued function for each coordinate
 * of each node of the coordinate field.  The function can be written as
 * f(x(xi)) where x is the Cartesian coordinate of the coordinate dof and
 * xi is the vector of parametric coordinates, which can be used to
 * compute x. On the interior of the domain, the xi coordinates are equal to
 * the xyz coordinates, for points on the boundary, the xi coordinates are the
 * parametric coordinates defined by the CAD system.
 * The finite difference for a given coordinate dof can be
 * computed as df/dxi_1 = (f( x(xi_1 + h) ) - f(x(xi)) )/h, where the
 * subscript _1 denotes the first element of the xi vector.  The CAD-based
 * derivative can be computed as df/dxi_1 = df/dx * dx/dxi_1.
 * It appears the xi -> x calculation (apf::Mesh::snapToModel) is accurate to
 * about 10^-9, however the derivative is only accurate to 10^-4 in the typical
 * case, and can have up to 20% error in some cases.
 */

// define a simple function f(x) for each coordinate of each node
// i is the coordinate (0 < i < m->getDimension())
double f_x(double x[3], int i)
{
  return 5*x[i] + 1;
}

double f_x(apf::Vector3 _x, int i)
{
  double x[3];
  _x.toArray(x);

  return f_x(x, i);
}

// compute the same function as f( x(xi) )
double f_xi(apf::Mesh* m, apf::MeshEntity* e, double xi[3], int i)
{
  apf::Vector3 x;
  apf::Vector3 _xi;
  apf::ModelEntity* me = m->toModel(e);
  int me_dim = m->getModelType(me);

  if (me_dim == m->getDimension())
  {
    for (int d=0; d < 3; ++d)
      x[d] = xi[d];
  } else
  {
    _xi.fromArray(xi);
    m->snapToModel(me, _xi, x);
  }

  return f_x(x, i);
}

// derivative of f
double df_dx(double x[3], int i)
{
  return 5;
}

double df_dx(apf::Vector3 _x, int i)
{
  double x[3];
  _x.toArray(x);
  return df_dx(x, i);
}


// returns number of parametric coordinates a given MeshEntity has
int countXiDofs(apf::Mesh* m, apf::MeshEntity* e)
{
  apf::ModelEntity* me = m->toModel(e);
  return m->getModelType(me);
}

void printArray(double* x, int n)
{
  std::cout << "[ ";
  for (int i=0; i < n; ++i)
    std::cout << x[i] << ", ";

  std::cout << "]" << std::endl;
}

void printArray(apf::Vector3 x)
{
  double _x[3];
  x.toArray(_x);
  printArray(_x, 3);
}


// Given df/dx, compute dx/dxi using the CAD provided derivative.
// node: specifies the node of the given MeshEntity (=0 for 1st and 2nd order
// coordinate fields).
// dx: contains df/dx (because there is a different function f for
// every coordinate of every node, only one value can be non-zero in this
// test case).
// dxi: overwritten with df/dxi
void dXTodXi(apf::Mesh* m, apf::MeshEntity* e, int node, double dx[3], double dxi[3])
{
  apf::ModelEntity* me = m->toModel(e);
  int me_dim = m->getModelType(me);
  assert(me_dim > 0); // geometric vertices have 0 dofs

  if (me_dim == m->getDimension())
  {
    // for interior points, use xyz coordinates
    for (int i=0; i < m->getDimension(); ++i)
      dxi[i] = dx[i];

  } else  // if not interior and not vertex, this point is on a bounding surface
  {
    // compute df/dxi = df/dx * dx/dxi + df/dy * dy/dxi
    apf::Vector3 x;  // xyz coordinates
    apf::Vector3 newx;
    apf::Vector3 xi;  // parametric coordinates
    m->getPoint(e, node, x);
    m->getClosestPoint(me, x, newx, xi);

    apf::Vector3 dx_dxi1, dx_dxi2;
    m->getFirstDerivative(me, xi, dx_dxi1, dx_dxi2);

    for (int d=0; d < 3; ++d)
      dxi[d] = 0;

    for (int d=0; d < m->getDimension(); ++d)
      dxi[0] += dx[d]*dx_dxi1[d];

    if (me_dim == 2)
      for (int d=0; d < m->getDimension(); ++d)
        dxi[1] += dx[d]*dx_dxi2[d];
  }

}


// Gets both the xyz and parametric coordinates of a given node.
void getNodeCoords(apf::Mesh* m, apf::MeshEntity* e, int node, double x[3],
                    double xi[3])
{
  apf::Vector3 _x, _newx, _xi;
  m->getPoint(e, node, _x);
  apf::ModelEntity* me = m->toModel(e);
  int me_dim = m->getModelType(me);

  // for geometrics regions, use xyz coordinates, for bounding surfaces, use
  // CAD parametric coordinates
  if (me_dim == m->getDimension())
  {
    for (int d=0; d < 3; ++d)
      xi[d] = x[d];

    _x.toArray(x);
  } else
  {
    apf::Vector3 x_snap;
    m->getClosestPoint(me, _x, _newx, _xi);
    m->snapToModel(me, _xi, x_snap);
    _xi.toArray(xi);
    std::cout << "x_orig = ";
    printArray(_x);
    std::cout << "closest x = ";
    printArray(_newx);
    std::cout << "x_snap = ";
    printArray(x_snap);
    std::cout << "xi = ";
    printArray(_xi);

    if (xi[0] < 1e-5)
    {
      double h = 1e-4;
      apf::Vector3 xp, _xi2;
      for (int d=0; d < 3; ++d)
        _xi2[d] = _xi[d];

      // print out coordinate at increments
      for (int i=0; i < 5; ++i)
      {
        _xi2[0] += h;
        m->snapToModel(me, _xi2, xp);
        std::cout << "point " << i << " coords = ";
        printArray(xp);
      }


      // check derivative at other end of loop
      double rng[2];
      m->getPeriodicRange(me, 0, rng);
      _xi2[0] = rng[1];
      m->snapToModel(me, _xi2, xp);
      std::cout << "at other end of loop, xi = ";
      printArray(_xi2);
      std::cout << "point = ";
      printArray(xp);

      apf::Vector3 t0, t1, t2, t3, t4, t5;
      m->getFirstDerivative(me, _xi, t0, t1);
      m->getFirstDerivative(me, _xi2, t2, t3);
      _xi2[0] = rng[0];
      m->getFirstDerivative(me, _xi2, t4, t5);

      std::cout << "derivative at xi of xyz point: ";
      printArray(t0);
      std::cout << "derivative at xi = ximax: ";
      printArray(t2);
      std::cout << "derivative at xi = 0: ";
      printArray(t4);

    }

    x_snap.toArray(x);  // try to return consistent x, xi
  }
}


// test CAD derivatives against finite difference
void test_fd(apf::Mesh* m)
{

  std::cout.precision(10);
  double h = 1e-6;  // finite difference step size
  double x[3] = {0, 0, 0};
  double xi[3] = {0, 0, 0};
  double dx[3] = {0, 0, 0};
  double dxi[3] = {0, 0, 0};
  double dxi_fd[3] = {0, 0, 0};

  double maxdiff = 0;
  double absdiff[3];

  apf::FieldShape* fshape = m->getShape();
  apf::Field* ferr = apf::createPackedField(m, "derivative error", m->getDimension(), fshape);
  for (int dim=0; dim <= m->getDimension(); ++dim)
  {
    if (!fshape->hasNodesIn(dim))
      continue;

    std::cout << "dimension " << dim << std::endl;
    apf::MeshIterator* it = m->begin(dim);
    apf::MeshEntity* e;

    while ( (e = m->iterate(it)) )
    {

      int ndof = countXiDofs(m, e);
      if (ndof == 0)
      {
        // write zeros to field
        for (int d=0; d < m->getDimension(); ++d)
          absdiff[d] = 0;

        for (int j=0; j < fshape->countNodesOn(m->getType(e)); ++j)
          apf::setComponents(ferr, e, j, absdiff);

        continue;
      }
      // The code works correctly for interior points, but the problem is
      // the boundary points, so do only those
      //if (ndof == m->getDimension())
      //  continue;

      for (int j=0; j < fshape->countNodesOn(m->getType(e)); ++j)
      {
        getNodeCoords(m, e, j, x, xi);


        for (int k=0; k < m->getDimension(); ++k)  // loop over xyz coords
        {
          // Note: f_x is actually a different function for every k value
          
          absdiff[k] = 0;

          std::cout << "k = " << k << std::endl;
          double val1 = f_xi(m, e, xi, k);  // value at initial point

          // loop over xi coordinates
          for (int d=0; d < ndof; ++d)
          {
            // finite difference
            xi[d] += h;
            double val2 = f_xi(m, e, xi, k);
            xi[d] -= h;
            dxi_fd[d] = (val2 - val1)/h;
          }
             
          // compute df/dxi = df/dx * dx/dxi
          dx[k] = df_dx(x, k);
          dXTodXi(m, e, j, dx, dxi);

          // compare results
          for (int d=0; d < ndof; ++d)
          {
            std::cout << "  d = " << d << ", df/dxi = " << dxi[d];
            std::cout << ", df/dxi fd = " << dxi_fd[d];
            std::cout << ", diff = " << dxi[d] - dxi_fd[d] << std::endl;

            double adiff = std::fabs(dxi[d] - dxi_fd[d]);
            if (adiff > maxdiff)
            {
              std::cout << "  updating maxdiff" << std::endl;
              maxdiff = adiff;
            }

            if (adiff > absdiff[k])
              absdiff[k] = adiff;
          }


          // zero things out
          dx[k] = 0;

        }  // end k

        apf::setComponents(ferr, e, j, absdiff);
      } // end j
    }  // end while

    m->end(it);
  }  // end dim

  std::cout << "maxdiff = " << maxdiff << std::endl;

}  // function


int main (int argc, char** argv)
{

  if (argc < 2 || argc > 3)
  {
    std::cerr << "Usage: test_init mesh_name.smb [geometry_file]" << std::endl;
    return 1;
  }

  const char* meshname = argv[1];
  char geoname[512];
  if (argc == 3)
    strcpy(geoname, argv[2]);
  else
    strcpy(geoname, ".null");

  std::cout << "Entered init\n" << std::endl;
  std::cout << "geoname = " << geoname << std::endl;
  
  MPI_Init(0,NULL);  // initilize MPI
  PCU_Comm_Init();   // initilize PUMI's communication

  // load mesh using null geometry
  gmi_register_null();
  gmi_register_mesh();
#ifdef HAVE_SIMMETRIX
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif

  std::cout << "loading geometric model" << std::endl;
  gmi_model* g = gmi_load(geoname);
  std::cout << "finished loading geometric model" << std::endl;
  apf::Mesh2* m = apf::loadMdsMesh(g, meshname );

  std::cout << "finished loading mesh" << std::endl;
  apf::reorderMdsMesh(m);

  // run the test
  test_fd(m);
//  apf::writeASCIIVtkFiles("output_check", m);
  apf::writeVtkFiles("output_check", m);

  std::cout << "finished writing paraview files" << std::endl;

  return 0;
}
     
