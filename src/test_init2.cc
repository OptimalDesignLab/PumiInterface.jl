#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>
#include <apf.h>
#include <gmi.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
//#include <apfNumbering.h>
//#include <apfShape.h>
#include <ma.h>
#include <crv.h>
#include "mpi.h"
#include <gmi_sim.h>
#include <SimUtil.h>

#include "dgSBPShape9.h"  // DGLagrange FieldShape
// prototype for constructing minimal reproducable examples, using gmi_sim


double computeValue(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3& xi, int order, apf::Vector3& x)
{
  apf::MeshElement* me = apf::createMeshElement(m, e);
  apf::mapLocalToGlobal(me, xi, x);
  apf::destroyMeshElement(me);

  return pow(x.x(), order) + pow(x.y(), order) + 1;
}


void writeField(apf::Mesh* m, apf::Field* f, int order)
{
  apf::FieldShape* fshape = apf::getShape(f);
  apf::MeshEntity* e;
  apf::Vector3 xi, x;
  apf::MeshIterator* it = m->begin(m->getDimension());
  double vals[1];

  while (( e = m->iterate(it)))
  {
    int type = m->getType(e);
    int nnodes = fshape->countNodesOn(type);
    std::cout << "element " << e << " has " << nnodes << " nodes" << std::endl;
    for (int i=0; i < nnodes; ++i)
    {
      fshape->getNodeXi(type, i, xi);
      std::cout << "xi (" << xi.x() << ", " << xi.y() << ")" << std::endl;
      vals[0] = computeValue(m, e, xi, order, x);
      apf::setComponents(f, e, i, vals);
    }
  }

  m->end(it);
}


void checkField(apf::Mesh* m, apf::Field* f, int order)
{
  apf::FieldShape* fshape = apf::getShape(f);
  apf::MeshEntity* e;
  apf::Vector3 xi, x;
  apf::MeshIterator* it = m->begin(m->getDimension());
  double vals[1];
  double val_exact;

  while (( e = m->iterate(it)))
  {
    int type = m->getType(e);
    int nnodes = fshape->countNodesOn(type);
    std::cout << "element " << e << " has " << nnodes << " nodes" << std::endl;
    for (int i=0; i < nnodes; ++i)
    {
      std::cout << "\nnode " << i << std::endl;
      fshape->getNodeXi(type, i, xi);
      std::cout << "xi (" << xi.x() << ", " << xi.y() << ")" << std::endl;
      val_exact = computeValue(m, e, xi, order, x);
      apf::getComponents(f, e, i, vals);

      double diff = std::abs(vals[0] - val_exact);
      if (diff > 1e-13)
      {
        std::cout << "point (" << x.x() << ", " << x.y() << ")" << std::endl;
        std::cout << "exact value = " << val_exact << ", field value = " << vals[0] << std::endl;
        std::cout << "diff = " << diff << std::endl;
      }
    }
  }

  m->end(it);
}



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

  std::cout << "initializing geo_sim" << std::endl;
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();

  std::cout << "loading geometric model" << std::endl;
  gmi_model* g = gmi_load(geoname);
  std::cout << "finished loading geometric model" << std::endl;
  apf::Mesh2* m = apf::loadMdsMesh(g, meshname );

  std::cout << "finished loading mesh" << std::endl;
  apf::reorderMdsMesh(m);

  //apf::writeASCIIVtkFiles("output_check", m);
  apf::writeVtkFiles("output_check", m);

  // create field
  std::cout << "creating field" << std::endl;
  int order = 3;
  apf::FieldShape* fshape = apf::getDG9SBPShape(order, m->getDimension());
  apf::Field* f = apf::createPackedField(m, "ho_field", 1, fshape);
  writeField(m, f, order);
  ma::SolutionTransfer* soltrans = ma::createFieldTransfer(f);

  std::cout << "adapting mesh" << std::endl;
  //ma::runUniformRefinement(m);
  ma::Input* in = ma::configureUniformRefine(m, 1, soltrans);
  ma::adapt(in);

  std::cout << "checking field" << std::endl;
  checkField(m, f, order);


  apf::writeVtkFiles("output_adapted", m);

  return 0;
}
     
