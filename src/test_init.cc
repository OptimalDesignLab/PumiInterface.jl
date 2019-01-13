#include <iostream>
#include <fstream>
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
#define HAVE_SIMMETRIX
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif


void printParams(apf::Mesh* m)
{

  apf::FieldShape* fshape = m->getShape();
  std::cout << "Printing CAD parametric coordinates for " << fshape->getName() << std::endl;

  apf::MeshEntity* e;
  apf::MeshIterator* it;
  apf::Vector3 p; p[2] = 0;
  apf::ModelEntity* me;
  for (int dim=0; dim <= m->getDimension(); ++dim)
  {
    if (!fshape->hasNodesIn(dim))
      continue;

    std::cout << "Dimension " << dim << std::endl;
    it = m->begin(dim);
    while ( (e = m->iterate(it)) )
    {
      // assume 1 node per entity
      m->getParam(e, p);
      me = m->toModel(e);
      int tag = m->getModelTag(me);
      int me_dim = m->getModelType(me);
      std::cout << "  entity " << e << " classified on model entity " << tag << " of dimension " << me_dim << " has parametric coords = " << p.x();
      std::cout << ", " << p.y() << ", " << p.z() << std::endl;

    }
    m->end(it);

  }
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
#ifdef HAVE_SIMMETRIX
  std::cout << "initializing geo_sim" << std::endl;
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#else
  std::cout << "not initializing geo_sim" << std::endl;
#endif
  std::cout << "loading geometric model" << std::endl;
  gmi_model* g = gmi_load(geoname);
  std::cout << "finished loading geometric model" << std::endl;
  apf::Mesh2* m = apf::loadMdsMesh(g, meshname );

  std::cout << "finished loading mesh" << std::endl;
  apf::reorderMdsMesh(m);

  apf::FieldShape* fshape_orig = m->getShape();
  std::cout << "fieldshape_orig name = " << fshape_orig->getName() << std::endl;

  // see if coordinates were moved into tags
  apf::MeshTag* coords_tag = m->findTag("coordinates_ver");
  std::cout << "coords_tag = " << coords_tag << std::endl;

  // force the coordinates into the field
  if (coords_tag != 0)
  {
    int order_orig = m->getShape()->getOrder();
    std::cout << "performing initial shape change" << std::endl;
    apf::changeMeshShape(m, apf::getLagrange(order_orig), false);
  }

  printParams(m);
  //apf::writeASCIIVtkFiles("output_check", m);
  apf::writeVtkFiles("output_check", m);

  std::cout << "finished writing paraview files" << std::endl;

  return 0;
}
     
