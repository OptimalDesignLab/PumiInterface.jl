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
//#include <apfNumbering.h>
//#include <apfShape.h>
#include "mpi.h"
#include "pumiInterface_config.h"
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif
#include <ma.h>

int main (int argc, char** argv)
{

  MPI_Init(0,NULL);  // initilize MPI
  PCU_Comm_Init();   // initilize PUMI's communication

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (argc < 3 || argc > 5)
  {
    if (myrank == 0)
      std::cerr << "Usage: " << argv[0] << " mesh_name.smb geometry_file outname [number of refinements]" << std::endl;
    return 1;
  }

  const char* meshname = argv[1];
  const char* geoname = argv[2];
  const char* outname = argv[3];
  int nrefines = 1;
  if (argc == 5)
  {
    nrefines = atoi(argv[4]);
    if (nrefines < 1)
    {
      if (myrank == 0)
        std::cerr << "Error: Number of refines requestions: " << nrefines << " < 1" << std::endl;
      return 1;
    }
  }

  gmi_register_null();
  gmi_register_mesh();
#ifdef HAVE_SIMMETRIX
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif

  double t1 = PCU_Time();
  gmi_model* g = gmi_load(geoname);
  apf::Mesh2* m = apf::loadMdsMesh(g, meshname);
  double t2 = PCU_Time();
  if (myrank == 0)
    std::cout << "loaded mesh and geometry in " << t2 - t1 << " seconds" << std::endl;

  double t3 = PCU_Time();
  ma::runUniformRefinement(m, nrefines);
  double t4 = PCU_Time();
  if (myrank == 0)
    std::cout << "mesh uniformly refined in " << t4 - t3 << " seconds" << std::endl;

  apf::writeVtkFiles("output_check", m);
  m->writeNative(outname);

  m->destroyNative();
  apf::destroyMesh(m);

#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
#endif


  return 0;
}
     
