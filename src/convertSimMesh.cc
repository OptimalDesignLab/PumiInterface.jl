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
#include <apfShape.h>
#include "mpi.h"

// Simmetrix files
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <apfSIM.h>


// add a field with the geometric parameters of the mid-edge nodes (if present)
void addParamField(apf::Mesh* m, pParMesh smesh)
{
  apf::FieldShape* fshape = m->getShape();
  if (fshape->hasNodesIn(1))
  {
    std::cout << "adding quadratic parameter field" << std::endl;
    apf::FieldShape* fshape_edges = apf::getConstant(1);
    apf::Field* fedges = apf::createPackedField(m, "edge_params", 2, fshape_edges);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(1);

    //pMesh smesh_local = PM_mesh(smesh, 0);
    double params[2];
    while ((e = m->iterate(it)))
    {

      params[0] = 0;
      params[1] = 0;

      // get the parameters using the Simmetrix API
      pEdge sedge = (pEdge)e;
      pPoint spoint = E_point(sedge, 0);
      int gtype = E_whatInType(sedge);
      if (gtype == Gedge)
        params[0] = P_param1(spoint);
      else if (gtype == Gface)
        P_param2(spoint, &params[0], &params[1], 0);
      // else this mid-edge node is classified on a region, leave params = 0

      // save to field so it will get converted later
      apf::setComponents(fedges, e, 0, params);
    }  // end while

    m->end(it);
  }
}

void fixMatches(apf::Mesh2* m)
{
  // this seems to re-arrange the adjacency structure of the MDS mesh, not
  // really sure why it is necessary, but the regular convert.cc does it.
  if (m->hasMatching())
    apf::alignMdsMatches(m);
}


int main(int argc, char** argv)
{

  // parse arguments
  if (argc < 4 || argc > 5)
  {
    std::cerr << "Usage: " << argv[0] << " <model file> [native model] <simmetrix mesh> <scorec mesh>" << std::endl;
    return 1;
  }

  const char* geofile = argv[1];
  const char* geofile_native = NULL;
  const char* smsfile = NULL;
  const char* smbfile = NULL;

  if (argc == 4)
  {
    smsfile = argv[2];
    smbfile = argv[3];
  } else
  {
    geofile_native = argv[2];
    smsfile = argv[3];
    smbfile = argv[4];
  }

  // initialize all the libraries
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc, &argv);
  gmi_sim_start();
  gmi_register_sim();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  // load geometric model
  gmi_model* g;
  if (geofile_native)
    g = gmi_sim_load(geofile_native, geofile);
  else
    g = gmi_load(geofile);

  pGModel smodel = gmi_export_sim(g);

  // load simmetrix mesh
  double t0 = PCU_Time();
  pParMesh smesh = PM_load(smsfile, smodel, progress);
  double t1 = PCU_Time();
  std::cout << "loaded simmetrix mesh in " << t1 - t0 << " seconds" << std::endl;

  // make apfSIM mesh
  apf::Mesh* apf_sim_mesh = apf::createMesh(smesh);
  addParamField(apf_sim_mesh, smesh);

  // make MDS mesh
  t0 = PCU_Time();
  apf::Mesh2* mesh = apf::createMdsMesh(g, apf_sim_mesh);
  t1 = PCU_Time();
  std::cout << "created MDS mesh in " << t1 - t0 << " seconds" << std::endl;

  apf::printStats(mesh);

  // free simmetrix data structures
  apf::destroyMesh(apf_sim_mesh);
  M_release(smesh);

  // verify the MDS mesh
  fixMatches(mesh);
  mesh->verify();
  mesh->writeNative(smbfile);

  mesh->destroyNative();
  apf::destroyMesh(mesh);
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}
