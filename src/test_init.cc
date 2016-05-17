#include <iostream>
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include "dgSBPShape1.h"

/*
void printModelClassification(apf::Mesh * m)
{
  apf::MeshEntity* e;
  apf::ModelEntity* me;
  int model_dim;
  int model_tag;
  int cnt = 0;
  char const* entity_names[4] = { "vertex ", "edge ", "face ", "region "};

  for (int i = 0; i < 4; ++ i)
  {
    apf::MeshIterator* it = m->begin(i);
    cnt = 0;
    while (e = m->iterate(it))
    {
      me = m->toModel(e);
      model_dim = m->getModelType(me);
      model_tag = m->getModelTag(me);
      std::cout << entity_names[i] << cnt << " is classified on model entity ";
      std::cout << model_tag << " of dimension " << model_dim << std::endl;
      ++cnt;
    }
    std::cout << std::endl;
  }
        

}  // end function
*/
int main ()
{


  std::cout << "Entered init\n" << std::endl;
  
  MPI_Init(0,NULL);  // initilize MPI
  PCU_Comm_Init();   // initilize PUMI's communication

  // load mesh using null geometry
  gmi_register_null();
  std::cout << "loading null geometric model" << std::endl;
  gmi_model* g = gmi_load(".null");
  std::cout << "finished loading geometric model" << std::endl;
  // using the mesh vortex3_1.smb works fine
  apf::Mesh2* m = apf::loadMdsMesh(g,"tri8l.smb" );

  std::cout << "finished loading mesh" << std::endl;
  apf::FieldShape* fshape = apf::getLagrange(1);
  apf::changeMeshShape(m, fshape, false);
  std::cout << "finished changing mesh shape" << std::endl;

  // create a Numbering with 3 nodes, all classified on faces, with two 
  // components on each node
  int ncomp = 2;
  fshape = apf::getDG1SBPShape(1);
  int nnodes_per_el = fshape->countNodesOn(apf::Mesh::TRIANGLE);
  apf::Numbering* n = apf::createNumbering(m, "dof numbers", fshape, ncomp);
  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;

  // populate the numbering sequentially
  int val = 0;
  while ( (e = m->iterate(it)) )
  {
    for (int node = 0; node < nnodes_per_el; ++node)
    {
      for (int comp = 0; comp < ncomp; ++comp)
      {
        apf::number(n, e, node, comp, val);
        ++val;
      }
    }
  }


//  apf::writeASCIIVtkFiles("output_check", m);
  apf::writeVtkFiles("output_check", m);

  std::cout << "finished writing paraview files" << std::endl;

//  printModelClassification(m);

  return 0;
}
     
