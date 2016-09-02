#include <iostream>
#include <fstream>
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include "mpi.h"
//#include "dgSBPShape1.h"
//#include "apfSBPShape.h"

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

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  int peer_rank = 1 - myrank;  // only 2 processes

  // load mesh using null geometry
  gmi_register_null();
  gmi_register_mesh();
  std::cout << "loading geometric model" << std::endl;
  gmi_model* g = gmi_load("parallel2.dmg");
  std::cout << "finished loading geometric model" << std::endl;
  // using the mesh vortex3_1.smb works fine
  apf::Mesh2* m = apf::loadMdsMesh(g,"parallel2.smb" );

  std::cout << "finished loading mesh" << std::endl;
  apf::reorderMdsMesh(m);

//  apf::FieldShape* fshape = apf::getSBPShape(1);
  apf::FieldShape* fshape = apf::getLagrange(1);
  apf::changeMeshShape(m, fshape, true);
  std::cout << "finished changing mesh shape" << std::endl;

  apf::Sharing* shr = getSharing(m);
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(1);
  std::ofstream f;
  char fname[256];
  sprintf(fname, "fout_%d.dat", myrank);
  f.open(fname);

  apf::MeshIterator* it = m->begin(3);
  apf::MeshEntity* e;
  apf::MeshEntity* verts[4];
  apf::MeshEntity* face_verts[3];
  apf::MeshEntity* remote_face_verts[3];
  apf::Vector3 point;

  while ( (e = m->iterate(it)) )
  {
    if ( e == 0x000000000000009f )
    {
      m->getDownward(e, 0, verts);

      // extract the face verts and print their coordinates
      for (int j=0; j < 3; ++j)
      {
        face_verts[j] = verts[j];
        m->getPoint(verts[j], 0, point);
        std::cout << "vert " << j << " coordinates = (" << point.x() << ", " << point.y() << ", " << point.z() << std::endl;

        // get remote pointer




      // print their coordinates





    }

  }

  for (int i=0; i <= 19; ++i)  // element 
  {

  }

  f.close();

//  apf::writeASCIIVtkFiles("output_check", m);
  apf::writeVtkFiles("output_check", m);

  std::cout << "finished writing paraview files" << std::endl;

//  printModelClassification(m);

  return 0;
}
     
