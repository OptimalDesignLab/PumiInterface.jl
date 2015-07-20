#include <iostream>
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>


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
  apf::Mesh2* m = apf::loadMdsMesh(g, ".smb");

  apf::writeVtkFiles("output_check", m);

  return 0;
}
     
