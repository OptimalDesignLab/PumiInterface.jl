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
#include "pumiInterface_config.h"
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif

//#include "dgSBPShape1.h"
//#include "apfSBPShape.h"


void printModelEntities(gmi_model* g)
{
//  double p[2];
//  double x[3];
  double r[2];
  for (int dim=0; dim <= 3; dim++)
  {
    std::cout << "dimension " << dim << " geometric entities" << std::endl;
    gmi_iter* it = gmi_begin(g, dim);
    gmi_ent* ge;
    while ( (ge = gmi_next(g, it)) )
    {
      int ge_tag = gmi_tag(g, ge);
      std::cout << "  geometry entity " << ge << " has tag " << ge_tag << std::endl;
      if (dim == 1)
      {
        gmi_range(g, ge, 0, r);
        std::cout << "    parametric range = " << r[0] << ", " << r[1] << std::endl;
        bool isPeriodic = gmi_periodic(g, ge, 0);
        std::cout << "    isPeriodic = " << isPeriodic << std::endl;
      }
    }
    gmi_end(g, it);
  }

}


int main (int argc, char** argv)
{

  if (argc != 2)
  {
    std::cerr << "Usage: printGeoEntities geometry_file" << std::endl;
    return 1;
  }

  const char* geoname = argv[1];


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

  printModelEntities(g);

  return 0;
}
     
