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
#include <ma.h>
#include "mpi.h"
#include "test_meshadapt.h"
//#include "dgSBPShape1.h"
//#include "apfSBPShape.h"

// print the names of all Numbering and Fields on the mesh
void printAssociatedData(apf::Mesh* m)
{

  std::cout << "Numberings:" << std::endl;
  // Numberings
  int nnumberings = m->countNumberings();
  for (int i=0; i < nnumberings; ++i)
  {
    apf::Numbering* num_i = m->getNumbering(i);
    std::cout << "  Numbering " << i << " has name " << apf::getName(num_i) << std::endl;
  }

  std::cout << "Fields:" << std::endl;
  int nfields = m->countFields();
  for (int i=0; i < nfields; ++i)
  {
    apf::Field* field_i = m->getField(i);
    std::cout << "  Field " << i << " has name " << apf::getName(field_i) << std::endl;

  }

}

int main (int argc, char** argv)
{

  if (argc != 2)
  {
    std::cerr << "Usage: test_init mesh_name.smb" << std::endl;
    return 1;
  }

  std::cout << "Entered init\n" << std::endl;
  
  MPI_Init(0,NULL);  // initilize MPI
  PCU_Comm_Init();   // initilize PUMI's communication

  // load mesh using null geometry
  gmi_register_null();
  gmi_register_mesh();
  std::cout << "loading geometric model" << std::endl;
  gmi_model* g = gmi_load(".null");
  std::cout << "finished loading geometric model" << std::endl;
  apf::Mesh2* m = apf::loadMdsMesh(g, argv[1] );

  std::cout << "finished loading mesh" << std::endl;
  apf::reorderMdsMesh(m);

  apf::FieldShape* fshape_orig = m->getShape();
  std::cout << "fieldshape_orig name = " << fshape_orig->getName() << std::endl;

  // see if coordinates were moved into tags
  apf::MeshTag* coords_tag = m->findTag("coordinates_ver");
  std::cout << "coords_tag = " << coords_tag << std::endl;
  apf::DynamicArray<apf::MeshTag*> tags;
  m->getTags(tags);

  for (unsigned int i=0; i < tags.getSize(); ++i)
  {
    std::cout << "tag " << i << " name = " << m->getTagName(tags[i]) << std::endl;
  }

  // force the coordinates into the field
  int order_orig = m->getShape()->getOrder();
  if (coords_tag != 0)
  {
    std::cout << "performing initial shape change" << std::endl;
    apf::changeMeshShape(m, apf::getLagrange(order_orig), true);
  }

  apf::writeASCIIVtkFiles("output_check", m);
//  apf::writeVtkFiles("output_check", m);

  std::cout << "finished writing paraview files" << std::endl;

  
  // create the numberings
  apf::Numbering* numberings[4];
  char name_buff[256];
  for (int i = 0; i <= m->getDimension(); ++i)
  {
    sprintf(name_buff, "entity%d", i);
    numberings[i] = apf::numberOwnedNodes(m, name_buff, apf::getConstant(i));
  }

  // write desired element size to field and then create a SolutionTransfer
  // for it
  apf::Field* f = apf::createFieldOn(m, "desired_size_field", apf::SCALAR);
  apf::Field* fsol = apf::createPackedField(m, "solution_field_interp", 4, m->getShape());

  for (int dim = 0; dim < m->getDimension(); ++dim)
  {
    if (!m->getShape()->hasNodesIn(dim))
      continue;

    apf::MeshIterator* it = m->begin(dim);
    apf::MeshEntity* e;
    apf::Vector3 coords;
    double comps[4];
    while ( (e = m->iterate(it)) )
    {
      m->getPoint(e, 0, coords);
      // go from element size 2/10 at inner radius to 1/10 at outer radius
      const double h_inner = 2.0/10.0;
      const double h_outer = 1.0/100.0;
      const double r_slope = (h_outer - h_inner)/2.0;
      const double r_in = 1;
      double r = std::sqrt(coords.x()*coords.x() + coords.y()*coords.y());
      double val = 2.0/10.0 + r_slope*(r - r_in);
      apf::setScalar(f, e, 0, val);
      
      comps[0] = r;
      comps[1] = r+1;
      comps[2] = r+2;
      comps[3] = r+3;

      apf::setComponents(fsol, e, 0, comps);

    }
  }  // end for loop


  ma::SolutionTransfers soltrans;
  soltrans.add(ma::createFieldTransfer(f));
  soltrans.add(ma::createFieldTransfer(fsol));

//  ma::SolutionTransfer* soltrans = ma::createFieldTransfer(f);

  // create the default IsotropicFunction
  IsotropicFunctionJ* isofunc = new IsotropicFunctionJ(m, numberings, f);

  std::cout << "isofunc address = " << isofunc << std::endl;

  // print fields before
  std::cout << "before mesh adaptation:" << std::endl;
  printAssociatedData(m);

  apf::writeASCIIVtkFiles("output_pre", m);
  for (int i=0; i < 1; ++i)
  {
    ma::adapt(m, isofunc, &soltrans);
    sprintf(name_buff, "mesh_iter%d", i);
    apf::writeVtkFiles(name_buff, m);
  }

  // print fields before
  std::cout << "after mesh adaptation:" << std::endl;
  printAssociatedData(m);
  apf::writeASCIIVtkFiles("output_post", m);

  // cleanup
  delete isofunc;



  return 0;
}
     
