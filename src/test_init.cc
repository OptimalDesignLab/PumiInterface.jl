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

void printRemoteInfo(apf::Mesh* m)
{

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // iostream
  std::fstream fs;
  char fname[50];
  sprintf(fname, "output_%d.dat", myrank);
  fs.open(fname, std::fstream::out);

  // mesh data
  apf::Sharing* shr = apf::getSharing(m);
  apf::Copies copies;  // used with getRemotes
  apf::CopyArray shares; // used with Sharing
  int dim = m->getDimension();
  apf::MeshIterator* it = m->begin(dim-1);
  apf::MeshEntity* e;

  apf::Parts parts;
  apf::getPeers(m, dim-1, parts);

  for (apf::Parts::iterator it = parts.begin(); it != parts.end(); ++it)
  {
    fs << "peer " << *it << std::endl;
  }



  int facenum = 0;
  while ( (e = m->iterate(it)) )
  {
    fs << "\nfacenum = " << facenum << std::endl;

    copies.clear();
    m->getRemotes(e, copies);

    fs << "printing remotes:" << std::endl;
    for (apf::Copies::iterator it = copies.begin(); it != copies.end(); ++it)
    {
      fs << "part num = " << it->first << ", entity = " << it->second << std::endl;
    }

    fs << "printing shares:" << std::endl;
    shares.setSize(0);
    shr->getCopies(e, shares);

    for (int i=0; i < int(shares.getSize()); ++i)
    {
      fs << "part num = " << shares[i].peer << ", entity = " << shares[i].entity << std::endl;
    }

  }

  m->end(it);
  fs.close();

} // function printRemoteInfo
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

void printCoordinates(apf::Mesh* m)
{
  apf::MeshEntity* e;
  apf::Vector3 coords;
  apf::FieldShape* fshape = m->getShape();
  apf::MeshIterator* it;

  for (int dim=0; dim <= m->getDimension(); dim++)
  {
    if (!fshape->hasNodesIn(dim))
      continue;

    std::cout << "Dimension " << dim << ":" << std::endl;
    it = m->begin(dim);
    int i = 0;
    while ( (e = m->iterate(it)) )
    {
      m->getPoint(e, 0, coords);
      std::cout << "  dimension " << dim << " entity " << i << " coords = " << coords.x() << ", " << coords.y() << coords.z() << std::endl;
      
    }
  }

} // function printCoordinates

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
  // using the mesh vortex3_1.smb works fine
  apf::Mesh2* m = apf::loadMdsMesh(g, argv[1] );

  std::cout << "finished loading mesh" << std::endl;
  apf::reorderMdsMesh(m);

  apf::FieldShape* fshape_orig = m->getShape();
  std::cout << "fieldshape_orig name = " << fshape_orig->getName() << std::endl;
//  int order_orig = fshape_orig->getOrder();


  // see if coordinates were moved into tags
  apf::MeshTag* coords_tag = m->findTag("coordinates_ver");
  std::cout << "coords_tag = " << coords_tag << std::endl;
  apf::DynamicArray<apf::MeshTag*> tags;
  m->getTags(tags);

  for (unsigned int i=0; i < tags.getSize(); ++i)
  {
    std::cout << "tag " << i << " name = " << m->getTagName(tags[i]) << std::endl;
  }

//  printCoordinates(m);
  printRemoteInfo(m);
/*
  if (coords_tag != 0)
  {
    std::cout << "performing initial shape change" << std::endl;
    apf::changeMeshShape(m, apf::getLagrange(order_orig), false);
  }
*/
//  apf::FieldShape* fshape = apf::getSBPShape(1);
//  apf::FieldShape* fshape = apf::getLagrange(1);
//  apf::changeMeshShape(m, fshape, true);
//  std::cout << "finished changing mesh shape" << std::endl;

//  apf::writeASCIIVtkFiles("output_check", m);
  apf::writeVtkFiles("output_check", m);

  std::cout << "finished writing paraview files" << std::endl;

//  printModelClassification(m);

  return 0;
}
     
