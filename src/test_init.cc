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
#endif

//#include "dgSBPShape1.h"
//#include "apfSBPShape.h"

/*
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
*/

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
      }
    }
    gmi_end(g, it);
  }

}


void printMidpointErrors(apf::Mesh* m)
{
  gmi_model* g = m->getModel();

  const int dim = 0;
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  gmi_ent* ge;
  apf::Downward verts;
  apf::Vector3 x1, x2;
  double midpoint[3], newpoint[3], newp[2];
  while ( (e = m->iterate(it)) )
  {
    ge = (gmi_ent*)m->toModel(e);
    if (dim == 0)
    {
      m->getPoint(e, 0, x1);
      x1.toArray(midpoint);
    } else
    {
      m->getDownward(e, 0, verts);
      m->getPoint(verts[0], 0, x1);
      m->getPoint(verts[1], 0, x2);

      midpoint[0] = (x1.x() + x2.x())/2;
      midpoint[1] =( x1.y() + x2.y())/2;
      midpoint[2] =( x1.z() + x2.z())/2;
    }

    gmi_closest_point(g, ge, midpoint, newpoint, newp);

    std::cout << "distance from midpoint to geometry is " << newpoint[0] - midpoint[0] << ", " << newpoint[1] - midpoint[1] << ", " << std::endl;
  }

  m->end(it);
}


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
    while ((e = m->iterate(it)))
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
      std::cout << "  dimension " << dim << " entity " << i << " coords = " << coords.x() << ", " << coords.y() << ", " << coords.z() << std::endl;
      
    }
  }

} // function printCoordinates

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
  apf::DynamicArray<apf::MeshTag*> tags;
  m->getTags(tags);

  for (unsigned int i=0; i < tags.getSize(); ++i)
  {
    std::cout << "tag " << i << " name = " << m->getTagName(tags[i]) << std::endl;
  }

  std::cout << "before changing shape, coords = " << std::endl;
  printCoordinates(m);
//  printRemoteInfo(m);

  // force the coordinates into the field
  if (coords_tag != 0)
  {
    int order_orig = m->getShape()->getOrder();
    std::cout << "performing initial shape change" << std::endl;
    apf::changeMeshShape(m, apf::getLagrange(order_orig), false);
  }

  // print coordinates again
  std::cout << "after changing shape, coords = " << std::endl;
  printCoordinates(m);

  std::cout << "printing model classification" << std::endl;
  printModelClassification(m);
  printModelEntities(g);
  printMidpointErrors(m);
  apf::writeASCIIVtkFiles("output_check", m);
//  apf::writeVtkFiles("output_check", m);

  std::cout << "finished writing paraview files" << std::endl;

  return 0;
}
     
