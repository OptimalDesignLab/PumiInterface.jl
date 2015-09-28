// this file contains the code necessary to perform the subtriangulatio of a 
// high order mesh into a low order one

#include "triangulation.h"

apf::Mesh2* createSubMesh(apf::Mesh* m, const int numtriangles, const int triangulation[][3], apf::Numbering* numberings[3])
{
// m is the existing (high order) mesh
// numtriangles is the number of triangles to break each large triangle into
// triangulation is a numtriangles x 3 array holding the indices of the nodes
//   to use as the vertices of the sub-triangles
//   numberings is an array holding a numbering for each entity dimension

  apf::FieldShape* mshape = m->getShape();
  // get number of different types of entities
  int entity_counts[3] = {m->count(0), m->count(1), m->count(2)};

  // get number of n odes on different types of entities
  int entity_nodes_on[3] = {mshape->countNodesOn(0), mshape->countNodesOn(1), mshape->countNodesOn(2)};

  // create new mesh
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m_new = apf::makeEmptyMdsMesh(g, 2, false);

  // step 1: create all vertices of sub-triangles
   
  // allocate arrays to hold vertices of each order entity
  apf::MeshEntity* verts[entity_counts[0]][entity_nodes_on[0]];
  apf::MeshEntity* edges[entity_counts[1]][entity_nodes_on[1]];
  apf::MeshEntity* faces[entity_counts[2]][entity_nodes_on[2]];

  // step 2: create elements from vertices
  
  int dim = 0;  // indicate vertices
  apf::MeshIterator* it = m->begin(dim);  // iterate over vertices
  apf::MeshEntity* e;  // uninitilized pointer?
  apf::Vector3 coords;
  while ( (e = m->iterate(it)) )
  {
   int idx = apf::getNumber(numberings[dim], e, 0, 0);
   m->getPoint(e, 0, coords);
   verts[idx][0] = m_new->createVert(0);
   m_new->setPoint(verts[idx][0], 0, coords);
  }


  // step 3: transfer fields



}


// function to determine the orientation of an edge in the high order mesh
// this assumes the convention used in PdePumiInterface that the the element
// with centroid x coordinate 
int getOrientation(apf::Mesh* m, apf::MeshEntity* e)
{ 
  int e_type = m->getType(e);

  // declare some static variables
  static apf::Up up_faces;
  static apf::Downward verts1;
  static apf::Downward verts2;
  static double coords1[3][3];
  static double coords2[3][3];

  if (m->typeDimenions[e_type] == 0 || m->typeDimensions[e_type] == 2)
    return 0;
  
  // else this must be an edge

  m->getUp(e, up_faces);
  apf::MeshElement* e1 = e[0];
  apf::MeshElement* e2 = e[1];

  // calculate centroid of two elements
  // get coordinates first, then calculate centroid
  m->getDownward(e1, 0, verts1);
  m->getDownward(e2, 0, verts2);

  for (int i = 0; i < 3; ++i)
  {


