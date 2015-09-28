// this file contains the code necessary to perform the subtriangulatio of a 
// high order mesh into a low order one

#include "triangulation.h"

apf::Mesh2* createSubMesh(apf::Mesh* m, int, numtriangles, int triangulation[][3], apf::Numbering* numberings[3])
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
  
  apf::MeshIterator* it = m->begin(0);  // iterate over vertices
  apf::MeshEntity* e;  // uninitilized pointer?
//  while (ap


  // step 3: transfer fields



}
