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
  std::size_t entity_counts[3] = {m->count(0), m->count(1), m->count(2)};

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
//  std::vector< std::vector<apf::MeshEntity*> > verts(entity_counts[0]);
//  std::vector< std::vector<apf::MeshEntity*> > edges(entity_counts[1]);
//  std::vector< std::vector<apf::MeshEntity*> > faces(entity_counts[2]);
    apf::MeshEntity* verts_all[3] = {verts[0][0], edges[0][0], faces[0][0]};
//  std::vector<std::array*> verts_all(3);
//  std::vector<std::vector*> verts_all(3);
//  verts_all[0] = &verts;
//  verts_all[1] = &edges;
//  verts_all[2] = &faces;

//  std::vector<apf::MeshEntity&> verts_all(3);
//  verts_all[0] = verts;
//  verts_all[1] = edges;
//  verts_all[2] = faces;
//  apf::MeshEntity* (*verts_all)[3] = {&verts[0], &edges[0], &faces[0]};  
  int dim = 0;  // indicate vertices
  apf::MeshIterator* it = m->begin(dim);  // iterate over vertices
  apf::MeshEntity* e;  // uninitilized pointer?
  apf::Vector3 coords;
     // because C is a stupid language with arbitrary and inconsistent rules 
     // about arrays
       
//      apf::MeshEntity* arr_i[entity_counts[dim]][entity_nodes_on[dim]] = verts_all[dim];

/*
      while ( (e = m->iterate(it)) )
      {
       int idx = apf::getNumber(numberings[dim], e, 0, 0);
       m->getPoint(e, 0, coords);
       arr_i[idx][0] = m_new->createVert(0); // should be 2d array
       m_new->setPoint(verts[idx][0], 0, coords);
      }

      dim = 1;  
      apf::MeshIterator* it = m->begin(dim);  // iterate over vertices
      while ( (e = m->iterate(it)) )
      {
       int idx = apf::getNumber(numberings[dim], e, 0, 0);
       m->getPoint(e, 0, coords);
       arr_i[idx][0] = m_new->createVert(0); // should be 2d array
       m_new->setPoint(verts[idx][0], 0, coords);
      }
*/


  // step 2: create elements from vertices
  // loop over elements on m, create element on m_new
  // here we use the edge orientation information
  

  // step 3: transfer fields
  // loop over elements on m, transfer values to m_new
  // this results in repeatedly copying some information, but
  // that is the price to pay for having a simple loop over element
  // which enables the code to know the relationship between 
  // elements on m and m_new
  // here we also use the element orientation information



  return m_new;
}


// function to determine the orientation of an edge in the high order mesh
// this assumes the convention used in PdePumiInterface that the the element
// with the smallest  centroid x coordinate is elementL, and therefore "owns"
// the edge, ie. node 0 of the edge is the first node from the perspective
// of elementL
int getOrientation(apf::Mesh* m, apf::MeshEntity* e)
{ 
  int e_type = m->getType(e);

  // declare some static variables
  static apf::Up up_faces;
  static apf::Downward verts1;
  static apf::Downward verts2;
//  static double coords1[3][3];
//  static double coords2[3][3];
  static apf::Vector3 coords;
  static apf::Vector3 centroid1(0, 0, 0);
  static apf::Vector3 centroid2(0, 0, 0);


  if (m->typeDimension[e_type] == 0 || m->typeDimension[e_type] == 2)
  {
    return 0;
  }
  // else this must be an edge

  m->getUp(e, up_faces);
  apf::MeshEntity* e1;
  apf::MeshEntity* e2;
  e1 = up_faces.e[0];
  e2 = up_faces.e[1];
  // calculate centroid of two elements
  // get coordinates first, then calculate centroid
  m->getDownward(e1, 0, verts1);
  m->getDownward(e2, 0, verts2);

  for (int i = 0; i < 3; ++i)
  {
    m->getPoint(verts1[i], 0, coords); 
    centroid1 += coords;

    m->getPoint(verts2[i], 0, coords); 
    centroid2 += coords;
  }

  centroid1 = centroid1/3.0;
  centroid2 = centroid2/3.0;
/*
  // if x cordinates are different enough
  if (abs(centroid1[0] - centroid2[0]) > 1e-10)
  {
    if (centroid1[0] < centroid2[1])
    {
*/
  return 42;  // need to change this
}
