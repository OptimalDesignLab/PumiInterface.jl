// this file contains the code necessary to perform the subtriangulatio of a 
// high order mesh into a low order one

#include "triangulation.h"

// given nodenum, the index of a node within an element (zero based), return the MeshEntity
// pointer for the corresponding vert on the new mesh
// typeOffsetsPerElement are the (1-based) indicies of the start index of the nodes
// of each entity type on an element
// el is the pointer to the element to which the node belongs
// entity_nodes_on is the number of nodes on verts, edges, faces
apf::MeshEntity* getVert(apf::Mesh* m, apf::MeshEntity* verts[], apf::MeshEntity* edges[], apf::MeshEntity* faces[], int typeOffsetsPerElement[], const int nodenum, uint8_t offset,  apf::MeshEntity* el, apf::Numbering* numberings[], int entity_nodes_on[])
{
  int entity_idx;  // local index of entity containing the node
  apf::MeshEntity* e;
  int entity_num;  // global number of entity
  apf::Downward down;
  int pos;  // linear address within verts, edges, faces
  int entity_node_idx;  // index of node from 0 to numnodes on all entity of this type
  int entity_node_idx_local; // index from 0 to nnodes on this entity
  int entity_node_offset_idx_local; // offset local idx
  int dim;  // dimension of entity containing the node
  if (nodenum < (typeOffsetsPerElement[1] - 1))
  {
    std::cout << "vertex node" << std::endl;
    dim = 0;
    entity_node_idx = nodenum - (typeOffsetsPerElement[dim] - 1);
    std::cout << "entity_node_idx = " << entity_node_idx << std::endl;
    entity_idx = entity_node_idx/entity_nodes_on[dim]; // integer division
    std::cout << "entity_idx = " << entity_idx << std::endl;
    entity_node_idx_local = entity_node_idx % entity_nodes_on[dim];
    std::cout << "entity_node_idx_local = " << entity_node_idx_local << std::endl;
    std::cout << "offset = " << (int)offset << std::endl;
    entity_node_offset_idx_local = abs(offset - (entity_node_idx_local + 1)) - 1;
    std::cout << "entity_node_offset_idx_local = " << entity_node_offset_idx_local << std::endl;

    // get the entity number
    m->getDownward(el, dim, down);
    e = down[entity_idx];
    entity_num = apf::getNumber(numberings[dim], e, 0, 0);
    std::cout << "entity_num = " << entity_num << std::endl;

    pos = entity_num*entity_nodes_on[dim] + entity_node_offset_idx_local;

    return verts[pos];

  } else if (nodenum < (typeOffsetsPerElement[2] - 1))
  {
    std::cout << "edge node" << std::endl;
    dim = 1;
    entity_node_idx = nodenum - (typeOffsetsPerElement[dim] - 1);
    std::cout << "entity_node_idx = " << entity_node_idx << std::endl;
    entity_idx = entity_node_idx/entity_nodes_on[dim]; // integer division
    std::cout << "entity_idx = " << entity_idx << std::endl;
    entity_node_idx_local = entity_node_idx % entity_nodes_on[dim];
    std::cout << "entity_node_idx_local = " << entity_node_idx_local << std::endl;
    std::cout << "offset = " << (int)offset << std::endl;
    entity_node_offset_idx_local = abs(offset - (entity_node_idx_local + 1)) - 1;
    std::cout << "entity_node_offset_idx_local = " << entity_node_offset_idx_local << std::endl;

    // get the entity number
    m->getDownward(el, dim, down);
    e = down[entity_idx];
    entity_num = apf::getNumber(numberings[dim], e, 0, 0);
    std::cout << "entity_num = " << entity_num << std::endl;

    pos = entity_num*entity_nodes_on[dim] + entity_node_offset_idx_local;


    
    return edges[pos];

  } else if (nodenum < (typeOffsetsPerElement[3] - 1))
  {

    std::cout << "face node" << std::endl;
    dim = 2;
    entity_node_idx = nodenum - (typeOffsetsPerElement[dim] - 1);
    std::cout << "entity_node_idx = " << entity_node_idx << std::endl;
    entity_idx = entity_node_idx/entity_nodes_on[dim]; // integer division
    std::cout << "entity_idx = " << entity_idx << std::endl;
    entity_node_idx_local = entity_node_idx % entity_nodes_on[dim];
    std::cout << "entity_node_idx_local = " << entity_node_idx_local << std::endl;
    std::cout << "offset = " << (int)offset << std::endl;
    entity_node_offset_idx_local = abs(offset - (entity_node_idx_local + 1)) - 1;
    std::cout << "entity_node_offset_idx_local = " << entity_node_offset_idx_local << std::endl;

    // get the entity number
    m->getDownward(el, dim, down);
    e = down[entity_idx];
    entity_num = apf::getNumber(numberings[dim], e, 0, 0);
    std::cout << "entity_num = " << entity_num << std::endl;

    pos = entity_num*entity_nodes_on[dim] + entity_node_offset_idx_local;
    std::cout << "pos = " << pos << std::endl;

    return faces[pos];

  } else {
    std::cerr << "Warning: in getVert,  nodenum too high, returning NULL" << std::endl;
    return NULL;
  }
}
  





apf::Mesh2* createSubMesh(apf::Mesh* m, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3])
{
// m is the existing (high order) mesh
// numtriangles is the number of triangles to break each large triangle into
// triangulation is a numtriangles x 3 array holding the indices of the nodes, 
// 1-based 
//   to use as the vertices of the sub-triangles
//   numberings is an array holding a numbering for each entity dimension

  apf::FieldShape* mshape = m->getShape();
  // get number of different types of entities
  std::size_t entity_counts[3] = {m->count(0), m->count(1), m->count(2)};

  // get number of nodes on different types of entities
  int entity_nodes_on[3] = {mshape->countNodesOn(0), mshape->countNodesOn(1), mshape->countNodesOn(2)};

  // create new mesh
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m_new = apf::makeEmptyMdsMesh(g, 2, false);

  // step 1: create all vertices of sub-triangles
   
  // allocate arrays to hold vertices of each order entity
  // these array should be in the Julia ordering, 
  // ie. getNumber(e_containing_edge) = first index
  apf::MeshEntity* verts[entity_counts[0]][entity_nodes_on[0]];
  apf::MeshEntity* edges[entity_counts[1]][entity_nodes_on[1]];
  apf::MeshEntity* faces[entity_counts[2]][entity_nodes_on[2]];
//  std::vector< std::vector<apf::MeshEntity*> > verts(entity_counts[0]);
//  std::vector< std::vector<apf::MeshEntity*> > edges(entity_counts[1]);
//  std::vector< std::vector<apf::MeshEntity*> > faces(entity_counts[2]);
//
//    apf::MeshEntity* verts_all[3] = {verts[0][0], edges[0][0], faces[0][0]};
//
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
       
     // create verts on the new mesh, store them in arrays in the same
     // order as the julia entities are stored
      while ( (e = m->iterate(it)) )
      {
       int idx = apf::getNumber(numberings[dim], e, 0, 0);
       m->getPoint(e, 0, coords);
       verts[idx][0] = m_new->createVert(0); // should be 2d array
       m_new->setPoint(verts[idx][0], 0, coords);

       std::cout << "creating vertex for vertex " << idx << " at coordinates " << coords << " pointer = " << verts[idx][0] << std::endl;
      }

      dim = 1;  
      it = m->begin(dim);  // iterate over vertices
      while ( (e = m->iterate(it)) )
      {
        int idx = apf::getNumber(numberings[dim], e, 0, 0);
        for (int i = 0; i < entity_nodes_on[dim]; i++)
        {
          m->getPoint(e, i, coords);
          edges[idx][i] = m_new->createVert(0); // should be 2d array
          m_new->setPoint(edges[idx][i], 0, coords);

          std::cout << "creating edge vertex for edge " << idx << ", node " << i << " at coordinates " << coords << " pointer = " << edges[idx][i]<< std::endl;
        }
      }

      dim = 2;  
      it = m->begin(dim);  // iterate over vertices
      while ( (e = m->iterate(it)) )
      {
        int idx = apf::getNumber(numberings[dim], e, 0, 0);
        for (int i = 0; i < entity_nodes_on[dim]; i++)
        {
          m->getPoint(e, i, coords);

          faces[idx][i] = m_new->createVert(0); // should be 2d array
          m_new->setPoint(faces[idx][i], 0, coords);

          std::cout << "creating edge vertex for face " << idx << ", node " << i << " at coordinates " << coords << " pointer = " << faces[idx][i] <<  std::endl;
        }
      }



    
  // step 2: create elements from vertices
  // loop over elements on m, create element on m_new
  // here we use the edge orientation information
  dim = 2;
  int elnum; // holds the element number
  int node;  // original node number
  int pos;  // linear address for elementNodeOffset
  const int nnodes_per_el = typeOffsetsPerElement[dim+1] - 1;
  apf::MeshEntity* el_verts[3]; // holds vertices of element to be created
  it = m->begin(dim);  // iterate over elements
  while ( (e = m->iterate(it)) )
  {
    elnum = apf::getNumber(numberings[dim], e, 0, 0); // zero base index
    std::cout << "\n\nsubtriangulating element " << elnum << std::endl;
    for (int i=0; i < numtriangles; ++i)  // loop over all subtriangles
    {
      std::cout << "\ncreating sub triangle " << i << std::endl;
      // get 3 vertices
      for (int j=0; j < 3; ++j)
      {
        node = triangulation[i][j] - 1; // one based index
        std::cout << "node = " << node << std::endl;
        // calculate linear offset for elementNode offsets
        pos = node + elnum*nnodes_per_el;
        std::cout << "offset pos = " << pos << std::endl;
        int offset_j = elementNodeOffsets[pos];
        
//        std::cout << "offset = " << offset_j << std::endl;
//        newnode = abs(elementNodeOffsets[pos] - node) - 1;
//        std::cout << "newnode = " << newnode << std::endl;
//        el_verts[j] = getVert(m, verts, edges, faces, typeOffsetsPerElement, newnode, e, numberings, entity_nodes_on );

        el_verts[j] = getVert(m, &verts[0][0], &edges[0][0], &faces[0][0], typeOffsetsPerElement, node, offset_j, e, numberings, entity_nodes_on );
        std::cout << "  el_vert " << j << " = " << el_verts[j] << std::endl;
      }
     
      // build the element 
      apf::buildElement(m_new, 0, apf::Mesh::TRIANGLE, el_verts);

    }  // end loop over subtriangles

  } // end loop over elements
   

  // build, verify  mesh
  std::cout << "deriving model" << std::endl;
  apf::deriveMdsModel(m_new);
  std::cout << "finished deriving model" << std::endl;
  m_new->acceptChanges();
  std::cout << "accepted changes" << std::endl;
  m_new->verify();
  std::cout << "verified" << std::endl;

  // write visualization file
  apf::writeVtkFiles("newmesh_linear", m_new);



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
