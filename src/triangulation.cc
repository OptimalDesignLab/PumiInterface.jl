// this file contains the code necessary to perform the subtriangulatio of a 
// high order mesh into a low order one

#include "triangulation.h"
#include "funcs1.h"

namespace tri {
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
//    std::cout << "vertex node" << std::endl;
    dim = 0;
    entity_node_idx = nodenum - (typeOffsetsPerElement[dim] - 1);
//    std::cout << "entity_node_idx = " << entity_node_idx << std::endl;
    entity_idx = entity_node_idx/entity_nodes_on[dim]; // integer division
//    std::cout << "entity_idx = " << entity_idx << std::endl;
    entity_node_idx_local = entity_node_idx % entity_nodes_on[dim];
//    std::cout << "entity_node_idx_local = " << entity_node_idx_local << std::endl;
//    std::cout << "offset = " << (int)offset << std::endl;
    entity_node_offset_idx_local = abs(offset - (entity_node_idx_local + 1)) - 1;
//    std::cout << "entity_node_offset_idx_local = " << entity_node_offset_idx_local << std::endl;

    // get the entity number
    m->getDownward(el, dim, down);
    e = down[entity_idx];
    entity_num = apf::getNumber(numberings[dim], e, 0, 0);
//    std::cout << "entity_num = " << entity_num << std::endl;

    pos = entity_num*entity_nodes_on[dim] + entity_node_offset_idx_local;

    return verts[pos];

  } else if (nodenum < (typeOffsetsPerElement[2] - 1))
  {
//    std::cout << "edge node" << std::endl;
    dim = 1;
    entity_node_idx = nodenum - (typeOffsetsPerElement[dim] - 1);
//    std::cout << "entity_node_idx = " << entity_node_idx << std::endl;
    entity_idx = entity_node_idx/entity_nodes_on[dim]; // integer division
//    std::cout << "entity_idx = " << entity_idx << std::endl;
    entity_node_idx_local = entity_node_idx % entity_nodes_on[dim];
//    std::cout << "entity_node_idx_local = " << entity_node_idx_local << std::endl;
//    std::cout << "offset = " << (int)offset << std::endl;
    entity_node_offset_idx_local = abs(offset - (entity_node_idx_local + 1)) - 1;
//    std::cout << "entity_node_offset_idx_local = " << entity_node_offset_idx_local << std::endl;

    // get the entity number
    m->getDownward(el, dim, down);
    e = down[entity_idx];
    entity_num = apf::getNumber(numberings[dim], e, 0, 0);
//    std::cout << "entity_num = " << entity_num << std::endl;

    pos = entity_num*entity_nodes_on[dim] + entity_node_offset_idx_local;


    
    return edges[pos];

  } else if (nodenum < (typeOffsetsPerElement[3] - 1))
  {

//    std::cout << "face node" << std::endl;
    dim = 2;
    entity_node_idx = nodenum - (typeOffsetsPerElement[dim] - 1);
//    std::cout << "entity_node_idx = " << entity_node_idx << std::endl;
    entity_idx = entity_node_idx/entity_nodes_on[dim]; // integer division
//    std::cout << "entity_idx = " << entity_idx << std::endl;
    entity_node_idx_local = entity_node_idx % entity_nodes_on[dim];
//    std::cout << "entity_node_idx_local = " << entity_node_idx_local << std::endl;
//    std::cout << "offset = " << (int)offset << std::endl;
    entity_node_offset_idx_local = abs(offset - (entity_node_idx_local + 1)) - 1;
//    std::cout << "entity_node_offset_idx_local = " << entity_node_offset_idx_local << std::endl;

    // get the entity number
    m->getDownward(el, dim, down);
    e = down[entity_idx];
    entity_num = apf::getNumber(numberings[dim], e, 0, 0);
//    std::cout << "entity_num = " << entity_num << std::endl;

    pos = entity_num*entity_nodes_on[dim] + entity_node_offset_idx_local;
//    std::cout << "pos = " << pos << std::endl;

    return faces[pos];

  } else {
    std::cerr << "Warning: in getVert,  nodenum too high, returning NULL" << std::endl;
    return NULL;
  }
}
 

// copy all of the Numbering values at a particular node from an old Numbering to a new one.
int copyNumberingNode(apf::Numbering* n_old, apf::Numbering* n_new, apf::MeshEntity* e_old, apf::MeshEntity* e_new, const int node_old, const int node_new, int ncomp)
{
  apf::Mesh* m = apf::getMesh(n_old);
  apf::FieldShape* nshape = apf::getShape(n_old);
  int type = m->getType(e_old);
  int dim = m->typeDimension[type];
//  std::cout << "    dimension of old entity = " << dim << std::endl;
//  std::cout << "type of old entity = " << type << std::endl;
//  std::cout << "has nodes in this dimension = " << mshape->hasNodesIn(dim) << std::endl;
//  bool isnumbered = apf::isNumbered(n_old, e_old, node_old, 0);
//  std::cout << "is numbered = " << isnumbered << std::endl;
  bool flag = false;
  int val;
  for (int i = 0; i < ncomp; ++i)
  {
    bool node_numbered = apf::isNumbered(n_old, e_old, node_old, i);
    bool has_nodes = nshape->hasNodesIn(dim);
//    std::cout << "    node_numbered = " << node_numbered << std::endl;
//    std::cout << "    has_nodes = " << has_nodes << std::endl;
    if (node_numbered && has_nodes)
    {
//      std::cout << "    copying value" << std::endl;
      val = apf::getNumber(n_old, e_old, node_old, i);
      apf::number(n_new, e_new, node_new, i, val);
      flag = true;
    } else
    {
//      std::cout << "skipping node" << std::endl;
      flag = false;
    }
  }

  if (flag)  { return 1;}
  else { return 0;}
}

// copy all of the field values at a particular node from an old Numbering to a new one.
void copyFieldNode(apf::Field* field_old, apf::Field* field_new, apf::MeshEntity* e_old, apf::MeshEntity* e_new, const int node_old, const int node_new, int ncomp, double vals[])
{

      apf::getComponents(field_old, e_old, node_old, vals);
/*      std::cout << "    retrieved components, vals = ";
      for (int i = 0; i < ncomp; ++i)
      {
        std::cout << vals[i] << ", ";
      }
      std::cout << std::endl;
*/
      apf::setComponents(field_new, e_new, node_new, vals);
//      std::cout << "    wrote components" << std::endl;
}




void printFields(apf::Mesh* m)
{
  int numfields = m->countFields();
  std::cout << "\n\nthis mesh has " << numfields << " fields" << std::endl;
  for (int i = 0; i < numfields; ++i)
  {
    apf::Field* field_i = m->getField(i);
    const char* fieldname_i = apf::getName(field_i);
    std::cout << "field " << i << " is named " << fieldname_i << std::endl;
  }

  int numnumberings = m->countNumberings();
  std::cout << "\nthis mesh has " << numnumberings << " numberings" << std::endl;
  for (int i = 0; i < numnumberings; ++i)
  {
    apf::Numbering* numbering_i = m->getNumbering(i);
    const char* name_i = apf::getName(numbering_i);
    std::cout << "numbering " << i << " is named " << name_i << std::endl;
    apf::Field* field_i = apf::getField(numbering_i);
    std::cout << "field_i* = " << field_i << std::endl;
    if (field_i != NULL)
    {
       const char* fieldname_i = apf::getName(field_i);
      std::cout << "  and has underlying field named " << fieldname_i << std::endl;
    }
  }

  std::cout << "finished printing Numbering names" << std::endl;
}  // end function



// copy a numbering defined on the elements of the old mesh to the elements
// of the new mesh
void transferNumberingToElements(apf::Mesh* m, apf::Mesh* m_new, const int numtriangles, apf::Numbering* numbering_old)
{

  apf::MeshIterator* itold = m->begin(2);
  apf::MeshIterator* itnew = m_new->begin(2);
  apf::FieldShape* fshape_new = apf::getConstant(2);
  int numcomp = apf::countComponents(numbering_old);
  const char* name_old = apf::getName(numbering_old);
//  std::cout << "name_old = " << name_old << std::endl;
  std::size_t name_length = strlen(name_old) + 1; // string lenth, including \0
//  std::cout << "name_length = " << name_length << std::endl;
  char name_new[name_length + 4]; // add space for "_old"
  strcpy(name_new, name_old);
//  std::cout << "after first strcat, name_new = " << name_new << std::endl;
  strcat(name_new, "_old");
//  std::cout << "name_new = " << name_new << std::endl;
  apf::Numbering* numbering_new = apf::createNumbering(m_new, name_new, fshape_new, numcomp);

  apf::MeshEntity* e;
  apf::MeshEntity* e_new;
  int val;
  while ( (e = m->iterate(itold)) )  // loop over old mesh elements
  {
      for (int i = 0; i < numtriangles; ++i)
      {
        for (int n = 0; n < numcomp; ++n)
        {
          val = apf::getNumber(numbering_old, e, 0, n);
          e_new = m_new->iterate(itnew);
          apf::number(numbering_new, e_new, 0, n, val);
        }
      }
  }


}  // end function


void transferNumberings(apf::Mesh* m, apf::Mesh* m_new, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3])
{

//  std::cout << "old mesh pointer = " << m << std::endl;
//  std::cout << "new mesh pointer = " << m_new << std::endl;
  // compute some quantities
  int nnodes_per_el = typeOffsetsPerElement[3] - 1;
  apf::FieldShape* mshape = m->getShape();
  // get number of different types of entities
//  std::size_t entity_counts[3] = {m->count(0), m->count(1), m->count(2)};

  // get number of nodes on different types of entities
  int entity_nodes_on[3] = {mshape->countNodesOn(0), mshape->countNodesOn(1), mshape->countNodesOn(2)};

  const int dimensionTypes[3] = {apf::Mesh::VERTEX, apf::Mesh::EDGE, apf::Mesh::TRIANGLE};

  int node_start[3][3]; // dim x max entities per dim (so 3 edges for triangle)
  int idx = 0; // current index, 

  for (int dim = 0; dim < 3; ++dim)
  {
    // get number of nodes on entities of this type 
    int num_nodes_d = mshape->countNodesOn(dimensionTypes[dim]);
    // loop over entities in this dimension
//    std::cout << "dim = " << dim << ", adjacentCount = " << m->adjacentCount[apf::Mesh::TRIANGLE][dim] << std::endl;
    for (int i = 0; i < m->adjacentCount[apf::Mesh::TRIANGLE][dim]; ++i)
    { 
//      std::cout << "idx = " << idx << std::endl;
      node_start[dim][i] = idx;
//      std::cout << "node_start[" << dim << "][" << i << "] = " << node_start[dim][i] << std::endl;
      idx = idx + num_nodes_d;
    }
  }



  // calculate the starting node index for every entity on the element


  // step 3: transfer Numberings
  // loop over elements on m, transfer values to m_new
  // this results in repeatedly copying some information, but
  // that is the price to pay for having a simple loop over element
  // which enables the code to know the relationship between 
  // elements on m and m_new
  // here we also use the element orientation information

  // any constant fieldshape gets mappted to getConstant(0) because all nodes
  // on the old mesh reside on a vertex of the new mesh, the only question
  // is how many components the numbering has
  // the undocumented function apf::countComponenets(Numbing* n); does it
  // Similarly for fields, the fieldshape is mapped to the equivalent linear
  // field
  // unfortunately there is no way to know what the function to get the 
  // equivalent linear fieldshape is other than parsing the FieldShape
  // name string and having a table of function pointers here
  // for now we will just stick with linear Lagrange
  // the number of components ofa  field can be found with
  // apf::countComponents(Field* f)



  // precalculate tables needed to lookup a given node on the new mesh
  int subelements[nnodes_per_el]; // hold index (0-based) of subtriangles
  int subelement_verts[nnodes_per_el]; // holds vertex index of subtriangle
  int el = 0;  // store found element
  int vert = 0; // store found node
  bool breaknow = false;  // flag for breaking second loopa
//  std::cout << "Precalculating sub-triangle lookup tables" << std::endl;
  for (int i = 0; i < nnodes_per_el; ++i)  // search for node i
  {
    // loop through the triangulation
    for ( int j = 0; j < numtriangles; ++j)  // loop over subtriangles
    {
      for ( int k = 0; k < 3; ++k)  // loop over nodes on the triangle
      {
        if ( (triangulation[j][k] - 1) == i)
        {
          el = j;
          vert = k;
//          std::cout << "node " << i << " is part of sub-triangle " << el << " node " << vert << std::endl;
          breaknow = true;
          break;
        }
      }

      if (breaknow)  // if inner loop breaks, break outer one too
      {
        breaknow = false;
        break;
      }
    }

    // break statements will end up here
//    std::cout << "node " << i << " is part of sub-triangle " << el << " node " << vert << std::endl;
    subelements[i] = el;
    subelement_verts[i] = vert;
  }

//  std::cout << "Finished precalculating lookup tables" << std::endl;

   // copy numberings
  int numnumberings = m->countNumberings();

  // store all the subtriangles corresponding to an old mesh element
  apf::MeshEntity* subtriangles[numtriangles];
  int nodes_on[3];
  apf::MeshEntity* node_entity;  // entity actually containing the node
  int nodenum; // node within element index
  int elnum;
  int new_elnum;
  int new_vertnum;
  int entity_node_idx; // index of node among all nodes on entities of same dimension
  int entity_idx;  // index of node containing entity in downward adjacencies
  int node_idx; // node index, pre offset
  int offset_node_idx; // node index, post offset
  int offset_i;
  int pos;
  apf::MeshEntity* e;
  apf::Downward down; // hold downward adjacent verts of old mesh
  apf::Downward down2; // hold downwarda adjacent verts of new mesh
  apf::MeshEntity* new_node_entity;
//  int num_entities_per_dim[3] = { 3, 3, 1};
//apf::MeshEntity* getEntityFromNode(apf::Mesh* m, int typeOffsetsPerElement[], const int nodenum,  apf::MeshEntity* el, int entity_nodes_on[])

  for (int i = 0; i < numnumberings; ++i)
  {
//    std::cout << "Copying Numbering " << i << std::endl;

    apf::Numbering* numbering_i = m->getNumbering(i);
    int numcomp = apf::countComponents(numbering_i);
    const char* name_i = apf::getName(numbering_i);
//    std::cout << "Numbering name is " << name_i << std::endl;
    apf::FieldShape* numbering_i_shape = apf::getShape(numbering_i);
//    const char* shape_name = numbering_i_shape->getName();
//    std::cout << "Numbering Fieldshape name = " << shape_name << std::endl;
//    int num_numbered_nodes = apf::countNodes(numbering_i);
/*
    for (int s_dim=0; s_dim < 3; ++s_dim)
    {
      std::cout << "Numbering Fieldshape has nodes in dimension " << s_dim << " = " << numbering_i_shape->hasNodesIn(s_dim) << std::endl;
    }
*/

    for (int j = 0; j < 3; ++j)
      nodes_on[j] = numbering_i_shape->countNodesOn(j);

    int numbered_count = 0;
//    std::cout << "Number of nodes with assigned numbers = " << num_numbered_nodes << std::endl;
    apf::Field* field_i = apf::getField(numbering_i);
    if (field_i != NULL)
    {
      std::cerr << "Warning:  Numbering with an underling Field detected.  Numbering will not be copied" << std::endl;
      continue;
    }

    // if Numbering is the element numbering, copy it specially too
    if ( strcmp(name_i, "faceNums") == 0 || strcmp(name_i, "coloring") == 0)
    {
//      std::cout << "Copying numbering " << name_i << " to elements of new mesh" << std::endl;
      transferNumberingToElements(m, m_new, numtriangles, numbering_i);
      continue;
    }

    apf::FieldShape* fshape_new = triDG::getNewShape(numbering_i_shape);

    // create new field and populate it
    apf::Numbering* n_new = apf::createNumbering(m_new, name_i, fshape_new, numcomp);

    if (nodes_on[0] == 1 && nodes_on[1] == 0 && nodes_on[2] == 0)
    {
//      std::cout << "copying vertices" << std::endl;
      triDG::transferVertices(m, m_new, numtriangles, triangulation, elementNodeOffsets, typeOffsetsPerElement, numbering_i, n_new);
      triDG::completeNumbering(m_new, n_new);
      continue;
    }

    if (nodes_on[0] == 0 && nodes_on[1] == 1 && nodes_on[2] == 0)
    {
//      std::cout << "copying edges" << std::endl;
      triDG::transferEdges(m, mshape, m_new, numtriangles, triangulation, elementNodeOffsets, typeOffsetsPerElement, numbering_i, n_new);
      triDG::completeNumbering(m_new, n_new);
      continue;
    }

//    std::cout << "copying everything" << std::endl;
    apf::MeshIterator* itold = m->begin(2);
    apf::MeshIterator* itnew = m_new->begin(2);
    while ( (e = m->iterate(itold)) ) // loop over elements in old mesh
    {

      elnum = getNumber(numberings[2], e, 0, 0);

//      std::cout << "  processing old mesh element number " << elnum << std::endl;
      // get all the subtriangles of the current element
      for (int i = 0; i < numtriangles; ++i)
      {
        subtriangles[i] = m_new->iterate(itnew);
//        std::cout << "  subtriangle " << i << " pointer = " << subtriangles[i] << std::endl;
      }

      nodenum = 0;  //reset nodenum to beginning of element

      for (int dim = 0; dim < 3; ++dim)  // loop over elements
      {

//        std::cout << "  looping over dimension " << dim << std::endl;
        m->getDownward(e, dim, down);
        
//        for (int tmp = 0; tmp < 12; ++tmp)
//        {
//          std::cout << "down[" << tmp << "] = " << down[tmp] << std::endl;
//        }

        //TODO: change to loop over entities of current dim, then 
        //      loop over nodes on each entity
        //      at start of entity loop, set nodenum to right value
        //      from precompuated values
          // loop over entities on this dimension
//          int nume = m -> adjacentCount[apf::Mesh::TRIANGLE][dim];
//          std::cout << "  number of entities on this dimension = " << nume << std::endl;
          for (int entity_d = 0; entity_d < m -> adjacentCount[apf::Mesh::TRIANGLE][dim]; ++entity_d)
          {
//            std::cout << "  on entity " << entity_d << std::endl;
            // get the starting node number for this entity
            nodenum = node_start[dim][entity_d];

            // get the type of the current entity

            // loop over nodes on this entity, in the fieldshape of
            // the numbering, not the mesh
//            int numnodes_e = numbering_i_shape->countNodesOn(dimensionTypes[dim]);
//            std::cout << "number of nodes on this entity = " << numnodes_e << std::endl;
            for (int node_e = 0; node_e < numbering_i_shape->countNodesOn(dimensionTypes[dim]); ++node_e)
            {

//              std::cout << "    on node " << node_e << std::endl;

    //          std::cout << "    processing node " << nodenum << std::endl;
              // get entity on old mesh containing the node
              entity_node_idx = nodenum - (typeOffsetsPerElement[dim] - 1);
              entity_idx = entity_node_idx/entity_nodes_on[dim]; // integer division
    //          std::cout << "entity_idx = " << entity_idx << std::endl;
              node_idx = entity_node_idx % entity_nodes_on[dim];
              node_entity = down[entity_idx]; // get the mesh entity pointer
    //          std::cout << "node_entity = " << node_entity << std::endl;

              // offset
    //          std::cout << "    calculating offset" << std::endl;
              pos = nodenum + elnum*nnodes_per_el;
              offset_i = elementNodeOffsets[pos]; 
              offset_node_idx = abs(offset_i - (node_idx + 1)) - 1;

              // get the entity on the new mesh containing the node
    //          std::cout << "    retrieving new mesh entity" << std::endl;
              new_elnum = subelements[nodenum];
              new_vertnum = subelement_verts[nodenum];
    //          std::cout << "    new_elnum = " << new_elnum << std::endl;
    //          std::cout << "    new_vertnum = " << new_vertnum << std::endl;
              m_new->getDownward(subtriangles[new_elnum], 0, down2);
              new_node_entity = down2[new_vertnum];
    //          std::cout << "    new_node_entity = " << new_node_entity << std::endl;

    //          std::cout << "    copying values" << std::endl;
              numbered_count += copyNumberingNode(numbering_i, n_new, node_entity, new_node_entity, offset_node_idx, 0, numcomp);

              ++nodenum;
              }  // end loop over nodes on this entity
           } // end loop over entities on this dimension
         } // end loop over dimensions

    } // end loop over old mesh elements
   
//    int num_numbered_nodes_new = apf::countNodes(n_new);
//    std::cout << "Actually copied " << num_numbered_nodes_new << " nodes" << std::endl; 
/*  // this test doesn't work because numbering of an entity dimension
 *      gets mapped to all the nodes on that entity dimension
    if (num_numbered_nodes != num_numbered_nodes_new)
    {
      std::cerr << "Warning: incorrect transfering of Numbering" << std::endl;
      std::cerr << "num_numbered_nodes = " << num_numbered_nodes << ", num_numered_nodes_new = " << num_numbered_nodes_new << std::endl;
    }
*/
  } // end loop over numberings

//  std::cout << "printing new mesh fields and numberings" << std::endl;
//  printFields(m_new);
}  // end function

} // end namespace tri


// this function transferes the specified field to the new mesh
void transferField(apf::Mesh* m, apf::Mesh* m_new, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3], apf::Field* field_old, apf::Field* field_new)
{

  // compute some quantities
  int nnodes_per_el = typeOffsetsPerElement[3] - 1;
  apf::FieldShape* mshape = m->getShape();
  // get number of different types of entities
//  std::size_t entity_counts[3] = {m->count(0), m->count(1), m->count(2)};

  // get number of nodes on different types of entities
  int entity_nodes_on[3] = {mshape->countNodesOn(0), mshape->countNodesOn(1), mshape->countNodesOn(2)};


  // step 3: transfer Fields
  // loop over elements on m, transfer values to m_new
  // this results in repeatedly copying some information, but
  // that is the price to pay for having a simple loop over element
  // which enables the code to know the relationship between 
  // elements on m and m_new
  // here we also use the element orientation information

  // any constant fieldshape gets mappted to getConstant(0) because all nodes
  // on the old mesh reside on a vertex of the new mesh, the only question
  // is how many components the numbering has
  // the undocumented function apf::countComponenets(Numbing* n); does it
  // Similarly for fields, the fieldshape is mapped to the equivalent linear
  // field
  // unfortunately there is no way to know what the function to get the 
  // equivalent linear fieldshape is other than parsing the FieldShape
  // name string and having a table of function pointers here
  // for now we will just stick with linear Lagrange
  // the number of components ofa  field can be found with
  // apf::countComponents(Field* f)



  // precalculate tables needed to lookup a given node on the new mesh
  int subelements[nnodes_per_el]; // hold index (0-based) of subtriangles
  int subelement_verts[nnodes_per_el]; // holds vertex index of subtriangle
  int el = 0;  // store found element
  int vert = 0; // store found node
  bool breaknow = false;  // flag for breaking second loopa
//  std::cout << "Precalculating sub-triangle lookup tables" << std::endl;
  for (int i = 0; i < nnodes_per_el; ++i)  // search for node i
  {
    // loop through the triangulation
    for ( int j = 0; j < numtriangles; ++j)  // loop over subtriangles
    {
      for ( int k = 0; k < 3; ++k)  // loop over nodes on the triangle
      {
        if ( (triangulation[j][k] - 1) == i)
        {
          el = j;
          vert = k;
//          std::cout << "node " << i << " is part of sub-triangle " << el << " node " << vert << std::endl;
          breaknow = true;
          break;
        }
      }

      if (breaknow)  // if inner loop breaks, break outer one too
      {
        breaknow = false;
        break;
      }
    }

    // break statements will end up here
    subelements[i] = el;
    subelement_verts[i] = vert;
  }

//  std::cout << "Finished precalculating lookup tables" << std::endl;

//  apf::FieldShape* mnew_shape = apf::getLagrange(1);

  // store all the subtriangles corresponding to an old mesh element
  apf::MeshEntity* subtriangles[numtriangles];
  apf::MeshEntity* node_entity;  // entity actually containing the node
  int nodenum; // node within element index
  int elnum;
  int new_elnum;
  int new_vertnum;
  int entity_node_idx; // index of node among all nodes on entities of same dimension
  int entity_idx;  // index of node containing entity in downward adjacencies
  int node_idx; // node index, pre offset
  int offset_node_idx; // node index, post offset
  int offset_i;
  int pos;
  apf::MeshEntity* e;
  apf::Downward down; // hold downward adjacent verts
  apf::Downward down_new; // hold downward adjacent verts of new mesh
  apf::MeshEntity* new_node_entity;
  int num_entities_per_dim[3] = { 3, 3, 1};
//apf::MeshEntity* getEntityFromNode(apf::Mesh* m, int typeOffsetsPerElement[], const int nodenum,  apf::MeshEntity* el, int entity_nodes_on[])


    int numcomp = apf::countComponents(field_old);
    const char* name_i = apf::getName(field_old);
    const char* name_new = apf::getName(field_new);
    int numcomp_new = apf::countComponents(field_new);

    if ( numcomp != numcomp_new)
    {
      std::cerr << "Warning: incompatable Field transfer requested" << std::endl;
      return;
    }
    
    std::cout << "\nCopying Field " << name_i << " to " << name_new << std::endl;


    double vals[numcomp];  // temporarily hold values

    apf::MeshIterator* itold = m->begin(2);
    apf::MeshIterator* itnew = m_new->begin(2);
    while ( (e = m->iterate(itold)) ) // loop over elements in old mesh
    {

      elnum = getNumber(numberings[2], e, 0, 0);

//      std::cout << "  processing old mesh element number " << elnum << std::endl;
      // get all the subtriangles of the current element
      for (int i = 0; i < numtriangles; ++i)
      {
        subtriangles[i] = m_new->iterate(itnew);
//        std::cout << "  subtriangle " << i << " pointer = " << subtriangles[i] << std::endl;
      }

      nodenum = 0;

      for (int dim = 0; dim < 3; ++dim)  // loop all dimensions
      {

//        std::cout << "  looping over nodes of dimension " << dim << std::endl;
        m->getDownward(e, dim, down);
/*        std::cout << "down = ";
        for (int k = 0; k < 12; ++k)
        {
          std::cout << "entry " << k << " = " << down[k] << std::endl;
        }
*/
        // loop over nodes of this dimension
        for (int j = 0; j < num_entities_per_dim[dim]*entity_nodes_on[dim]; ++j)
        {
//          std::cout << "    processing node " << j << std::endl;
          // get entity on old mesh containing the node
          entity_node_idx = nodenum - (typeOffsetsPerElement[dim] - 1);
          entity_idx = entity_node_idx/entity_nodes_on[dim]; // integer division
          node_idx = entity_node_idx % entity_nodes_on[dim];
          node_entity = down[entity_idx]; // get the mesh entity pointer

          // offset
//          std::cout << "    calculating offset" << std::endl;
          pos = nodenum + elnum*nnodes_per_el;
          offset_i = elementNodeOffsets[pos]; 
          offset_node_idx = abs(offset_i - (node_idx + 1)) - 1;
//          std::cout << "    offset_node_idx = " << offset_node_idx << std::endl;
          // get the entity on the new mesh containing the node
//          std::cout << "    retrieving new mesh entity" << std::endl;
          new_elnum = subelements[nodenum];
          new_vertnum = subelement_verts[nodenum];
//          std::cout << "    new_elnum = " << new_elnum << std::endl;
//          std::cout << "    new_vertnum = " << new_vertnum << std::endl;
          m_new->getDownward(subtriangles[new_elnum], 0, down_new);
          new_node_entity = down_new[new_vertnum];
//          std::cout << "    new_node_entity = " << new_node_entity << std::endl;
//          std::cout << "    node_entity = " << node_entity << std::endl;

//          std::cout << "    copying values" << std::endl;
          tri::copyFieldNode(field_old, field_new, node_entity, new_node_entity, offset_node_idx, 0, numcomp, vals);

          ++nodenum;
        }
      }

    } // end loop over old mesh elements
    
} // end function




apf::Mesh2* createSubMesh(apf::Mesh* m, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3])
{
// m is the existing (high order) mesh
// numtriangles is the number of triangles to break each large triangle into
// triangulation is a numtriangles x 3 array holding the indices of the nodes, 
// 1-based 
//   to use as the vertices of the sub-triangles
//   numberings is an array holding a numbering for each entity dimension

//  return NULL;
  std::cout << "\ncreating sub-triangulated mesh" << std::endl;

  apf::FieldShape* mshape = m->getShape();
  // get number of different types of entities
  std::size_t entity_counts[3] = {m->count(0), m->count(1), m->count(2)};

  // get number of nodes on different types of entities
  int entity_nodes_on[3] = {mshape->countNodesOn(0), mshape->countNodesOn(1), mshape->countNodesOn(2)};

  // create new mesh
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m_new = apf::makeEmptyMdsMesh(g, 2, false);
  pushMeshRef(m_new);

  // step 1: create all vertices of sub-triangles
   
  // allocate arrays to hold vertices of each order entity
  // these array should be in the Julia ordering, 
  // ie. getNumber(e_containing_edge) = first index
  apf::MeshEntity* verts[entity_counts[0]][entity_nodes_on[0]];
  apf::MeshEntity* edges[entity_counts[1]][entity_nodes_on[1]];
  apf::MeshEntity* faces[entity_counts[2]][entity_nodes_on[2]];

  std::cout << "creating vertices on new mesh" << std::endl;
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

//       std::cout << "creating vertex for vertex " << idx << " at coordinates " << coords << " pointer = " << verts[idx][0] << std::endl;
      }

      dim = 1;  
      it = m->begin(dim);  // iterate over edges
      while ( (e = m->iterate(it)) )
      {
        int idx = apf::getNumber(numberings[dim], e, 0, 0);
        for (int i = 0; i < entity_nodes_on[dim]; i++)
        {
          m->getPoint(e, i, coords);
          edges[idx][i] = m_new->createVert(0); // should be 2d array
          m_new->setPoint(edges[idx][i], 0, coords);

//          std::cout << "creating edge vertex for edge " << idx << ", node " << i << " at coordinates " << coords << " pointer = " << edges[idx][i]<< std::endl;
        }
      }

      dim = 2;  
      it = m->begin(dim);  // iterate over faces
      while ( (e = m->iterate(it)) )
      {
        int idx = apf::getNumber(numberings[dim], e, 0, 0);
        for (int i = 0; i < entity_nodes_on[dim]; i++)
        {
          m->getPoint(e, i, coords);

          faces[idx][i] = m_new->createVert(0); // should be 2d array
          m_new->setPoint(faces[idx][i], 0, coords);

//          std::cout << "creating edge vertex for face " << idx << ", node " << i << " at coordinates " << coords << " pointer = " << faces[idx][i] <<  std::endl;
        }
      }



    
  // step 2: create elements from vertices
  // loop over elements on m, create element on m_new
  // here we use the edge orientation information
  std::cout << "creating elements on new mesh" << std::endl;
  dim = 2;
  int elnum; // holds the element number
  int node;  // original node number
  int pos;  // linear address for elementNodeOffset
  const int nnodes_per_el = typeOffsetsPerElement[dim+1] - 1;
  int offset_j;
  apf::MeshEntity* el_verts[3]; // holds vertices of element to be created
  it = m->begin(dim);  // iterate over elements
  while ( (e = m->iterate(it)) )
  {
    elnum = apf::getNumber(numberings[dim], e, 0, 0); // zero base index
//    std::cout << "\n\nsubtriangulating element " << elnum << std::endl;
    for (int i=0; i < numtriangles; ++i)  // loop over all subtriangles
    {
//      std::cout << "\ncreating sub triangle " << i << std::endl;
      // get 3 vertices
      for (int j=0; j < 3; ++j)
      {
        node = triangulation[i][j] - 1; // one based index
//        std::cout << "\nnode = " << node << std::endl;
        // calculate linear offset for elementNode offsets
        pos = node + elnum*nnodes_per_el;
//        std::cout << "offset pos = " << pos << std::endl;
        offset_j = elementNodeOffsets[pos];
        
//        std::cout << "offset = " << offset_j << std::endl;

        el_verts[j] = tri::getVert(m, &verts[0][0], &edges[0][0], &faces[0][0], typeOffsetsPerElement, node, offset_j, e, numberings, entity_nodes_on );
//        std::cout << "  el_vert " << j << " = " << el_verts[j] << std::endl;
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

  tri::transferNumberings(m, m_new, numtriangles, triangulation, elementNodeOffsets, typeOffsetsPerElement, numberings);

  // write visualization file
  apf::writeVtkFiles("newmesh_linear", m_new);

  return m_new;
}


