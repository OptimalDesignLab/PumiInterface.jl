// this file contains the code necessary to perform the subtriangulatio of a 
// high order mesh into a low order one

#include "triangulationDG.h"

namespace triDG {
#define N_VERT_PER_EL 3  // triangles
//-----------------------------------------------------------------------------
// Function declarations (private functions)
//-----------------------------------------------------------------------------

/*
 * Gets a pointer to the vertex MeshEntity* on the new mesh
 * Inputs:
 *   m : old mesh
 *   verts[]: array of vertex MeshEntity* on the new mesh for this element
 *   edges[]: array of edge MeshEntity* on the new mesh for this element
 *   faces[]: array of face MeshEntity* on the new mesh for this element
 *   typeOffsetsPerElement: the starting index of each type of MeshEntity in 
 *                          an array of of all the entities on the element
 *   nodenum: the index on the current element of the node you want to get the 
 *            vertex MeshEntity* on the new mesh of
 *   offset: the offset of the node (used for remapping shared nodes on a CG mesh)
 *   el: the MeshEntity* of the element on the old mesh
 *   numberings[]: array of the Numbering* [vertnums, edgenums, elnums] on the old
 *                 mesh
 *   entity_nodes_on[]: array of number of nodes on each type of MeshEntity*
 *
 * Outputs:
 *   the MeshEntity* of the vertex on the new mesh
*/
apf::MeshEntity* getVert(apf::Mesh* m, apf::MeshEntity* verts[], apf::MeshEntity* edges[], apf::MeshEntity* faces[], int typeOffsetsPerElement[], const int nodenum, uint8_t offset,  apf::MeshEntity* el, apf::Numbering* numberings[], int entity_nodes_on[]);

/*
 * Copy the values of a Numbering at a node from old mesh to the new mesh, checking
 * if the old mesh has the specified node and has  assigned a value to it
 * Inputs:
 *   n_old: the Numbering* on the old mesh
 *   n_new: the Numbering* on the new mesh
 *   e_old: the MeshEntity* containing the node on the old mesh
 *   e_new: the MeshEntity* containing the node on the new mesh
 *   node_old: the index of the node on the old mesh entity
 *   node_new: the index of hte node on the new mesh entity
 *   ncomp: number of components in the numbering
 * Outputs:
 *   1 if values were copied, zero otherwise
*/   
int copyNumberingNode(apf::Numbering* n_old, apf::Numbering* n_new, apf::MeshEntity* e_old, apf::MeshEntity* e_new, const int node_old, const int node_new, int ncomp);


/*
 * Copy the values of a field at a node from the old mesh to the new mesh
 * Inputs:
 *   f_old: the Field* on the old mesh
 *   f_new: the Field* on the new mesh
 *   e_old: the MeshEntity* containing the node on the old mesh
 *   e_new: the MeshEntity* containing the node on the new mesh
 *   node_old: the index of the node on the old mesh entity
 *   ncomp: number of components in the numbering
*/
void copyFieldNode(apf::Field* field_old, apf::Field* field_new, apf::MeshEntity* e_old, apf::MeshEntity* e_new, const int node_old, const int node_new, int ncomp, double vals[]);


/* Prints all the fields of a mesh to STDOUT
 * Inputs:
 *   m: a Mesh*
*/
void printFields(apf::Mesh* m);

/* Copies a Numbering with FieldShape apf::getConstant(2) (ie. values only on
 *        elements) to all of the corresponding elements of the new mesh
 * Inputs:
 *   m: the old Mesh*
 *   m_new: the new Mesh*
 *   numtriangles: the number of triangles each old mesh element is divided into
 *   numbering_old: the Numbering* on the old mesh
*/
void transferNumberingToElements(apf::Mesh* m, apf::Mesh* m_new, const int numtriangles, apf::Numbering* numbering_old);


/* Copies all Numberings on the old mesh to the new mesh.  Certain element Numberings
 *        are copied using transferNumberingToElements
 * Inputs:
 *   m: the old mesh
 *   m_new: the new mesh
 *   numtriangles: the number of triangles each old mesh element is divided into
 *   triangulation: the numtriangles x 3 array that specifies which vertices and 
 *                  nodes of an old mesh element are used as the vertices of the
 *                  new mesh element (one based indexing)
 *   elementNodesOffsets[]: the standard thing
 *   numberings[]: the array of Numbering* for each dimension entity
*/
void transferNumberings(apf::Mesh* m, apf::FieldShape* mshape, apf::Mesh* m_new, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3]);

/* Calculates the linear index (0-based) for an element of an array
 * Inputs:
 *   i: the row number (0-based)
 *   j: the column number (1-based)
 *   si: the number of rows
 *   sj: the number of columns
 * Outputs:
 *   the index
*/
int getindex(int i, int j, int si, int sj);

/* Calculates the linear index (0-based) for an element of a 3D array stored
 *   column major
 * Inputs:
 *   i, j, k: indices, zero based
 *   si, sj, sk: array dimensions
 * Output:
 *   the index
*/
int getindex3c(const int i, const int j, const int k, const int si, const int sj, const int sk);

/*
 * Gets all the values of a Field on an element and puts them in an array
 * Inputs:
 *   field: the Field*
 *   el:  the MeshEntity* element
 * Inputs/Outputs:
 *   vals: a num_nodes_per_el x ncomp array to store the values in, where ncomp
 *         is the number of components of the field
*/
void getFieldElementVals(apf::Field* field, apf::MeshEntity* el, double* vals);


/* small matrix-matrix multiplication Ax = b
 * Inputs
 *   A: m x n matrix
 *   x: n x p matrix
 *   m: array dimension 
 *   n: array dimension
 *   p: array dimension
 * Inputs/Outputs:
 *   b: m x p matrix to store result in
*/
void smallMatMat(const double* A, const double* x, double* b, const int m, const int n, const int p);

/* Interpolates field values to the vertices of the new mesh
 * Inputs:
 *   m: the old mesh
 *   m_new: the new mesh
 *   numtriangles: see other explanations
 *   trianglulation: see other explanations
 *   elementNodeOffsets: the offsets used to remap nodes on shared entities
 *   typeOffsetsPerElement[]: see other explanations
 *   numberings[]: see other explanations
 *   field_old: the field on the old mesh
 *   field_new: the field on the new mesh
*/
void transferVertices(apf::Mesh* m, apf::Mesh* m_new, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3], apf::Field* field_old, double* interp_op, apf::Field* field_new);


// calculate the lookup tables that map:
// vertex number on old mesh element -> subtriangle number, vert index
// on new mesh
// for vertices only
//
// Inputs:
//    triangulation: a numtriangles x 3 array containing the indices of the 
//                   nodes on the old mesh that will be used to create the
//                   vertices of the new mesh
//    numtriangles:  the number of triangles each old mesh element is divided
//                  into 
// Outputs:
//    subelements: array to be populated with the subelement index (0-based)
//    subelement_verts: array to be populated with the index of the vertex
//                     on the sub element
void getVertLookupTables(const int triangulation[][3], const int numtriangles, int* subelements, int* subelement_verts);


// precalculate tables
// node on old mesh  -> local element number, vertex index on new mesh
// where local elemlent number tells which element it is of all the 
// subelements created for the current old mesh element
// This lookup table skips any nodes that might be on vertices, because
// they are handled through getVertLookupTables and transferVertices
//
/*
  Inputs:
    nnodes_per_element: number of nodes on an element
    triangulation: a numtriangles x 3 array containing the indices of the 
                   nodes on the old mesh that will be used to create the
                   vertices of the new mesh
    numtriangles:  the number of triangles each old mesh element is divided
                   into 
  Outputs:
    subelements: array to be populated with the indices of the elements on the
                 new mesh that are part of the element on the old mesh
    subelement_verts: array to be popultes with the indices of the vertex of
                      the new mesh vertex correspondong to the old mesh node
*/
void getFieldLookupTables(const int nnodes_per_el, const int triangulation[][3], const int numtriangles, int* subelements, int* subelement_verts);


/* Gets the maximum number present in a Numbering
 * Inputs:
 *   m_new: a mesh
 *   n_new: a Numbering*
 */
int getMaxNumber(apf::Mesh* m_new, apf::Numbering* n_new);



/* Get the mapping old edge local number -> subelement, subelement edge local
 * number
 * Inputs:
 *   triangulation: the triangulation
 *   numtriangles: the number of triangles
 *   subelements: array of length 3 to be populated with the subelement numbers
 *   subelement_edges: array of length 3 to be populated with the subelement
 *                     edge local numbers
 *   num_edge_nodes: number of nodes on an edge (on the FieldShape of the 
 *                   triangulation)
 */
void getEdgeLookupTables(const int triangulation[][3], const int numtriangles, int* subelements, int* subelement_edges, const int num_edge_nodes=0);

/* Determines if an array contains 2 specified values
 * Inputs:
 *   arr: the array
 *   len: length of the array
 *   val1: the first value
 *   val2: the second value
 */
bool contains2( int* arr, int len, const int val1, const int val2);

/* Determines the local edge number of an edge defined by 2 vertices
 * Inputs:
 *   arr: array of vertiex indices
 *   val1: first vertex
 *   val2: second vertex
 * Outputs:
 *   -1 if edge not found, otherwise the local edge number
 */
int getEdgePos(const int arr[3], const int val1, const int val2);


/* Determine the orientation of an edge relative to the element
 * Inputs:
 *   m: the mesh
 *   el: the element
 *   edge_idx: the local index of the edge
 * Outputs:
 *   bool: true if edge is reversed, false otherwise
 */
bool getEdgeOrientation(apf::Mesh* m, apf::MeshEntity* el, int edge_idx);

//-----------------------------------------------------------------------------
// Function definitions
//-----------------------------------------------------------------------------


// given nodenum, the index of a node within an element (zero based), return the MeshEntity (make sure the N_VERT_PER_EL is removed!)
// pointer for the corresponding vert on the new mesh
// typeOffsetsPerElement are the (1-based) indicies of the start index of the nodes
// of each entity type on an element
// el is the pointer to the element to which the node belongs
// entity_nodes_on is the number of nodes on verts, edges, faces
apf::MeshEntity* getVert(apf::Mesh* m, apf::MeshEntity*** verts, apf::MeshEntity*** edges, apf::MeshEntity*** faces, int typeOffsetsPerElement[], const int nodenum, uint8_t offset,  apf::MeshEntity* el, apf::Numbering* numberings[], int entity_nodes_on[])
{
  int entity_idx;  // local index of entity containing the node
  apf::MeshEntity* e;
  int entity_num;  // global number of entity
  apf::Downward down;
//  int pos;  // linear address within verts, edges, faces
  int entity_node_idx;  // index of node from 0 to numnodes on all entity of this type
  int entity_node_idx_local; // index from 0 to nnodes on this entity
  int entity_node_offset_idx_local; // offset local idx
  int dim;  // dimension of entity containing the node

//  std::cout << "getting new mesh vertex for nodenum " << nodenum << std::endl;
  // this assumes there are nodes on vertices
  if (nodenum < ( typeOffsetsPerElement[1] - 1))  // first 3 "nodes" are verts
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

//    pos = entity_num*entity_nodes_on[dim] + entity_node_offset_idx_local;
    return verts[entity_num][entity_node_offset_idx_local];
//    return verts[pos];

  } else if (nodenum < ( typeOffsetsPerElement[2] - 1))
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

//    pos = entity_num*entity_nodes_on[dim] + entity_node_offset_idx_local;

    return edges[entity_num][entity_node_offset_idx_local];  
//    return edges[pos];

  } else if (nodenum < ( typeOffsetsPerElement[3] - 1))
  {

//    std::cout << "face node" << std::endl;
    dim = 2;
    entity_node_idx = nodenum - (typeOffsetsPerElement[dim] - 1);
//    std::cout << "entity_node_idx = " << entity_node_idx << std::endl;
//    std::cout << "entity_nodes_on[dim] = " << entity_nodes_on[dim] << std::endl;
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
//    std::cout << "e = " << e << std::endl;
    entity_num = apf::getNumber(numberings[dim], e, 0, 0);
//    std::cout << "entity_num = " << entity_num << std::endl;

//    pos = entity_num*entity_nodes_on[dim] + entity_node_offset_idx_local;
//    std::cout << "pos = " << pos << std::endl;
    return faces[entity_num][entity_node_offset_idx_local];
//    return faces[pos];

  } else {
    std::cerr << "Warning: in getVert,  nodenum too high, returning NULL" << std::endl;
    return NULL;
  }
}
 

// copy all of the Numbering values at a particular node from an old Numbering to a new one.
int copyNumberingNode(apf::Numbering* n_old, apf::Numbering* n_new, apf::MeshEntity* e_old, apf::MeshEntity* e_new, const int node_old, const int node_new, int ncomp)
{
  apf::Mesh* m = getMesh(n_old);
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
//      std::cout << "copying value from entity " << e_old << " node " << node_old<< " to entity " << e_new << " node " << node_new << std::endl;
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
//  std::cout << "copying from " << e_old << " node " << node_old;
//  std::cout << " to " << e_new << " node " << node_new << std::endl;

      apf::getComponents(field_old, e_old, node_old, vals);
//      std::cout << "    retrieved components, vals = ";
//      for (int i = 0; i < ncomp; ++i)
//      {
//        std::cout << vals[i] << ", ";
//      }
//      std::cout << std::endl;

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
//    std::cout << "processing element " << e << std::endl;
      for (int i = 0; i < numtriangles; ++i)
      {
        e_new = m_new->iterate(itnew);
//          std::cout << "    e_new = " << e_new << std::endl;

//        std::cout << "  triangle " << i << std::endl;
        for (int n = 0; n < numcomp; ++n)
        {
//          std::cout << "    component " << n << std::endl;
          val = apf::getNumber(numbering_old, e, 0, n);
          apf::number(numbering_new, e_new, 0, n, val);
        }
      }
  }

}  // end function


void transferNumberings(apf::Mesh* m, apf::FieldShape* mshape, apf::Mesh* m_new, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3])
{
  // compute some quantities
  int nnodes_per_el = typeOffsetsPerElement[3] - 1;
//  apf::FieldShape* mshape = m->getShape();
  // get number of different types of entities
//  std::size_t entity_counts[3] = {m->count(0), m->count(1), m->count(2)};

  // get number of nodes on different types of entities
  int entity_nodes_on[3] = {mshape->countNodesOn(0), mshape->countNodesOn(1), mshape->countNodesOn(2)};

  const int dimensionTypes[3] = {apf::Mesh::VERTEX, apf::Mesh::EDGE, apf::Mesh::TRIANGLE};

  // get the node index of the first node on each entity of each dimension
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
  // name string and having a table of function pointers here`
  // for now we will just stick with linear Lagrange
  // the number of components ofa  field can be found with
  // apf::countComponents(Field* f)


  // precalculate tables needed to lookup a given node on the new mesh
  int subelements[nnodes_per_el]; // hold index (0-based) of subtriangles
  int subelement_verts[nnodes_per_el]; // holds vertex index of subtriangle

  getFieldLookupTables(nnodes_per_el, triangulation, numtriangles, subelements, subelement_verts);


   // copy numberings
  int numnumberings = m->countNumberings();

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
  apf::Downward down; // hold downward adjacent verts of old mesh
  apf::Downward down2; // hold downwarda adjacent verts of new mesh
  apf::MeshEntity* new_node_entity;
//  int num_entities_per_dim[3] = { 3, 3, 1};
//apf::MeshEntity* getEntityFromNode(apf::Mesh* m, int typeOffsetsPerElement[], const int nodenum,  apf::MeshEntity* el, int entity_nodes_on[])

  for (int i = 0; i < numnumberings; ++i)
  {
    std::cout << "Copying Numbering " << i << std::endl;

    apf::Numbering* numbering_i = m->getNumbering(i);
    int numcomp = apf::countComponents(numbering_i);
    const char* name_i = apf::getName(numbering_i);
    std::cout << "Numbering name is " << name_i << std::endl;
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
      std::cout << "transferring numbering to elements" << std::endl;
      transferNumberingToElements(m, m_new, numtriangles, numbering_i);
      continue;
    }

    // create new field and populate it
    apf::FieldShape* fshape_new = getNewShape(numbering_i_shape);
    std::cout << "creating new Numbering with name: " << fshape_new->getName() << std::endl;
    apf::Numbering* n_new = apf::createNumbering(m_new, name_i, fshape_new, numcomp);

    // if Numbering has values on vertices, transfer them here
    if (numbering_i_shape->hasNodesIn(0)) 
    {
      std::cout << "transfering vertices" << std::endl;
      transferVertices(m, m_new, numtriangles, triangulation, elementNodeOffsets, typeOffsetsPerElement, numbering_i, n_new);
    }

    if (numbering_i_shape->hasNodesIn(1))
    {
      std::cout << "transferring edges" << std::endl;
      transferEdges(m, mshape, m_new, numtriangles, triangulation, elementNodeOffsets, typeOffsetsPerElement, numbering_i, n_new);
    }

    // now transfer all non-vertex values
    std::cout << "transferring elements" << std::endl;
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

      for (int dim = 2; dim < 3; ++dim)  // loop over elements
      {

//        std::cout << "  looping over dimension " << dim << std::endl;
        m->getDownward(e, dim, down);

        // skip this dimension if no nodes on the old mesh Numbering
        if (!numbering_i_shape->hasNodesIn(dim))
        {
          continue;
        }



//        for (int tmp = 0; tmp < 12; ++tmp)
//        {
//          std::cout << "down[" << tmp << "] = " << down[tmp] << std::endl;
//        }

          // loop over entities on this dimension
//          int nume = m -> adjacentCount[apf::Mesh::TRIANGLE][dim];
//          std::cout << "  number of entities on this dimension = " << nume << std::endl;
          for (int entity_d = 0; entity_d < m -> adjacentCount[apf::Mesh::TRIANGLE][dim]; ++entity_d)
          {
//            std::cout << "  on entity " << entity_d << std::endl;
            // get the starting node number for this entity on the old mesh
            // in the fieldshape of the mesh
            // This has to be reset here because the fieldshape of the numbering
            // and the mesh might be different
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
//              std::cout << "entity_nodes_on[dim] = " << entity_nodes_on[dim] << std::endl;
              entity_idx = entity_node_idx/entity_nodes_on[dim]; // integer division
//              std::cout << "entity_idx = " << entity_idx << std::endl;
              node_idx = entity_node_idx % entity_nodes_on[dim];
              node_entity = down[entity_idx]; // get the mesh entity pointer
//              std::cout << "node_entity = " << node_entity << std::endl;

              // offset
//              std::cout << "    calculating offset" << std::endl;
              pos = nodenum + elnum*nnodes_per_el;
              offset_i = elementNodeOffsets[pos]; 
              offset_node_idx = abs(offset_i - (node_idx + 1)) - 1;

              // get the entity on the new mesh containing the node
//              std::cout << "    retrieving new mesh entity" << std::endl;
              new_elnum = subelements[nodenum];
              new_vertnum = subelement_verts[nodenum];
//              std::cout << "    new_elnum = " << new_elnum << std::endl;
//              std::cout << "    new_vertnum = " << new_vertnum << std::endl;
              m_new->getDownward(subtriangles[new_elnum], 0, down2);
              new_node_entity = down2[new_vertnum];
//              std::cout << "    new_node_entity = " << new_node_entity << std::endl;

//              std::cout << "    copying values" << std::endl;
              numbered_count += copyNumberingNode(numbering_i, n_new, node_entity, new_node_entity, offset_node_idx, 0, numcomp);

              ++nodenum;
              }  // end loop over nodes on this entity
           } // end loop over entities on this dimension
         } // end loop over dimensions

    } // end loop over old mesh elements
 
//    std::cout << "Complete numbering" << std::endl; 
    completeNumbering(m_new, n_new);

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

  std::cout << "printing old mesh fields and numberings" << std::endl;
  triDG::printFields(m);

  std::cout << "printing new mesh fields and numberings" << std::endl;
  triDG::printFields(m_new);

}  // end function


// computes the linear index of a 2D array from 2 indices and the dimensions
// of the array
// assumes row major storage
int getindex(int i, int j, int si, int sj)
{
  return j + i*sj;
}

// compute the linear index of an element in a 3D array, stored column major
// all indices are zero based
int getindex3c(const int i, const int j, const int k, const int si, const int sj, const int sk)
{
  return k*si*sj + j*si + i;
}


// get all the values of a field on a particular element
// el is the the element 
// vals is a 2D array that must be large enough to store all the values
// vals must be nnodes per element x ncomp
void getFieldElementVals(apf::Field* field, apf::MeshEntity* el, double* vals)
{
  apf::FieldShape* fshape = apf::getShape(field);
  int ncomp = apf::countComponents(field);
  const int entity_nodes_on[3] = {fshape->countNodesOn(0), fshape->countNodesOn(1), fshape->countNodesOn(2)};

//  int nnodes = entity_nodes_on[0]*3 + entity_nodes_on[1]*3 + entity_nodes_on[2];

  int curr_node = 0; // current node (of all nodes on the element)
  for (int dim = 0; dim < 3; ++dim)  // loop over dimensions
  {
    // loop over nodes of this dimension
    for (int node = 0; node < entity_nodes_on[dim]; ++node)
    {
      apf::getComponents(field, el, node, &vals[curr_node]);
      curr_node += ncomp;
    }
  }

} // end function

// small matrix-matrix multiplication Ax = b
// A must be m x n
// x must be n x p
// b must be m x p
void smallMatMat(const double* A, const double* x, double* b, const int m, const int n, const int p)
{
  // dimensions of the matrices
  // these should get constant folded
  const int Arows = m;
  const int Acols = n;
  const int xrows = n;
  const int xcols = p;
  const int brows = m;
  const int bcols = p;

  double tmp = 0;  // stack allocate the tmp
  int A_idx = 0;  // linear index in A
  int x_idx = 0;  // linear index in x
  int b_idx = 0;  // linear index in b

  for (int k = 0; k < xcols; ++k)  // loop over columns of x
  {
    for (int i = 0; i < m; ++i)  // loop over rows of A
    {
      // dot product of this row of A with this column of x
      for ( int j = 0; j < n; ++j)
      {
        A_idx = getindex(i, j, Arows, Acols);
        x_idx = getindex(j, k, xrows, xcols);
        tmp += A[A_idx]*x[x_idx];
      }
      b_idx = getindex(i, k, brows, bcols);
      b[b_idx] = tmp;
      tmp = 0.0;
    }
  }

} // end function


// transfer (possibly interpolated) field values from the old mesh to the
// vertices of the new mesh
// not interpolating is currently not supported
void transferVertices(apf::Mesh* m, apf::Mesh* m_new, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3], apf::Field* field_old, double* interp_op, apf::Field* field_new)
{

  std::cout << "entered transferVertices" << std::endl;
  // the lookup tables
  int subelements[N_VERT_PER_EL];
  int subelement_verts[N_VERT_PER_EL];
  getVertLookupTables(triangulation, numtriangles, subelements, subelement_verts);
  apf::MeshEntity* subtriangles[numtriangles];

  const int ncomp = apf::countComponents(field_old);
  const int nnodes_per_el = typeOffsetsPerElement[3] - 1;
  double node_vals[nnodes_per_el][ncomp];  // hold values at node locations
  double vert_vals[N_VERT_PER_EL][ncomp];  // hold values interpolated to verts
  apf::MeshIterator* itold = m->begin(2);
  apf::MeshIterator* itnew = m_new->begin(2);
  apf::MeshEntity* el;  // hold old mesh element
  apf::MeshEntity* el_new;  // hold new mesh element
  apf::MeshEntity* vert_new;  // hold vertex of new element;
  apf::Downward down;  // hold vertices of new mesh element

  int e_cnt = 0;
//  std::cout << "about to start looping over elements" << std::endl;
  while ( (el = m->iterate(itold)) )  // loop over elements of old mesh
  {
//    std::cout << "processing element " << e_cnt << std::endl;
    // do interpolation, store result in vert_vals
    getFieldElementVals(field_old, el, &node_vals[0][0]);
    smallMatMat(interp_op, &node_vals[0][0], &vert_vals[0][0], N_VERT_PER_EL, nnodes_per_el, ncomp);

    // get subelements on new mesh
    for (int i = 0; i < numtriangles; ++i)
    {
      subtriangles[i] = m_new->iterate(itnew);
    }

    for (int i = 0; i < N_VERT_PER_EL; ++i)  // loop over vertices of old mesh
    {
//      std::cout << "  vert " << i << std::endl;

      // get the element, vertex on new mesh
      el_new = subtriangles[subelements[i]];
      m_new->getDownward(el_new, 0, down);
      vert_new = down[ subelement_verts[i] ];

      // copy values to new mesh
//      std::cout << "copying to " << vert_new << std::endl;
      apf::setComponents(field_new, vert_new, 0, vert_vals[i]);
    }  // ned loop over vertices of old mesh
    ++e_cnt;
  }  // end loop over old mesh elements
        
}

// transfer Numbering values from the old mesh to the
// vertices of the new mesh
// This function does not use elementNodeOffsets, because there are not
// nodes on vertices in the old mesh
void transferVertices(apf::Mesh* m, apf::Mesh* m_new, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* n_old, apf::Numbering* n_new)
{
  // the lookup tables
  int subelements[N_VERT_PER_EL];
  int subelement_verts[N_VERT_PER_EL];
  getVertLookupTables(triangulation, numtriangles, subelements, subelement_verts);
  apf::MeshEntity* subtriangles[numtriangles];
  // now do the transfer
  const int ncomp = apf::countComponents(n_old);
  apf::MeshIterator* itold = m->begin(2);
  apf::MeshIterator* itnew = m_new->begin(2);
  apf::MeshEntity* el;  // hold old mesh element
  apf::MeshEntity* el_new;  // hold new mesh element
  apf::MeshEntity* vert_new;  // hold vertex of new element;
  apf::MeshEntity* vert_old; // hold vertex on old mesh element
  apf::Downward down2;  // hold vertices of old mesh element
  apf::Downward down;  // hold vertices of new mesh element
  while ( (el = m->iterate(itold)) )  // loop over elements of old mesh
  {
    // get the old mesh vertex
    m->getDownward(el, 0, down2);

    // get subelements on new mesh
    for (int i = 0; i < numtriangles; ++i)
    {
      subtriangles[i] = m_new->iterate(itnew);
    }

    for (int i = 0; i < N_VERT_PER_EL; ++i)  // loop over vertices of old mesh
    {
      // get the vertex on the old mesh
      vert_old = down2[i];

      // get the element, vertex on new mesh
      el_new = subtriangles[subelements[i]];
      m_new->getDownward(el_new, 0, down);
      vert_new = down[ subelement_verts[i] ];

      // copy values to new mesh
      copyNumberingNode(n_old, n_new, vert_old, vert_new, 0, 0, ncomp);
    }  // ned loop over vertices of old mesh
  }  // end loop over old mesh elements
        
}

// transfer Numbering values from the old mesh to the
// edges of the new mesh
// This function does not use elementNodeOffsets, because there are not
// nodes on edges in the old mesh
// the new Numbering must have the same number of nodes per edge as the 
// old Numbering
void transferEdges(apf::Mesh* m, apf::FieldShape* mshape, apf::Mesh* m_new, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* n_old, apf::Numbering* n_new)
{

  int subelements[3];
  int subelement_edges[3];
  int num_edge_nodes = mshape->countNodesOn(1);
  getEdgeLookupTables(triangulation, numtriangles, subelements, subelement_edges, num_edge_nodes);

  apf::FieldShape* nshape = apf::getShape(n_old);
  const int nnodes_per_edge = nshape->countNodesOn(m->simplexTypes[1]);
  const int ncomp = apf::countComponents(n_old);
  apf::MeshEntity* subtriangles[numtriangles];
  apf::MeshIterator* itold = m->begin(2);
  apf::MeshIterator* itnew = m_new->begin(2);
  apf::MeshEntity* el;  // hold old mesh element
  apf::MeshEntity* el_new;  // hold new mesh element
  apf::MeshEntity* edge_new;  // hold edge of new element;
  apf::MeshEntity* edge_old; // hold edge on old mesh element
  apf::Downward down2;  // hold edges of old mesh element
  apf::Downward down;  // hold edges of new mesh element
  bool old_edge_reversed; // is the old edge reversed
  bool new_edge_reversed; // is the new edge reversed
  int orient_diff; 
  int offset;  // offset (used to reverse the edge node ordering


  int el_cnt = 0;
  while ( (el = m->iterate(itold)) )
  {
    // get old mesh edges
    m->getDownward(el, 1, down2);

    for (int i = 0; i < numtriangles; ++i)
    {
      subtriangles[i] = m_new->iterate(itnew);
    }

    for (int i = 0; i < 3; ++i)  // loop over edges of old mesh element
    {
      edge_old = down2[i]; // get old mesh edge

      // get new mesh edge
      el_new = subtriangles[subelements[i]];

      m_new->getDownward(el_new, 1, down);
      edge_new = down[subelement_edges[i]];

//      std::cout << "got new mesh edge" << std::endl;

      // get the edge orientations
      old_edge_reversed = getEdgeOrientation(m, el, i);
      new_edge_reversed = getEdgeOrientation(m_new, el_new, subelement_edges[i]);

//      std::cout << "got edge orientations" << std::endl;

      orient_diff = old_edge_reversed - new_edge_reversed;
      offset = orient_diff*(nnodes_per_edge - 1);

//      std::cout << "orient_diff = " << orient_diff << std::endl;
//      std::cout << "offset = " << offset << std::endl;

      for (int node = 0; node < nnodes_per_edge; ++node)
      {
        copyNumberingNode(n_old, n_new, edge_old, edge_new, node, offset - node, ncomp);
      }
    }  // end loop over element edges
    ++el_cnt;
  }  // end loop over old mesh elements

}  // end function
// calculate the lookup tables that map:
// vertex number on old mesh element -> subtriangle number, vert index
// on new mesh
// for vertices only
//
// Inputs:
//    triangulation: a numtriangles x 3 array containing the indices of the 
//                   nodes on the old mesh that will be used to create the
//                   vertices of the new mesh
//    numtriangles:  the number of triangles each old mesh element is divided
//                  into 
// Outputs:
//    subelements: array to be populated with the subelement index (0-based)
//    subelement_verts: array to be populated with the index of the vertex
//                     on the sub element
void getVertLookupTables(const int triangulation[][3], const int numtriangles, int* subelements, int* subelement_verts)
{
  int el = 0;  // store found element
  int vert = 0; // store found node
  bool breaknow = false;  // flag for breaking second loop
//  std::cout << "Precalculating sub-triangle lookup tables" << std::endl;
  for (int i = 0; i < N_VERT_PER_EL; ++i)  // search for node i
  {
    // loop through the triangulation
    for ( int j = 0; j < numtriangles; ++j)  // loop over subtriangles
    {
      for ( int k = 0; k < N_VERT_PER_EL; ++k)  // loop over nodes on the triangle
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
}

// precalculate tables
// node on old mesh  -> local element number, vertex index on new mesh
// where local element number tells which element it is of all the 
// subelements created for the current old mesh element
// This lookup table skips any nodes that might be on vertices, because
// they are handled through getVertLookupTables and transferVertices
//
/*
  Inputs:
    nnodes_per_element: number of nodes on an element
    triangulation: a numtriangles x 3 array containing the indices of the 
                   nodes on the old mesh that will be used to create the
                   vertices of the new mesh
    numtriangles:  the number of triangles each old mesh element is divided
                   into 
  Outputs:
    subelements: array to be populated with the indices of the elements on the
                 new mesh that are part of the element on the old mesh
    subelement_verts: array to be popultes with the indices of the vertex of
                      the new mesh vertex correspondong to the old mesh node
*/
void getFieldLookupTables(const int nnodes_per_el, const int triangulation[][3], const int numtriangles, int* subelements, int* subelement_verts)
{
//  int subelements[nnodes_per_el]; // hold index (0-based) of subtriangles
//  int subelement_verts[nnodes_per_el]; // holds vertex index of subtriangle
  int el = 0;  // store found element
  int vert = 0; // store found node
  bool breaknow = false;  // flag for breaking second loop
//  std::cout << "Precalculating sub-triangle lookup tables" << std::endl;
  for (int i = 0; i < nnodes_per_el; ++i)  // search for node i
  {
    // loop through the triangulation
    for ( int j = 0; j < numtriangles; ++j)  // loop over subtriangles
    {
      for ( int k = 0; k < 3; ++k)  // loop over nodes on the triangle
      {
        // check if we have found the node we are looking for
        // -1 converts triangulation to 1 based indices
        // - N_VET_PER_EL accounts for the fact that the triangulation
        // numbers the "vertices" as the first 3 nodes, even though
        // vertices don't exist
        if ( (triangulation[j][k] - 1 - N_VERT_PER_EL) == i)
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
}


// map an local edge number on the old mesh to the local element number and
// local edge number
void getEdgeLookupTables(const int triangulation[][3], const int numtriangles, int* subelements, int* subelement_edges, const int num_edge_nodes)
{
  int vert1;
  int vert2;
  int ret_val;

  for (int edge_old = 0; edge_old < 3; ++edge_old)
  {
    vert1 = edge_old;
    if (num_edge_nodes > 0)
    {
      vert2 = edge_old + (num_edge_nodes-1)*edge_old + N_VERT_PER_EL;
    } else
    {
      vert2 = (edge_old + 1) % 3;
    }

    // check all triangles for this edge
    for (int subel = 0; subel < numtriangles; ++ subel)
    {
      // +1 because triangulation is 1-based
      ret_val = getEdgePos( triangulation[subel], vert1 +1, vert2 +1);
      if (ret_val >= 0)  // found edge
      {
        subelements[edge_old] = subel;
        subelement_edges[edge_old] = ret_val;
        continue;
      }
    }
  }


}  // end function


// does an array contain 2 specified values
bool contains2(const int* arr, int len, const int val1, const int val2)
{
  bool found1 = false;
  bool found2 = false;
  int tmp;
  for (int i=0; i < len; ++i)
  {
    tmp = arr[i];
    if (tmp == val1) {found1 = true;}
    if (tmp == val2) {found2 = true;}
  }

  return found1 && found2;
}

// get the local edge number of the edge defined by two values
int getEdgePos(const int arr[3], const int val1, const int val2)
{
  int pos1=-1;
  int pos2=-1;
  int tmp;

  for (int i=0; i < 3; ++i)
  {
    tmp = arr[i];
    if (tmp == val1) {pos1 = i;}
    if (tmp == val2) {pos2 = i;}
  }

  if (pos1 == -1 || pos2 == -1)
  {
    return -1;  // edge not found
  }

  if (pos1 == 0 && pos2 == 1)
  {
    return 0;
  }
  else if (pos1 == 1 && pos2 == 0)
  {
    return 0;
  }
  else if (pos1 == 1 && pos2 == 2)
  {
    return 1;
  }
  else if (pos1 == 2 && pos2 == 1)
  { 
    return 1;
  }
  else if (pos1 == 2 && pos2 == 0) 
  {
    return 2;
  }
  else if (pos1 == 0 && pos2 == 2)
  { 
    return 2;
  }
  else 
  { 
    return -1;
  }

} // end function

// get the edge orientation of a specified edge (true = revered relative to
// element, false = same as element)
bool getEdgeOrientation(apf::Mesh* m, apf::MeshEntity* el, int edge_idx)
{
  apf::Downward downverts;
  apf::Downward downedges;
  m->getDownward(el, 0, downverts);
  m->getDownward(el, 1, downedges);
  apf::MeshEntity* edge = downedges[edge_idx];

  m->getDownward(edge, 0, downedges); // reuse downedges

  apf::MeshEntity* vert1 = downedges[1];
  apf::MeshEntity* vert2 = downedges[2];

  int vert1pos = 0;
  int vert2pos = 1;

  // find positions of verts within downverts
  for (int i = 0; i < 3; ++i)
  {
    if (downverts[i] == vert1) { vert1pos = i;}
    if (downverts[i] == vert2) { vert2pos = i;}
  }

  // figure out what orientation the edge is in by comparing vert positions

  if (edge_idx == 1 || edge_idx == 2)
  {
    if (vert2pos > vert1pos)
    {
      return false;
    } else
    {
      return true;
    }
  } else // edge_idx  == 3
  {
    if (vert2pos > vert1pos)
    {
      return true;
    } else
    {
      return false;
    }
  }  

} // end function
void printArray(std::ostream& os, apf::MeshEntity** arr, const int si, const int sj)
{
  int idx = 0;
  for (int i = 0; i < si; ++i)  // loop rows
  {
    for (int j = 0; j < sj;  ++j)  // loop columns
    {
      idx = getindex(i, j, si, sj);
      os << arr[idx] << " ";
    }
    os << std::endl;
  }

}  // end function


void completeNumbering(apf::Mesh* m_new, apf::Numbering* n_new)
{
  const int maxval = getMaxNumber(m_new, n_new) + 1;
//  std::cout << "maxval = " << maxval << std::endl;
  apf::FieldShape* nshape = apf::getShape(n_new);
  const int ncomp = apf::countComponents(n_new);
  apf::MeshEntity* e;

  for (int dim = 0; dim < 3; ++dim)
  {
    // skip dimensions without nodes
    if ( !(nshape->hasNodesIn(dim)) )
    {
      continue;
    }

    // otherwise loop over them
    apf::MeshIterator* it = m_new->begin(dim);
    int entity_cnt = 0;
    while ( (e = m_new->iterate(it)) )  // loop over entities of this dimension
    {
      for (int node = 0; node < nshape->countNodesOn(m_new->simplexTypes[dim]); ++node)
      {
        for (int comp = 0; comp < ncomp; ++comp)
        {
          if (!apf::isNumbered(n_new, e, node, comp))
          {
//            std::cout << "assigning value to entity " << entity_cnt << " of dimension " << dim << " node " << node << " comp " << comp << std::endl;
            apf::number(n_new, e, node, comp, maxval);
          }  // end outer if
        }  // end loop over components
      }  // end loop over nodes
      ++entity_cnt;
    }  // end loop over entities
  }  // end loop over dimensions

}  // end function


// get the maximum number present in a numbering on the new mesh (safely)
int getMaxNumber(apf::Mesh* m_new, apf::Numbering* n_new)
{
//  std::cout << "entered getMaxNumber" << std::endl;
  apf::FieldShape* nshape = apf::getShape(n_new);
//  std::cout << "name of Numbering Fieldshape is: " << nshape->getName() << std::endl;
  const int ncomp = apf::countComponents(n_new);
  apf::MeshEntity* e;
  int val;
  int maxval = 0;

  for (int dim = 0; dim < 3; ++dim)
  {
    if ( !(nshape->hasNodesIn(dim)) )
    {
//      std::cout << "skipping dimensions " << dim << std::endl;
      continue;
    }


//    std::cout << "checking dimension " << dim << std::endl;
    int e_cnt = 0;
    apf::MeshIterator* it = m_new->begin(dim);
    while ( (e = m_new->iterate(it)) )  // loop over entities of this dimension
    {
//      std::cout << "  entity " << e_cnt << std::endl;
      for (int node = 0; node < nshape->countNodesOn(m_new->simplexTypes[dim]); ++node)
      {
//        std::cout << "    node " << node << std::endl;
        for (int comp = 0; comp < ncomp; ++comp)
        {
          if (apf::isNumbered(n_new, e, node, comp))
          {
            val = apf::getNumber(n_new, e, node, comp);
//            std::cout << "      has value = " << val << std::endl;
            if (val > maxval)
            {
              maxval = val;
            }  // end inner if
          }  // end outer if
        }  // end loop over components
      }  // end loop over nodes
      ++e_cnt;
    }  // end loop over entities
  }  // end loop over dimensions

//  std::cout << "  returning max value = " << maxval << std::endl;
  return maxval;

}  // end function


// figure out what the FieldShape on the new mesh should be
// The mapping is nonlinear, so this is an approximation
apf::FieldShape* getNewShape(apf::FieldShape* fshape_old)
{
  int nnodes_vert = fshape_old->countNodesOn(0);
  int nnodes_edge = fshape_old->countNodesOn(1);
  int nnodes_face = fshape_old->countNodesOn(2);

  // edge only shape
  if (nnodes_edge > 0 && nnodes_vert == 0 && nnodes_face == 0)
  {
    return apf::getConstant(1);
  } 
  // one node per face
  else if (nnodes_face == 1 && nnodes_edge == 0 && nnodes_vert == 0)
  {
    return apf::getConstant(2);
  } else
  {
    return apf::getLagrange(1);
  }

}

} // namespace triDG


// this function transferes the specified field to the new mesh
void transferFieldDG(apf::Mesh* m, apf::Mesh* m_new, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3], double* interp_op, apf::Field* field_old, apf::Field* field_new)
{

  // compute some quantities
  int nnodes_per_el = typeOffsetsPerElement[3] - 1;
  apf::FieldShape* mshape = apf::getShape(field_old);
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

  // transfer field values from the vertices of the old mesh to the new one
  triDG::transferVertices(m, m_new, numtriangles, triangulation, elementNodeOffsets, typeOffsetsPerElement, numberings, field_old, interp_op, field_new);

  // get lookup tables
  // this assumes there are no nodes on vertices of the old mesh
  int subelements[nnodes_per_el];
  int subelement_verts[nnodes_per_el];
  triDG::getFieldLookupTables(nnodes_per_el, triangulation, numtriangles, subelements, subelement_verts);

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

      // the node number on the old mesh that we are trying to transfer
      // to the new mesh
      nodenum = 0;

      for (int dim = 1; dim < 3; ++dim)  // loop all dimensions except verts
      {
        if (entity_nodes_on[dim] == 0)
        {
          continue;
        }

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
//          std::cout << "entity_node_idx = " << entity_node_idx << std::endl;
//          std::cout << "entity_idx = " << entity_idx << std::endl;
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

          triDG::copyFieldNode(field_old, field_new, node_entity, new_node_entity, offset_node_idx, 0, numcomp, vals);

          ++nodenum;
        }
      }

    } // end loop over old mesh elements
    
} // end function

void freeMemory(apf::MeshEntity*** verts, apf::MeshEntity*** edges, apf::MeshEntity*** faces, std::size_t entity_counts[3])
{
  apf::MeshEntity*** entities[3] = {verts, edges, faces};
  for (int dim = 0; dim < 3; ++dim)
  {
    apf::MeshEntity*** entities_dim = entities[dim];
    for (std::size_t i = 0; i < entity_counts[dim]; ++i)
    {
      free(entities_dim[i]);
    }
    free(entities_dim);
  }
}


// coords is a 2 x numNodesPerElement x numEl array of x,y coordinates of
// each node
// nodemap_ptos is the nodemap from the Pumi ordering to the SBP ordering ( 1 based)
apf::Mesh2* createSubMeshDG(apf::Mesh* m, apf::FieldShape* mshape, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], uint8_t* nodemap_ptos, apf::Numbering* numberings[3], double* coords_arr)
{
// m is the existing (high order) mesh
// numtriangles is the number of triangles to break each large triangle into
// triangulation is a numtriangles x 3 array holding the indices of the nodes, 
// 1-based 
//   to use as the vertices of the sub-triangles
//   numberings is an array holding a numbering for each entity dimension

//  return NULL;
  std::cout << "\ncreating sub-triangulated mesh" << std::endl;

//  apf::FieldShape* mshape_coords = m->getShape();
  // get number of different types of entities
  std::size_t entity_counts[3] = {m->count(0), m->count(1), m->count(2)};
  const int numEl = entity_counts[2];

  // get number of nodes on different types of entities
  int entity_nodes_on[3] = {mshape->countNodesOn(0), mshape->countNodesOn(1), mshape->countNodesOn(2)};
//  std::cout << "number of nodes on a face = " << entity_nodes_on[2] << std::endl;
  const int numNodesPerElement = entity_nodes_on[0] + entity_nodes_on[1] + entity_nodes_on[2];
  // create new mesh
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m_new = apf::makeEmptyMdsMesh(g, 2, false);

  // step 1: create all vertices of sub-triangles
   
  // allocate arrays to hold vertices of each order entity
  // these array should be in the Julia ordering, 
  // ie. getNumber(e_containing_edge) = first index
  std::cout << "allocating memory" << std::endl;
  apf::MeshEntity*** verts = (apf::MeshEntity***) calloc(entity_counts[0], sizeof(apf::MeshEntity**));  // 2 d array
  for (std::size_t i = 0; i < entity_counts[0]; ++i)
  {
    verts[i] = (apf::MeshEntity**)calloc(1, sizeof(apf::MeshEntity*));
  }
//  apf::MeshEntity* verts[entity_counts[0]][1]; // always create vertices, even
                                               // if no nodes
  apf::MeshEntity*** edges = (apf::MeshEntity***)calloc(entity_counts[1], sizeof(apf::MeshEntity**));
  for (std::size_t i = 0; i < entity_counts[1]; ++i)
  {
    edges[i] = (apf::MeshEntity**)calloc(entity_nodes_on[1], sizeof(apf::MeshEntity*));
  }
//  apf::MeshEntity* edges[entity_counts[1]][entity_nodes_on[1]];

  apf::MeshEntity*** faces = (apf::MeshEntity***)calloc(entity_counts[2], sizeof(apf::MeshEntity**));
  for (std::size_t i = 0; i < entity_counts[2]; ++i)
  {
    faces[i] = (apf::MeshEntity**)calloc(entity_nodes_on[2], sizeof(apf::MeshEntity*));
  }
//  apf::MeshEntity* faces[entity_counts[2]][entity_nodes_on[2]];

  std::cout << "finished allocating memory" << std::endl;

//  std::cout << "creating vertices on new mesh" << std::endl;
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
//       std::cout << "creating new vertex for old vertex " << idx << std::endl;
       m->getPoint(e, 0, coords);
       verts[idx][0] = m_new->createVert(0); // should be 2d array
       m_new->setPoint(verts[idx][0], 0, coords);
//       std::cout << "new vert pointer = " << verts[idx][0]; 
//       std::cout << " at coordinates " << coords[0] << " " << coords[1] << std::endl;
//       std::cout << "creating vertex for vertex " << idx << " at coordinates " << coords << " pointer = " << verts[idx][0] << std::endl;
      }

      // this probably doesn't work for DG meshes: there are no edge nodes
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

//        std::cout << "creating vertices on element " << idx << std::endl;
        for (int i = 0; i < entity_nodes_on[dim]; i++)
        {

//          std::cout << "  node " << i << std::endl;
          // get the coordinates
//          m->getPoint(e, i, coords);
          int coords_idx = triDG::getindex3c(0, nodemap_ptos[i] - 1, idx, 2,  numNodesPerElement, numEl);
//          std::cout << "component " << 0 << " node " << nodemap_ptos[i] - 1 << " element " << idx << " has offset " << coords_idx << std::endl;
          coords[0] = coords_arr[coords_idx];
          coords_idx = triDG::getindex3c(1, nodemap_ptos[i] - 1, idx, 2,  numNodesPerElement, numEl);
          coords[1] = coords_arr[coords_idx];

          faces[idx][i] = m_new->createVert(0); // should be 2d array
//          std::cout << "  new vert pointer = " << faces[idx][i] << " at coordinates " << coords[0] << " " << coords[1] << std::endl;
          m_new->setPoint(faces[idx][i], 0, coords);

//          std::cout << "creating edge vertex for face " << idx << ", node " << i << " at coordinates " << coords << " pointer = " << faces[idx][i] <<  std::endl;
        }
      }


/*
  std::cout << "Vert MeshEntity* :" << std::endl;
  triDG::printArray(std::cout, &verts[0][0], entity_counts[0], 1);

  std::cout << "Edge MeshEntity* :" << std::endl;
  triDG::printArray(std::cout, &edges[0][0], entity_counts[1], entity_nodes_on[1]);

  std::cout << "Face MeshEntity* :" << std::endl;
  triDG::printArray(std::cout, &faces[0][0], entity_counts[2], entity_nodes_on[2]);
*/

  // step 2: create elements from vertices
  // loop over elements on m, create element on m_new
  // here we use the edge orientation information
//  std::cout << "creating elements on new mesh" << std::endl;
  dim = 2;
  int elnum; // holds the element number
  int node;  // original node number
  int pos;  // linear address for elementNodeOffset
  const int nnodes_per_el = typeOffsetsPerElement[dim+1] - 1;
  int offset_j;
  apf::MeshEntity* el_verts[3]; // holds vertices of element to be created
  apf::MeshEntity* e_tmp;
  apf::Downward down;
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

        // handle verts as a special case because for DG meshes
        // there might not be nodes on vertices
        if ( node < N_VERT_PER_EL)
        {
//          std::cout << "this node is a vertex" << std::endl;
          const int dim = 0; // this is a vertex
          m->getDownward(e, dim, down);
          e_tmp = down[node];   // node = 1, 2, or 3, and is therefore the index
                            // of the vertex we want  
          int entity_num = apf::getNumber(numberings[dim], e_tmp, 0, 0);
      //    std::cout << "entity_num = " << entity_num << std::endl;

          el_verts[j] = verts[entity_num][0];
        } else {
//          std::cout << "this node is a regular node" << std::endl;
        //  not a vertex, follow normal procedure because there is definitely
        // a node on the vertex

          // calculate linear offset for elementNode offsets
          pos = node - N_VERT_PER_EL + elnum*nnodes_per_el;
//          std::cout << "offset pos = " << pos << std::endl;
          offset_j = elementNodeOffsets[pos];
          
//          std::cout << "offset = " << offset_j << std::endl;

          el_verts[j] = triDG::getVert(m, verts, edges, faces, typeOffsetsPerElement, node - N_VERT_PER_EL, offset_j, e, numberings, entity_nodes_on );
//          std::cout << "  el_vert " << j << " = " << el_verts[j] << std::endl;
        }
      }

//      std::cout << "creating element with verts " << el_verts[0] << " ";
//      std::cout << el_verts[1] << " " << el_verts[2] << std::endl;
     
      // build the element 
      apf::buildElement(m_new, 0, apf::Mesh::TRIANGLE, el_verts);

    }  // end loop over subtriangles

  } // end loop over elements
  
  freeMemory(verts, edges, faces, entity_counts);
  verts = NULL;
  edges = NULL;
  faces = NULL;

  // build, verify  mesh
  std::cout << "deriving model" << std::endl;
  apf::deriveMdsModel(m_new);
  std::cout << "finished deriving model" << std::endl;
  m_new->acceptChanges();
  std::cout << "accepted changes" << std::endl;
  m_new->verify();
  std::cout << "verified" << std::endl;

  triDG::transferNumberings(m, mshape, m_new, numtriangles, triangulation, elementNodeOffsets, typeOffsetsPerElement, numberings);

  // write visualization file
//  apf::writeVtkFiles("newmesh_linear", m_new);

  return m_new;
}



