// this file is the interface between julia and PUMI
// all functions defined in this file are usable using only standard C data 
// types
//
    // Every Julia process that calls this function will see its own "copy" of this
// library, having its own fully unique state
// this library is only tested to work on 3d meshes

// for now, this library is intended *only* to iterate over mesh entities in order
// functions (except those that loop over the type named in their declaration)
// do not increment the iterators


#include "funcs1.h"
//#include "a2.h"

//=============================================================================
//declare global variables (persistent state of library)

apf::Mesh2* m;
apf::FieldShape* mshape;
apf::MeshEntity* entity_global;  // token mesh entity used for EntityShape
apf::Numbering* elNums; // element numbering
apf::Numbering* faceNums; // face numbering
apf::Numbering* edgeNums; // edge numbering
apf::Numbering* vertNums; // vertex numbers

bool meshloaded = false;  // record whether or not a mesh has been loaded

apf::MeshEntity* e_tmp;

// arrays to hold variables for each entity type
// array[0] = vertex, array[1] = edge, array[2] = face, array[3] = element
apf::MeshIterator* its[4];  // mesh iterators
int numEntity[4];  // number of each type of entity
int numDown[4][4];  // number of downward adjacencies entity i has of type j
                 // define vertices to have zero
                 // this assumes only one type of element per mesh
apf::Numbering* numberings[4];  // numberings of each type of entity
                                // numberings are 1 index (not zero)
apf::MeshTag* globalVertNums;    // tag each vertex with global vertex number
const char *names[] = { "vertex", "edge", "face", "element"};  // array of strings, used for printing output inside loops
IsotropicFunctionJ isofunc;  // declare isotropic function at global scope
AnisotropicFunctionJ anisofunc; // declare anisotropic function at global scope


// init for 3d mesh
// order = order of shape functions to use
// load_mesh = load mesh from files or not (for reinitilizing after mesh adaptation, do not load from file)
// shape_type: type of shape functions, 0 =  lagrange, 1 = SBP
int initABC(char* dmg_name, char* smb_name, int number_entities[4], apf::Mesh2* m_ptr_array[1], apf::FieldShape* mshape_ptr_array[1], int order, int load_mesh, int shape_type )
{
//  std::cout << "Entered init\n" << std::endl;

  // various startup options

  // initilize communications if needed
  if (!PCU_Comm_Initialized())
  {
    MPI_Init(0,NULL);  // initilize MPI 
    PCU_Comm_Init();   // initilize PUMI's communication
  }
 

  if (load_mesh)  // if the user said to load a new mesh
  {
    if ( meshloaded)  // if a mesh has been loaded before
    {
      std::cout << "Performing cleanup before loading new mesh" << std::endl;
      cleanup(m);
    }

    if (strcmp(dmg_name, ".null") == 0)
    {
      gmi_register_null();
//      std::cout << "loading null geometric model" << std::endl;
      gmi_model* g = gmi_load(".null");
//      std::cout << "finished loading geometric model" << std::endl;
      m = apf::loadMdsMesh(g, smb_name);
//    apf::changeMeshShape(m, apf::getLagrange(2), true);
//    apf::changeMeshShape(m, apf::getLagrange(1), false); // for linear meshes
//    apf::changeMeshShape(m, apf::getSerendipity(), true);
//        apf::changeMeshShape(m, m->getShape(), false);
      m->verify();
    } else {
      gmi_register_mesh();
//      std::cout << "loading geometric model from file" << std::endl;
      m = apf::loadMdsMesh(dmg_name, smb_name);
    }

    meshloaded = true;  // record the fact that a mesh is now loaded
    apf::reorderMdsMesh(m);
    // apply shape functions to newly loaded mesh
    
  apf::FieldShape* fshape; // variable to hold field shape
  if ( shape_type == 0)
  {
    fshape = apf::getLagrange(order);
  } else if ( shape_type == 1)  // use SBP shape functions
  {
      fshape = apf::getSBP3Shape(order);

  } else  // default to lagrange shape functions
  {
    std::cout << "Warning: unrecognized shape_type, defaulting to Lagrange" << std::endl;
    fshape = apf::getLagrange(order);
  }



    if ( order == 1)
    {
      apf::changeMeshShape(m, fshape, false);
    } else
    {
      apf::changeMeshShape(m, fshape, true);
    }

//    std::cout << "finished loading mesh, changing shape" << std::endl;

  } else {  // if not loading a mesh
    destroyNumberings(3);  // destroy the numberings before creating new ones
//    std::cout << "finished destroying numberings of existing mesh" << std::endl;
  }

  // this should have a new filename everytime
//  apf::writeASCIIVtkFiles("output_check", m);


  m_ptr_array[0] = m;
  mshape_ptr_array[0] = m->getShape();
  its[0] = m->begin(0);
  entity_global = m->iterate(its[0]);  // get token mesh entity
  
  // initilize  number of each type of entity
  for (int i = 0; i < 4; ++i)
  {
    numEntity[i] = apf::countOwned(m, i);
    number_entities[i] = numEntity[i];
  }

  // initilize numberings
  numberings[0] = numberOwnedDimension(m, "vertNums", 0);
  numberings[1] = numberOwnedDimension(m, "edgeNums", 1);
  numberings[2] = numberOwnedDimension(m, "faceNums", 2);
  numberings[3] = numberOwnedDimension(m, "elNums", 3);



  // initilize iterators
  its[0] = m->begin(0);
  its[1] = m->begin(1);
  its[2] = m->begin(2);
  its[3] = m->begin(3);

  resetFaceIt();

//  apf::writeASCIIVtkFiles("output_init", m);
  return 0;
}






// init for 2d mesh or 3d mesh
// order = order of shape functions to use
// load_mesh = load mesh from files or not (for reinitilizing after mesh adaptation, do not load from file)
// PCU appears to not want a Communicator object?
int initABC2(const char* dmg_name, const char* smb_name, int number_entities[], apf::Mesh2* m_ptr_array[1], apf::FieldShape* mshape_ptr_array[1], int dim_ret[1], int order, int load_mesh, int shape_type )
{
//  std::cout << "Entered init2\n" << std::endl;

  // various startup options
  // initilize communications if needed
  int flag;
  MPI_Initialized(&flag);
  if (!flag)
  {
//    std::cout << "initializing MPI" << std::endl;
    MPI_Init(0,NULL);  // initilize MPI 
  }
  if (!PCU_Comm_Initialized())
  {
//    std::cout << "initializing PCU" << std::endl;
    PCU_Comm_Init();   // initilize PUMI's communication
  }
 
  int dim;
  if (load_mesh)  // if the user said to load a new mesh
  {
    if ( meshloaded)  // if a mesh has been loaded before
    {
      std::cout << "Performing cleanup before loading new mesh" << std::endl;
      cleanup(m);
    }

    if (strcmp(dmg_name, ".null") == 0)
    {
      gmi_register_null();
//      std::cout << "loading null geometric model" << std::endl;
      gmi_model* g = gmi_load(".null");
//      std::cout << "finished loading geometric model" << std::endl;
      m = apf::loadMdsMesh(g, smb_name);
//    apf::changeMeshShape(m, apf::getLagrange(2), true);
//    apf::changeMeshShape(m, apf::getLagrange(1), false); // for linear meshes
//    apf::changeMeshShape(m, apf::getSerendipity(), true);
//        apf::changeMeshShape(m, m->getShape(), false);
      m->verify();
    } else {
      gmi_register_mesh();
//      std::cout << "loading geometric model from file" << std::endl;
      m = apf::loadMdsMesh(dmg_name, smb_name);
    }

    dim = m->getDimension();
    meshloaded = true;  // record the fact that a mesh is now loaded
    apf::reorderMdsMesh(m);

    int order_orig = m->getShape()->getOrder();
    bool change_shape;
    apf::FieldShape* fshape = getFieldShape(shape_type, order, dim, change_shape);
    // check if the coordinates have been moved into tags
    // this is a result of the mesh shape having been changed previously
    apf::MeshTag* coords_tag = m->findTag("coordinates_ver");
    
    // if the tag exists and if original shape is linear
    // this is a workaround for old meshes where the coordinate field was
    // not saved correctly.  Presumably they have linear fields only, so
    // only force the tag data back into the regular field for linear meshes.
    // For higher order meshes, keep the data.
    if (coords_tag != 0 && order_orig == 1)
    {
      apf::changeMeshShape(m, apf::getLagrange(1), false);
      std::cout << "finished first mesh shape change" << std::endl;
    }
    
    if ( change_shape && order_orig > order)
    {
      std::cerr << "Warning: changing mesh coordinate field from " << order_orig << " to " << order << " will result in loss of resolution" << std::endl;
    }

    if (change_shape)
    {
      std::cout << "about to perform final mesh shape change" << std::endl;
      if ( order == 1 )
      {
        apf::changeMeshShape(m, fshape, true);
      } else
      {
        apf::changeMeshShape(m, fshape, true);
      }

//      std::cout << "finished loading mesh, changing shape" << std::endl;
//      std::cout << "new shape name is " << m->getShape()->getName() << std::endl;
    } else
    {
//      std::cout << "finished loading mesh" << std::endl;
    }

  } else {  // if not loading a mesh
    dim = m->getDimension();
    destroyNumberings(2);  // destroy the numberings before creating new ones
  }

  dim_ret[0] = dim; // return to julia

  // this should have a new filename everytime
//  apf::writeASCIIVtkFiles("output_check", m);


  m_ptr_array[0] = m;
  mshape_ptr_array[0] = m->getShape();
  its[0] = m->begin(0);
  entity_global = m->iterate(its[0]);  // get token mesh entity
  
  // initilize  number of each type of entity
  for (int i = 0; i < (dim+1); ++i)
  {
    numEntity[i] = (m->count(i));
//      apf::countOwned(m, i);
    number_entities[i] = numEntity[i];
  }

  // create numberings
  numberings[0] = apf::createNumbering(m, "vertNums", apf::getConstant(0), 1);
  numberings[1] = apf::createNumbering(m, "edgeNums", apf::getConstant(1), 1);
  numberings[2] = apf::createNumbering(m, "faceNums", apf::getConstant(2), 1);
  if (dim == 3)
    numberings[3] = apf::createNumbering(m, "regionNums", apf::getConstant(3), 1);




  // initilize iterators
  its[0] = m->begin(0);
  its[1] = m->begin(1);
  its[2] = m->begin(2);
  if (dim == 3)
    its[3] = m->begin(3);


  for (int i = 0; i < (dim+1); ++i)  // loop over dimensions
  {
    int curr_num = 0;
    apf::MeshEntity* e_local;
    while ( (e_local = m->iterate(its[i])) )
    {
//      std::cout << "curr_num = " << curr_num << std::endl;
      apf::number(numberings[i], e_local, 0, 0, curr_num);
      ++curr_num;
    }
    its[i] = m->begin(i);  // reset iterator
  }
   
  resetFaceIt();

  apf::writeASCIIVtkFiles("output_init", m);
/*
  // write curved mesh visualization file
  apf::FieldShape* mshape = m->getShape();
  crv::writeControlPointVtuFiles(m, "houtput");
  crv::writeCurvedVtuFiles(m, apf::Mesh::EDGE, mshape->countNodesOn(apf::Mesh::EDGE), "houtput");
  crv::writeCurvedVtuFiles(m, apf::Mesh::TRIANGLE, mshape->countNodesOn(apf::Mesh::TRIANGLE), "houtput");
//  ma::writePointSet(m, 2, 21, "pointcloud");
*/
  return 0;
}

// perform cleanup activities, making it safe to load a new mesh
void cleanup(apf::Mesh* m_local)
{
  m_local->destroyNative();
  apf::destroyMesh(m_local);
  meshloaded = false;
//  PCU_Comm_Free();
//  MPI_Finalize();
}

// destroy numberings after mesh adaptation
// dim specifies dimension of existing mesh (2 or 3d)
void destroyNumberings(int dim)
{
  apf::destroyNumbering(numberings[0]);
  apf::destroyNumbering(numberings[1]);
  apf::destroyNumbering(numberings[2]);

  if (dim > 2)
  {
   apf::destroyNumbering(numberings[2]);
  }
 
}

// this function defines the mapping from the shape_type integer to the
// FieldShape itself
apf::FieldShape* getFieldShape(int shape_type, int order, int dim, bool& change_shape)
{
//  std::cout << "dimension " << dim << " mesh requests shape type " << shape_type << " of order " << order << std::endl;
  apf::FieldShape* fshape;
  if (dim == 2)
  {
    if ( shape_type == -1)  // use existing
    {
      fshape = NULL;
      change_shape = false;
    }
    if ( shape_type == 0)  // use lagrange
    {
      fshape = apf::getLagrange(order);
      change_shape = true;
    } else if ( shape_type == 1)  // use SBP shape functions
    {
        fshape = apf::getSBPShape(order);
        change_shape = true;
    } else if ( shape_type == 2)  // use SBP DG1 shape functions
    {
      fshape = apf::getDG1SBPShape(order);
      change_shape = true;
    } else if ( shape_type == 3) // use SBP DG2 shape functions
    {
      fshape = apf::getDG2SBPShape(order);
      change_shape = true;
    } else  // default to lagrange shape functions
    {
      std::cout << "Warning: unrecognized shape_type, not changing mesh shape" << std::endl;
      fshape = apf::getLagrange(1); // unused, but avoids compiler warning
      change_shape = false;
    }
  } else if (dim == 3)
  {
    if (shape_type == -1)
    {
      fshape = NULL;
      change_shape = false;
    }
    if (shape_type == 0) // use lagrange
    {
      fshape = apf::getLagrange(order);
      change_shape = true;
    } else if (shape_type == 2)  // use SBP DG1 shape functions
    {
      fshape = apf::getDG1SBP3Shape(order);
      change_shape = true;
    } else if (shape_type == 3)
    {
      fshape = apf::getDG2SBP3Shape(order);
      change_shape = true;
    } else  // default to lagrange shape functions
    {
      std::cout << "Warning: unrecognizes shape_type, not changing mesh shape" << std::endl;
      fshape = apf::getLagrange(1);
      change_shape = false;
    }
  } else
  {
    std::cout << "Warning: unrecognized dimension, not changing mesh shape" << std::endl;
    fshape = apf::getLagrange(1);
    change_shape = false;
  }

//  std::cout << "returning FieldShape " << fshape->getName() << std::endl;
     
  return fshape;
}

apf::Mesh2* getMeshPtr()
{
  return m;
}

apf::FieldShape* getMeshShapePtr(apf::Mesh* m)
{
  return m->getShape();
}

// get constant field shape of given dimension
apf::FieldShape* getConstantShapePtr(int dimension)
{
  return apf::getConstant(dimension);
}



apf::Numbering* getVertNumbering()
{
  return numberings[0];
}

apf::Numbering* getEdgeNumbering()
{
  return numberings[1];
}

apf::Numbering* getFaceNumbering()
{
  return numberings[2];
}

apf::Numbering* getElNumbering()
{
  return numberings[3];
}

void resetVertIt()
{
  its[0] = m->begin(0);
}

void resetEdgeIt()
{
  its[1] = m->begin(1);
}

void resetFaceIt()
{
  its[2] = m->begin(2);
}

void resetElIt()
{
  its[3] = m->begin(3);
}

void resetIt(int dim)
{
  its[dim] = m->begin(dim);
}

void incrementVertIt()
{
  m->iterate(its[0]);
}

void incrementIt(int dim)
{
  m->iterate(its[dim]);
}

// increment vertex iterator n times
void incrementVertItn(int n)
{
  for (int i=0; i < n; ++i)
    m->iterate(its[0]);
}

void incrementEdgeIt()
{
  m->iterate(its[1]);
}

// increment edge iterator n times
void incrementEdgeItn(int n)
{
  for (int i=0; i < n; ++i)
    m->iterate(its[1]);
}


void incrementFaceIt()
{
  m->iterate(its[2]);
}


// increment face iterator n times
void incrementFaceItn(int n)
{
  for (int i=0; i < n; ++i)
    m->iterate(its[2]);
}

void incrementElIt()
{
  m->iterate(its[3]);
}


// increment element iterator n times
void incrementElItn(int n)
{
  for (int i=0; i < n; ++i)
    m->iterate(its[3]);
}


int count(apf::Mesh2* m_local, int dimension)
{
  return m_local->count(dimension);
}

void writeVtkFiles(char* name, apf::Mesh2* m_local)
{

//   apf::writeASCIIVtkFiles(name, m_local);
    
   apf::writeVtkFiles(name, m_local);
/*
    apf::FieldShape* mshape = m->getShape();
    crv::writeControlPointVtuFiles(m, name);
    crv::writeCurvedVtuFiles(m, apf::Mesh::EDGE, mshape->countNodesOn(apf::Mesh::EDGE), name);
    crv::writeCurvedVtuFiles(m, apf::Mesh::TRIANGLE, mshape->countNodesOn(apf::Mesh::TRIANGLE), name);

  }
  */
}


// tag current vertex with val, part of MeshTag globalNodeNumbers
//
// these functions no longer work
void setGlobalVertNumber(int val)
{
  apf::MeshEntity* e = m->deref(its[0]);  // get current vertex
  // tag the vertex with a value
  // setIntTag expects an array, so pass it a pointer to a single integer
  m->setIntTag(e, globalVertNums, &val);  // tag the vertex with the a value
}
 
// get the global vertex number of the current vertex
int getGlobalVertNumber()
{
  apf::MeshEntity* e = m->deref(its[0]);
  int val;
  m->getIntTag(e, globalVertNums, &val);
  return val;
}


int getVertNumber()
{
  apf::MeshEntity*e = m->deref(its[0]);
  int num = apf::getNumber(numberings[0], e, 0, 0);
  return num;
}

int getEdgeNumber()
{
  apf::MeshEntity*e = m->deref(its[1]);
  int num = apf::getNumber(numberings[1], e, 0, 0);
  return num;
}

int getFaceNumber()
{
  apf::MeshEntity*e = m->deref(its[2]);
  int num = apf::getNumber(numberings[2], e, 0, 0);
  return num;
}


int getElNumber()
{
  apf::MeshEntity*e = m->deref(its[3]);
  int num = apf::getNumber(numberings[3], e, 0, 0);
  return num;
}

// return MeshEntity pointer to current vertex
apf::MeshEntity* getVert() 
{
  apf::MeshEntity* e = m->deref(its[0]);
  return e;
}

// return MeshEntity pointer to current edge
apf::MeshEntity* getEdge() 
{
  apf::MeshEntity* e = m->deref(its[1]);
  return e;
}


// return MeshEntity pointer to current face
apf::MeshEntity* getFace() 
{
  apf::MeshEntity* e = m->deref(its[2]);
  return e;
}


// return MeshEntity pointer to current element
apf::MeshEntity* getEl() 
{
  apf::MeshEntity* e = m->deref(its[3]);
  return e;
}

apf::MeshEntity* getEntity(int dim)
{
  apf::MeshEntity* e = m->deref(its[dim]);
  return e;
}

// get number of the vertex MeshEntity* e
int getVertNumber2(apf::MeshEntity* e)
{
  int i = apf::getNumber(numberings[0], e, 0 , 0);
  return i;
}

// get number of the edge MeshEntity* e
int getEdgeNumber2(apf::MeshEntity* e)
{
  int i = apf::getNumber(numberings[1], e, 0 , 0);
  return i;
}


// get number of the face MeshEntity* e
int getFaceNumber2(apf::MeshEntity* e)
{
  int i = apf::getNumber(numberings[2], e, 0 , 0);
  return i;
}

// get number of the element MeshEntity* e
int getElNumber2(apf::MeshEntity* e)
{
  int i = apf::getNumber(numberings[3], e, 0 , 0);
  return i;
}


// get model info of a mesh element
apf::ModelEntity* toModel(apf::Mesh* m_local, apf::MeshEntity* e)
{
  return m_local->toModel(e);
}

// get the *dimension* of the model entity
// not its type like the name implies
int getModelType(apf::Mesh* m_local, apf::ModelEntity* e)
{
  return m_local->getModelType(e);
}

int getModelTag(apf::Mesh* m_local, apf::ModelEntity* e)
{
  return m_local->getModelTag( e);
}




// get dimension of mesh
int getMeshDimension(apf::Mesh2* m_local)
{
  int i = m_local->getDimension();
//  std::cout << "Mesh Dimesion = " << i << std::endl;
  return i;
}

// untested
int getType(apf::Mesh2* m_local, apf::MeshEntity* e)
{
  return m_local->getType(e);
}

// get the downward adjacencies of a MeshEntity
// supply a 12 element array of Ptr{Void}
// the function return an integer telling the number of entities put into the
// array
int getDownward(apf::Mesh2* m_local, apf::MeshEntity* e, int dimension, apf::MeshEntity* downwards[12])
{
//  std::cout << "in C++ getDownwards, dimension = " << dimension << std::endl;
  int numDown = m->getDownward(e, dimension, downwards);
  return numDown;
}


apf::DynamicArray<apf::MeshEntity*> adjacencies;
bool adjacent_ready;  // stores whether countAdjacent has been called
                              // before getAdjacent
// gets the mesh entities adjacenct to the given meshentity e, stores them in
// the global variable, and returns the number of entities
// the getAdjacent  function only returns the adjcencies fetched by this function
int countAdjacent(apf::Mesh2* m_local, apf::MeshEntity* e, int dimension)
{
  m_local->getAdjacent(e, dimension, adjacencies);
  int num_adjacent = adjacencies.getSize();
  adjacent_ready = true;
//  std::cout << " set adjcent_read = true" << std::endl;
  return num_adjacent;
}

// returns the adjacencies fetched by countAdjacent, copying them into
// adjacencies_ret[]
void getAdjacent(apf::MeshEntity* adjacencies_ret[])
{
  if (!adjacent_ready)
  {
    std::cout << "Warning, adjacencies might not  be ready" << std::endl;
  }

  int num_adjacent = adjacencies.getSize();
  for (int i = 0; i < num_adjacent; ++i)
  {
    adjacencies_ret[i] = adjacencies[i];
  }

  adjacent_ready = false;
 // std::cout << "set adjacent_ready = false" << std::endl;

}


// second order adjacencies, similar to countAdjacent/getAdjacent
apf::DynamicArray<apf::MeshEntity*> adjacencies2;
bool adjacent_ready2 = false;

int countBridgeAdjacent(apf::Mesh* m_local, apf::MeshEntity* origin, int bridge_dimension, int target_dimension)
{
  apf::getBridgeAdjacent(m_local, origin, bridge_dimension, target_dimension, adjacencies2);

  int num_adjacent2 = adjacencies2.getSize();
  adjacent_ready2 = true;
  return num_adjacent2;
}


void getBridgeAdjacent(apf::MeshEntity* adjacencies_ret2[])
{

  if (!adjacent_ready2)
  {
    std::cout << "Warning, adjacencies might not be ready" << std::endl;
  }

  int num_adjacent2 = adjacencies2.getSize();
  for (int i = 0; i < num_adjacent2; ++i)
  {
    adjacencies_ret2[i] = adjacencies2[i];
  }

  adjacent_ready2 = false;

}



void getAlignment(apf::Mesh* m_local, apf::MeshEntity* elem, apf::MeshEntity* boundary, int which[1], bool flip[1], int rotate[1])
{
//  std::cout << "in getAlignment, which = " << which << " flip = " << flip << " rotate = " << rotate << std::endl;
  apf::getAlignment(m_local, elem, boundary, which[0], flip[0], rotate[0]);

//  std::cout << "after apf::getAlignment, which = " << which << " flip = " << flip << " rotate = " << rotate << std::endl;
}



// apf::Fieldshape Functions
// check whether the given field as nodes on entities of the specified dimension
bool hasNodesIn(apf::FieldShape* fshape_local, int dimension)
{
  return fshape_local->hasNodesIn(dimension);
}

// determine the number of nodes on entities of a given type (apf::Mesh::Type)
int countNodesOn(apf::FieldShape* fshape_local, int type)
{
  return fshape_local->countNodesOn(type);
}

// get the EntityShape (object describing shape functions) of a type of entity
apf::EntityShape* getEntityShape(apf::FieldShape* mshape_local, int type)
{
  return mshape_local->getEntityShape(type);
}

int getOrder(apf::FieldShape* fshape)
{
  return fshape->getOrder();
}


// MeshElement related functions
extern apf::MeshElement* createMeshElement( apf::Mesh2* m_local, apf::MeshEntity* e)
{
  return apf::createMeshElement(m_local, e);
}

// count the number of integration points needed to achieve specified order of 
// accuracy.  MeshElement can be edges as well as 2D and 3D regions.
// not sure what happens if it is a vertex.
int countIntPoints(apf::MeshElement* e, int order)
{
  return apf::countIntPoints(e, order);
}

extern void getIntPoint(apf::MeshElement* e, int order, int point, double coords[3])
{
  apf::Vector3 vec;
  apf::getIntPoint(e, order, point, vec); // coordinates are now in vec
  coords[0] = vec[0];
  coords[1] = vec[1];
  coords[2] = vec[2];
}


extern double getIntWeight(apf::MeshElement* e, int order, int point)
{
  return apf::getIntWeight(e, order, point);
}

// gets the jacobian of a mesh element at a specified location in parametric coordinates
extern void getJacobian(apf::MeshElement* e, double coords[3], double jac[3][3])
{
  // copy coordinates to vec
  apf::Vector3 vec (coords[0], coords[1], coords[2]);

  // create matrix to hold jacobian
  apf::Matrix3x3 jac_local;

  // poppulate jac_local with the jacobian
  apf::getJacobian(e, vec, jac_local);

  // copy the jacobian to jac to be returned
  jac_local.toArray(jac);

//  std::cout << "jac_local = " << jac_local << std::endl;
}


// functions involving EntityShape

int countNodes(apf::EntityShape* eshape_local)
{
  return eshape_local->countNodes();
}

// get shape function values
// vals had better be the right size
void getValues(apf::EntityShape* eshape_local, double xi[3], double vals[])
{
  apf::Vector3 xi_vec (xi[0], xi[1], xi[2]); // create vector of coordinates
  apf::NewArray<double> vals_array;  // array to store retrieved values in
  int numN_local = eshape_local->countNodes();  // get number of nodes
  eshape_local->getValues(m, entity_global, xi_vec, vals_array);

  // copy vals_array into vals to be returned
  for (int i = 0; i < numN_local; ++i)
  {
    vals[i] = vals_array[i];
  }

}

// get shape function gradient values
// vals had better be the right size
void getLocalGradients(apf::EntityShape* eshape_local, double xi[3], double vals[][3])
{
  apf::Vector3 xi_vec (xi[0], xi[1], xi[2]); // create vector of coordinates
  apf::NewArray<apf::Vector3> vals_array;  // array to store retrieved values in
  int numN_local = eshape_local->countNodes();  // get number of nodes
  eshape_local->getLocalGradients(m, entity_global, xi_vec, vals_array);

  // copy vals_array into vals to be returned
  for (int i = 0; i < numN_local; ++i)
  {
    apf::Vector3 tmp = vals_array[i];  // get the ith vector
    for (int j = 0; j < 3; ++j)
    {
      vals[i][j] = tmp[j];  // copy the vector into a column of vals
    }
  }

}

// get the array that transforms from local element order to canonical order
// the length of order must be the number of nodes classified on
// apf::MeshEntity* shared
void alignSharedNodes(apf::EntityShape* eshape_local, apf::Mesh* m_local, apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
{
  eshape_local->alignSharedNodes(m_local, elem, shared, order);
}

// check that global variable are persisent
// doesn't work in 2d
void checkVars()
{
  std::cout << "Entered checkVars()" << std::endl;
  for (int i = 0; i < 4; ++i)
  {
    int type = m->getType(m->deref(its[i]));
    std::cout << "type of its[" << i << "] = " << type;
    int index = apf::getMdsIndex(m, m->deref(its[i]));
    std::cout << "index its[" << i << "] = " << index << std::endl;
    m->iterate(its[i]);
  }

  std::cout << "iterated all iterators stored in array" << std::endl;
  for (int i = 0; i < 4; ++i)
  {
    int index = apf::getMdsIndex(m, m->deref(its[i]));
    std::cout << "index its[" << i << "] = " << index << std::endl;

  }

}


// print numberings of all mesh entitites
// doesn't work in 2d
void checkNums()
{
  std::cout << "Entered checkNums" << std::endl;
  apf::MeshEntity* e;
  int i = 0;  // counter
  int id = 0; // store numbering output
  std::cout << "about to start loop" << std::endl;
  apf::MeshIterator* it = m->begin(0);
  for (int j = 0; j < 4; ++j)  // loop over mesh entity types
  {
    it = m->begin(j);
    while ( ( e = m->iterate(it) ) ) 
//    while (( e = m->iterate(its[j]) ))  // loop over entities of type j
    {
 //       int dim = m->getType(e);
//        std::cout << "type of current entity = " << dim << std::endl;
  //    std::cout << "About to fetch numbering" << std::endl;
      id = apf::getNumber(numberings[j], e, 0, 0);
      std::cout << names[j] << " number " << i << " has number " << id << std::endl;
  //    std::cout << "i = " << i << std::endl;
       ++i;
    }
    // reset for next outer loop
    std::cout << std::endl;
    i = 0;
  }

  // reset iterators
  resetVertIt();
  resetEdgeIt();
  resetFaceIt();
  resetElIt();

}

// get the coordinates of a vertex
// coords is 2d array to populate, sx and sy are dimension of the array
void getVertCoords(double coords[][3], int sx, int sy)
{
//  std::cout << "in getVertCoords" << std::endl;
//  std::cout << "sx = " << sx << " ,sy = " << sy << std::endl;
  apf::MeshEntity* e = m->deref(its[0]);
  apf::Vector3 vec;
  m->getPoint(e, 0, vec);
//  std::cout << "coords = " << vec << std::endl;
  coords[0][0] = vec[0];
  coords[0][1] = vec[1];
  coords[0][2] = vec[2];

}

// get the coordinates of a vertex
// coords is 2d array to populate, sx and sy are dimension of the array
void getVertCoords2(apf::MeshEntity* e, double coords[][3], int sx, int sy)
{
//  std::cout << "in getVertCoords" << std::endl;
//  std::cout << "sx = " << sx << " ,sy = " << sy << std::endl;
  apf::Vector3 vec;
  m->getPoint(e, 0, vec);
//  std::cout << "coords = " << vec << std::endl;
  coords[0][0] = vec[0];
  coords[0][1] = vec[1];
  coords[0][2] = vec[2];

}

// get the coordinates of the two points that define, iterate edge iterator
// an edge
int getEdgeCoords(double coords[2][3], int sx, int sy)
{
//  std::cout <<"in getEdgeCoords" << std::endl;
  if (sx < 2 || sy != 3)
  {
    std::cout << "Warning: coords array too small" << std::endl;
    return -1;
  }

  apf::MeshEntity* e = m->deref(its[1]);

  apf::Vector3 vec;
  apf::MeshEntity* verts[2];
//  std::cout << " about to get downward entities" << std::endl;
  m->getDownward(e, 0, verts);
//  std::cout << "got vertices of edge" << std::endl;
  // get coordinates of each vertex
  m->getPoint(verts[0], 0, vec);
  coords[0][0] = vec[0];
  coords[0][1] = vec[1];
  coords[0][2] = vec[2];

//  std::cout << "first point coordinates = " << vec << std::endl;

  m->getPoint(verts[1], 0, vec);
  coords[1][0] = vec[0];
  coords[1][1] = vec[1];
  coords[1][2] = vec[2];

//  std::cout << "second point coordinates = " << vec << std::endl;

  return 0;
}

// get the coordinates of the two points that define and edge
int getEdgeCoords2(apf::MeshEntity* e, double coords[2][3], int sx, int sy)
{
//  std::cout <<"in getEdgeCoords" << std::endl;
  if (sx < 2 || sy != 3)
  {
    std::cout << "Warning: coords array too small" << std::endl;
    return -1;
  }

  apf::Vector3 vec;
  apf::MeshEntity* verts[2];
//  std::cout << " about to get downward entities" << std::endl;
  m->getDownward(e, 0, verts);
//  std::cout << "got vertices of edge" << std::endl;
  // get coordinates of each vertex
  m->getPoint(verts[0], 0, vec);
  coords[0][0] = vec[0];
  coords[0][1] = vec[1];
  coords[0][2] = vec[2];

//  std::cout << "first point coordinates = " << vec << std::endl;

  m->getPoint(verts[1], 0, vec);
  coords[1][0] = vec[0];
  coords[1][1] = vec[1];
  coords[1][2] = vec[2];

//  std::cout << "second point coordinates = " << vec << std::endl;

  return 0;
}


// get populate coords with coordinates of vertices that define the current face
//  and iterates of face iterator
// sx = x dimension of array, sy = y dimension of array
int getFaceCoords(double coords[][3], int sx, int sy)
{
  std::cout << "Entered getFaceCoords" << std::endl;
  apf::MeshEntity* e = m->deref(its[2]);
  apf::Downward verts;  // hold vertices
  int numDownward = m->getDownward(e, 0, verts);  // populate verts
  std::cout << "number of downward verts = " << numDownward << std::endl;
  if ( sx < numDownward || sy != 3)
  {
    std::cout << " Warning, array not right size ";
    std::cout << " Array must be " << numDownward << " by 3, is ";
    std::cout << sx << " by " << sy << std::endl;
    return -1;
  }

  apf::Vector3 vec;  // vector to hold coordinates of a point
  for ( int i = 0; i < numDownward; ++i) // loop over vertices
  {
    std::cout << "getting coordinates for vertex " << i << " of face" << std::endl;
    m->getPoint(verts[i], 0, vec); // populate vec with coordinates of vertex
//    std::cout << "coordinates of vertex " << i << " = " << vec << std::endl;
    coords[i][0] = vec[0];
    coords[i][1] = vec[1];
    coords[i][2] = vec[2];
  }

  return 0;
}

// get populate coords with coordinates of vertices that define the current face
// sx = x dimension of array, sy = y dimension of array
int getFaceCoords2(apf::MeshEntity* e, double coords[][3], int sx, int sy)
{
//  std::cout << "Entered getFaceCoords22" << std::endl;
//  std::cout << " element = " << e << std::endl;
  apf::Downward verts;  // hold vertices
  int numDownward = m->getDownward(e, 0, verts);  // populate verts

//  std::cout << "number of downward verts = " << numDownward << std::endl;
  if ( sx < numDownward || sy != 3)
  {
    std::cout << " Warning, array not right size ";
    std::cout << " Array must be " << numDownward << " by 3, is ";
    std::cout << sx << " by " << sy << std::endl;
    return -1;
  }

  apf::Vector3 vec;  // vector to hold coordinates of a point
  for ( int i = 0; i < numDownward; ++i) // loop over vertices
  {

//    std::cout << "getting coordinates for vertex " << i << " of face" << std::endl;
    m->getPoint(verts[i], 0, vec); // populate vec with coordinates of vertex
//    std::cout << "coordinates of vertex " << i << " = " << vec << std::endl;
    coords[i][0] = vec[0];
    coords[i][1] = vec[1];
    coords[i][2] = vec[2];
  }

  return 0;
}




int getElCoords(double coords[][3], int sx, int sy)
{
//  std::cout << "Entered getElCoords" << std::endl;

  apf::MeshEntity* e = m->deref(its[3]);
  apf::Downward verts; // hold vertices
  int numDownward = m->getDownward(e, 0, verts); // populate verts

  if ( sx < numDownward || sy != 3)
  {
    std::cout << " Warning, array not right size ";
    std::cout << " Array must be " << numDownward << " by 3, is ";
    std::cout << sx << " by " << sy << std::endl;
    return -1;
  }

  apf::Vector3 vec;  // vector to hold coordinates of a point
  for ( int i = 0; i < numDownward; ++i) // loop over vertices
  {
    m->getPoint(verts[i], 0, vec); // populate vec with coordinates of vertex
//    std::cout << "coordinates of vertex " << i << " = " << vec << std::endl;
    coords[i][0] = vec[0];
    coords[i][1] = vec[1];
    coords[i][2] = vec[2];
  }

  return 0;
}



int getElCoords2(apf::MeshEntity* e, double coords[][3], int sx, int sy)
{
  apf::Downward verts; // hold vertices
  int numDownward = m->getDownward(e, 0, verts); // populate verts
//  std::cout << "entity e is of type " << m->getType(e) << std::endl;
  if ( sx < numDownward || sy != 3)
  {
    std::cout << " Warning, array not right size ";
    std::cout << " Array must be " << numDownward << " by 3, is ";
    std::cout << sx << " by " << sy << std::endl;
    return -1;
  }

  apf::Vector3 vec;  // vector to hold coordinates of a point
  for ( int i = 0; i < numDownward; ++i) // loop over vertices
  {
    m->getPoint(verts[i], 0, vec); // populate vec with coordinates of vertex
//    std::cout << "coordinates of vertex " << i << " = " << vec << std::endl;
    coords[i][0] = vec[0];
    coords[i][1] = vec[1];
    coords[i][2] = vec[2];
  }

  return 0;
}

// this function gets the coordinates of all nodes on a mesh entity (including
// all verts, edges, faces, and regions) and puts them in the array coords
// Coords should be a dim x numNodesPerElement (in the FieldShape of the
// coordinate field).
// The nodes are ordered the same as the downward adjacencies, from lowest
// dimensions to highest (ie. verts, edges, faces, regions)
void getAllEntityCoords(apf::Mesh* m, apf::MeshEntity* e, double* coords)
{

  apf::FieldShape* coordshape = m->getShape();
  int max_dim = apf::Mesh::typeDimension[m->getType(e)];


  int nentities=0;
  apf::Downward down_entities;
  apf::Vector3 coords_vec;
  apf::MeshEntity* down_entity;

  int idx=0;
  for ( int dim=0; dim <= max_dim; ++dim)
  {
    if ( coordshape->hasNodesIn(dim) )
    {
      nentities = m->getDownward(e, dim, down_entities);
      // loop over entities of current dimensions
      for ( int entitynum=0; entitynum < nentities; entitynum++)
      {
        down_entity = down_entities[entitynum];
        int entity_type = m->getType(down_entity);
        int nnodes = coordshape->countNodesOn(entity_type);
        // loop over nodes on current entity
        for (int nodenum=0; nodenum < nnodes; nodenum++)
        {
          m->getPoint(down_entity, nodenum, coords_vec);
          for (int i=0; i < max_dim; ++i)
          {
            coords[idx] = coords_vec[i];
            idx++;
          }  // end loop i
        }  // end loop over nodes
      }  // end loop over entities
    } // end if nodes in this dimension
  }   // end loop over dimensions
}  // function getAllElementCoords


//TODO: create a similar function to the above calling getNodeXi
// look at apf.h for vertex xi coordinates


// create a generally defined numbering from julia
apf::Numbering* createNumberingJ(apf::Mesh2* m_local, char* name, apf::FieldShape* field, int components)
{
//    return apf::createNumbering(m, "num1", m->getShape(), 1);
//    return apf::createNumbering(m_local, "num1", m->getShape(), 1);
//    return apf::createNumbering(m_local, "num1", m_local->getShape(), components);
    return apf::createNumbering(m_local, name, field, components);
}

apf::FieldShape* getNumberingShape(apf::Numbering* n)
{
  return apf::getShape(n);
}

// number an entity in a given numbering from julia
int numberJ(apf::Numbering* n, apf::MeshEntity* e, int node, int component, int number)
{
  apf::number(n, e, node, component, number);
//  int i = apf::getNumber(n,e,node,component);
  return 0;
}

// retrieve a number from julia
int getNumberJ(apf::Numbering* n, apf::MeshEntity* e, int node, int component)
{
  int i = apf::getNumber(n, e, node, component);
  return i;
}


// this is a non elementary function, here for performance reasons only
// gets the dof numbers for all nodes affected by the element, in order
// n is the dof numbering
// entities is an array of node entities that have nodes on the elements, including repeats, in the Pumi order
// element is the element
// dofnums is the output array, in the SBP order
// the output array will be seen by julia as being num_comp by length(entities)
// nodemap is typically the map from the Pumi node order to the desired
// ordering, 1 based indexing
// collect dofs in Pumi ordering, write them to array in SBP order
int getDofNumbers(apf::Numbering* n, apf::MeshEntity* entities[], uint8_t node_offsets[], uint8_t nodemap[], apf::MeshEntity* element, int dofnums[])
{
//  std::cout << "Entered getDofNumbers" << std::endl;
  // declare and initialize static variables
  static apf::Mesh* m_local = apf::getMesh(n); 
  static int el_type = m_local->getType(element);
  static int el_dim = m_local->typeDimension[el_type];
  static int num_comp = apf::countComponents(n);
  static apf::FieldShape* fshape_local = apf::getShape(n);
  static apf::MeshEntity* e = entities[0];
  static int col = 0;  // current node
  static int ptr = 0;  // pointer to linear address in entities

  // populate them
  m_local = apf::getMesh(n); 
  el_type = m_local->getType(element);
  el_dim = m_local->typeDimension[el_type];
  num_comp = apf::countComponents(n);
  fshape_local = apf::getShape(n);
  e = entities[0];

//  std::cout << "fshape = " << fshape_local->getName() << std::endl;
  col = 0;
  ptr = 0;  // pointer to linear address in entities
  uint8_t offset_k = 0; // narrowing conversion error?
  int new_node;  // store new node values (calculated using offsets

  for (int i = 0; i <= el_dim; i++)   // loop over verts, edges, faces, regions
  {
//    std::cout << "looping over dimension " << i << std::endl;
    for (int j=0; j < m_local->adjacentCount[el_type][i]; j++)  // loop over all entities of this type
    {
      int type = apf::Mesh::simplexTypes[i];
//      std::cout << "  entity number " << j << std::endl;
//      std::cout << "  this entity has " << fshape_local->countNodesOn(type) << " nodes " << std::endl;

      for (int k=0; k < fshape_local->countNodesOn(type); k++)
      {
//        std::cout << "    node number " << k << std::endl;
        e = entities[col];  // get current entity
//        std::cout << "    col = " << col << std::endl;
//        std::cout << "    e = " << e << std::endl;
        offset_k = node_offsets[col];
//        std::cout << "    offset_k = " << offset_k + 0 << std::endl;
         
          // calculate new node index 
          // convert to 1 based indexing, do offsets, then convert back
          // to zero based indexing
          new_node = abs(offset_k - (k+1)) - 1;
//          std::cout << "      new_node = " << new_node << std::endl;
 
        for (int p = 0; p < num_comp; p++)
        {
          ptr = (nodemap[col] - 1)*num_comp + p;

//          std::cout << "      component number " << p << std::endl;
//          std::cout << "      ptr = " << ptr << std::endl;

          dofnums[ptr] = apf::getNumber(n, e, new_node, p);
        }
        col++;
      }  // end loop over nodes
    }  // end loops over entities of this type
  } // end loop over types


  return 0;
}





void setNumberingOffset(apf::Numbering* num, int off)
{
  apf::SetNumberingOffset(num, off);
}

// get node numbers in canonical ordering
// n is the Numbering
// e is the entity that contains the nodes
void getElementNumbers(apf::Numbering* n, apf::MeshEntity*e, int num_dof, int nums[])
{
  static apf::NewArray<int> nums_pumi;
  // get element numbers in canonical order
  apf::getElementNumbers(n, e, nums_pumi); 

  // copy to returned array
  for (int i = 0; i < num_dof; ++i)
  {
    nums[i] = nums_pumi[i];
  }
}


apf::Mesh* getMesh(apf::Numbering* n)
{
  return apf::getMesh(n);
}

void printNumberingName(apf::Numbering* n)
{
  std::cout << "name of numbering = " << apf::getName(n) << std::endl;
}

// initilizes tags on *mesh entities* (not the mesh field)
apf::MeshTag* createDoubleTag(apf::Mesh2 * m_local, char* name, int size)
{
  return m_local->createDoubleTag(name, size);
}

void setDoubleTag(apf::Mesh2 * m_local, apf::MeshEntity* e, apf::MeshTag* tag,  double data[])
{
//  std::cout << "in C++, about to set tag" << std::endl;
  m_local->setDoubleTag( e, tag, data);
//  std::cout << "finished seeting tag" << std::endl;
}

void getDoubleTag(apf::Mesh2 * m_local, apf::MeshEntity* e, apf::MeshTag* tag,  double* data)
{
  m_local->getDoubleTag( e, tag, data);
}





// mesh adapatation function
extern void createIsoFunc(apf::Mesh2* m_local, double(*sizefunc)(apf::MeshEntity*vert, apf::Mesh2* m_local, double* u), double *u)
{
  std::cout << "in createIsoFunc" << std::endl;
  IsotropicFunctionJ newisofunc(m_local, sizefunc, u); // create new function
  isofunc = newisofunc; // copy to global isofunc
}

// using a double* for the operator only works on 64 bit systems
void createAnisoFunc(apf::Mesh2* m_local,  void (*sizefunc)(apf::MeshEntity* vert, double r[3][3], double h[3], apf::Mesh2* m_ptr, void *f_ptr, double *operator_ptr), apf::Field *f_ptr, double *operator_ptr)
{
//  std::cout << " in c++, operator_ptr = " << operator_ptr << std::endl;

  AnisotropicFunctionJ newanisofunc(m_local, sizefunc, f_ptr, operator_ptr);  // create new function
  anisofunc = newanisofunc;  // copy to global anisofunc
}
// run mesh adaptation using isofunc
void runIsoAdapt(apf::Mesh2* m_local)
{
  IsotropicFunctionJ* isofunc_ptr = &isofunc;
  ma::Input* inputconfig = ma::configure(m_local, isofunc_ptr);
  ma::adapt(inputconfig);
}

void runAnisoAdapt(apf::Mesh2* m_local)
{
  AnisotropicFunctionJ* anisofunc_ptr = &anisofunc;
  ma::Input* inputconfig = ma::configure(m_local, anisofunc_ptr);
  ma::adapt(inputconfig);
}


// apf::Field functions (needed for automagical solution transfer)
apf::Field* createPackedField(apf::Mesh* m, char* fieldname, int numcomponents,  apf::FieldShape* fshape)
{
  return apf::createPackedField(m, fieldname, numcomponents, fshape);
}

// set the values of all components on a node
void setComponents(apf::Field* f, apf::MeshEntity* e, int node,  double const components[])
{
  apf::setComponents(f, e, node, components);
}

// copy the field values into the array components
// the array had better be the right size
void getComponents(apf::Field* f, apf::MeshEntity*e, int node, double components[])
{
  apf::getComponents(f, e, node, components);
}

void zeroField(apf::Field* f)
{
  apf::zeroField(f);
}

//-----------------------------------------------------------------------------
// Parallelization function
//-----------------------------------------------------------------------------

apf::Parts parts;
// count the number of peer parts
std::size_t countPeers(apf::Mesh* m, int dim)
{
  parts.clear();  // always empty the container before getting new elements
  apf::getPeers(m, dim, parts);
  return parts.size();
}

void getPeers(apf::Mesh*m, int part_nums[])
{
  int i = 0;

  for (apf::Parts::iterator it = parts.begin(); it != parts.end(); ++it)
  {
    part_nums[i] = *it;
    ++i;
  }

} // end function


int isShared(apf::Mesh* m, apf::MeshEntity* e)
{
  return m->isShared(e);
}

apf::Copies copies;
std::size_t countRemotes(apf::Mesh* m, apf::MeshEntity* e)
{
  copies.clear();  // always empty container before getting new elements
  m->getRemotes(e, copies);
  return copies.size();
}

void getRemotes(int part_nums[], apf::MeshEntity* entities[])
{
  int i = 0;
  for (apf::Copies::iterator it = copies.begin(); it != copies.end(); ++it)
  {
    part_nums[i] = it->first;
    entities[i] = it->second;
    ++i;
  }
}

//-----------------------------------------------------------------------------
// Mesh warping functions
//-----------------------------------------------------------------------------

// coords must be of length 3
void setPoint(apf::Mesh2* m, apf::MeshEntity* e, int node, double* coords)
{
  apf::Vector3 vec(coords);
  m->setPoint(e, node, vec);
}

void acceptChanges(apf::Mesh2* m)
{
  m->acceptChanges();
}

void Verify(apf::Mesh* m)
{
  m->verify();
}

void getPoint(apf::Mesh* m, apf::MeshEntity* e,  int node, double* coords)
{
  static apf::Vector3 vec;
  m->getPoint(e, node, vec);
  vec.toArray(coords);
}

//-----------------------------------------------------------------------------
// Matched entity functions
//-----------------------------------------------------------------------------

bool hasMatching(apf::Mesh* m)
{
  return m->hasMatching();
}

apf::Sharing* getSharing(apf::Mesh* m)
{
  return apf::getSharing(m);
}

bool isOwned(apf::Sharing* shr, apf::MeshEntity* e)
{
  return shr->isOwned(e);
}

apf::CopyArray copies1;
std::size_t countCopies(apf::Sharing* shr, apf::MeshEntity* e)
{
  // zero out the array
  // TODO: fix DynamicArray so reducing the sizes doesn't deallocate the memory
  //       and add a shrink function to shrink the memory to the current size
  copies1.setSize(0);
  shr->getCopies(e, copies1);
  return copies1.getSize();
}

void getCopies(int part_nums[], apf::MeshEntity* entities[])
{
  for (int i=0; i < int(copies1.getSize()); ++i)
  {
    part_nums[i] = copies1[i].peer;
    entities[i] = copies1[i].entity;
  }
}

apf::Matches matches1;
std::size_t countMatches(apf::Mesh* m, apf::MeshEntity* e)
{
  matches1.setSize(0);
  m->getMatches(e, matches1);
  return matches1.getSize();
}

void getMatches(int part_nums[], apf::MeshEntity* entities[])
{
  for (int i =0; i < int(matches1.getSize()); ++i)
  {
    part_nums[i] = matches1[i].peer;
    entities[i] = matches1[i].entity;
  }
}

