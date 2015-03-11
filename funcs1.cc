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


#include <iostream>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <stdlib.h>   // malloc, free, etc.
#include <math.h>

#include "funcs1.h"

//=============================================================================
//declare global variables (persistent state of library)

apf::Mesh2* m;
apf::FieldShape* mshape;

apf::Numbering* elNums; // element numbering
apf::Numbering* faceNums; // face numbering
apf::Numbering* edgeNums; // edge numbering
apf::Numbering* vertNums; // vertex numbers


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



// initilize global variables, used by all fucntions
// downward_counts = numDown
// number_entities = numEntity
// m_ptr = mesh pointer
// mshape_ptr = pointer to mesh shape (m->getShape())
int initABC(char* dmg_name, char* smb_name, int downward_counts[4][4], int number_entities[4], apf::Mesh2* m_ptr_array[1], apf::FieldShape* mshape_ptr_array[1] )
{
  std::cout << "Entered init\n" << std::endl;

  MPI_Init(0,NULL);  // initilize MPI 
  PCU_Comm_Init();   // initilize PUMI's communication
  gmi_register_mesh();

  // load mesh
//  m = apf::loadMdsMesh("cube.dmg", "tet-mesh-1.smb");
  m = apf::loadMdsMesh(dmg_name, smb_name);
  m_ptr_array[0] = m;
  mshape_ptr_array[0] = m->getShape();
  std::cout << std::endl;
/* 
  // initilize iterators
  elIt = m->begin(3);
  its[2] = m->begin(2);
  its[1] = m->begin(1);
  its[0] = m->begin(0);
  
*/
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

  // initalize tags
  globalVertNums = m->createIntTag("globalNodeNumber", 1);


  // initilize iterators
  its[0] = m->begin(0);
  its[1] = m->begin(1);
  its[2] = m->begin(2);
  its[3] = m->begin(3);
  for (int i = 0; i < 4; ++i)
  {
    int type = m->getType(m->deref(its[i]));
    std::cout << "type of its[" << i << "] = " << type;
    std::cout << "index its[" << i << "] = " << type << std::endl;
  }

    
  // initilize number of downward adjacencies each type has, assuming all
  // elements have same number of downward adjacencies
  apf::MeshEntity* e_tmp;
  apf::Downward tmp;
  for (int i = 0; i < 4 ; ++i)
  {
    e_tmp = m->deref(its[i]);
    for (int j = 0; j < 4; ++j)
    {
      if ( j < i)
      {
        numDown[i][j] = m->getDownward(e_tmp, j, tmp);
        downward_counts[j][i] = numDown[i][j];
      }
      else
      {
        numDown[i][j] = 0;
        downward_counts[j][i] = 0;
      }
      std::cout << names[i] << " has " << numDown[i][j];
      std::cout << " downward adjacencies of type " << names[j] << std::endl;
    }
  }

  std::cout << std::endl;
/*
  numberings[0] = numberOwnedDimension(m, "elNums", 0);
  numberings[1] = numberOwnedDimension(m, "elNums", 1);
  numberings[2] = numberOwnedDimension(m, "elNums", 2);
  numberings[3] = numberOwnedDimension(m, "elNums", 3);
*/

  std::cout << "numV = " << numEntity[0] << " , numEdge = " << numEntity[1];
  std::cout << " , numFace = " << numEntity[2] << " , numEl = " << numEntity[3] << std::endl;
  std::cout << std::endl;

  return 0;
}

apf::Mesh2* getMeshPtr()
{
  return m;
}

apf::FieldShape* getMeshShapePtr()
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


void incrementVertIt()
{
  m->iterate(its[0]);
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


// tag current vertex with val, part of MeshTag globalNodeNumbers
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

// get dimension of mesh
int getMeshDimension(apf::Mesh2* m_local)
{
  int i = m_local->getDimension();
  std::cout << "Mesh Dimesion = " << i << std::endl;
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
  std::cout << "in C++ getDownwards, dimension = " << dimension << std::endl;
  int numDown = m->getDownward(e, dimension, downwards);
  return numDown;
}



// check that global variable are persisent
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
  std::cout << "in getVertCoords" << std::endl;
  std::cout << "sx = " << sx << " ,sy = " << sy << std::endl;
  apf::MeshEntity* e = m->deref(its[0]);
  apf::Vector3 vec;
  m->getPoint(e, 0, vec);
  std::cout << "coords = " << vec << std::endl;
  coords[0][0] = vec[0];
  coords[0][1] = vec[1];
  coords[0][2] = vec[2];

}

// get the coordinates of the two points that define, iterate edge iterator
// an edge
int getEdgeCoords(double coords[2][3], int sx, int sy)
{
  std::cout <<"in getEdgeCoords" << std::endl;
  if (sx < 2 || sy != 3)
  {
    std::cout << "Warning: coords array too small" << std::endl;
    return -1;
  }

  apf::MeshEntity* e = m->deref(its[1]);

  apf::Vector3 vec;
  apf::MeshEntity* verts[2];
  std::cout << " about to get downward entities" << std::endl;
  m->getDownward(e, 0, verts);
  std::cout << "got vertices of edge" << std::endl;
  // get coordinates of each vertex
  m->getPoint(verts[0], 0, vec);
  coords[0][0] = vec[0];
  coords[0][1] = vec[1];
  coords[0][2] = vec[2];

  std::cout << "first point coordinates = " << vec << std::endl;

  m->getPoint(verts[1], 0, vec);
  coords[1][0] = vec[0];
  coords[1][1] = vec[1];
  coords[1][2] = vec[2];

  std::cout << "second point coordinates = " << vec << std::endl;

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
    std::cout << "coordinates of vertex " << i << " = " << vec << std::endl;
    coords[i][0] = vec[0];
    coords[i][1] = vec[1];
    coords[i][2] = vec[2];
  }

  return 0;
}

int getElCoords(double coords[][3], int sx, int sy)
{
  std::cout << "Entered getElCoords" << std::endl;

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
    std::cout << "coordinates of vertex " << i << " = " << vec << std::endl;
    coords[i][0] = vec[0];
    coords[i][1] = vec[1];
    coords[i][2] = vec[2];
  }

  return 0;
}

/*
int countNodesOn(apf::FieldShape* field, int type)
{
//  apf::FieldShape* field2 = m->getShape();
//  int i = field2->countNodesOn(2);
  int i = field->countNodesOn(type);
  return i;
}
*/

// create a generally defined numbering from julia
apf::Numbering* createNumberingJ(apf::Mesh2* m_local, char* name, apf::FieldShape* field, int components)
{
//    return apf::createNumbering(m, "num1", m->getShape(), 1);
//    return apf::createNumbering(m_local, "num1", m->getShape(), 1);
//    return apf::createNumbering(m_local, "num1", m_local->getShape(), components);
    return apf::createNumbering(m_local, name, field, components);
}

// number an entity in a given numbering from julia
int numberJ(apf::Numbering* n, apf::MeshEntity* e, int node, int component, int number)
{
  apf::number(n, e, node, component, number);
  int i = apf::getNumber(n,e,node,component);
  return i;
}

// retrieve a number from julia
int getNumberJ(apf::Numbering* n, apf::MeshEntity* e, int node, int component)
{
  int i = apf::getNumber(n, e, node, component);
  return i;
}


int countNodesOn(apf::FieldShape* mshape_ptr, int type)
{
  return mshape_ptr->countNodesOn(type);
}

extern void printNumberingName(apf::Numbering* n)
{
  std::cout << "name of numbering = " << apf::getName(n) << std::endl;
}

