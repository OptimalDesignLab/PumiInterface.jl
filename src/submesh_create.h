// header file for submesh creation

#ifndef SUBMESH_CREATE_H
#define SUBMESH_CREATE_H
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <unordered_set>

// Pumi headers
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <ma.h>
//#include <stdlib.h>   // malloc, free, etc.
//#include <math.h>
//#include <string.h>

class SubMeshData
{
  public:
    SubMeshData() {};  // default constructor

    // create some data structures needed to construct the new mesh
    // this constructor also creates the new mesh object, but does not create
    // any mesh entities
    SubMeshData(apf::Mesh* _m_old, apf::Numbering* _numberings[], 
                int* _el_list, int numel);  // useful constructor


    // default destructor works here
    
    // data members
    apf::Mesh* m_old = NULL;
    std::vector<apf::Numbering*> numberings;
  //  apf::Numbering* numberings[] = {NULL};
    apf::Mesh2* m_new = NULL;
    std::vector<int> el_list; // 0-based
    int dim = 0;  // dimension of mesh

    // map of number on the old mesh to the MeshEntity* on the new mesh
    // length = number of entities on old mesh, unused values set to NULL
    std::vector<apf::MeshEntity*> verts;
    std::vector<apf::MeshEntity*> edges;
    std::vector<apf::MeshEntity*> faces;
    std::vector<apf::MeshEntity*> regions;

    // array of the above for easy looping
    std::vector<apf::MeshEntity*> entities[4];

    // vector of element MeshEntity* on old mesh
    std::vector<apf::MeshEntity*> el_entities;
    //  std::unordered_set<apf::MeshEntity*> el_set; // el_list with O(1) test for containment
   
    // member functions
    void writeNewMesh(const char* fname)
    {
      m_new->writeNative(fname);
    }

    apf::Mesh2* getNewMesh()
    {
      return m_new;
    }
};  // class SubMeshData

// vertify input data is valid
void checkInput(SubMeshData* sdata);
void createVertices(SubMeshData* sdata);
void createEntities(SubMeshData* sdata);
// create a vertex on the new mesh given a vertex on the old mesh
// uses same geometry classification as original mesh
apf::MeshEntity* createVert(SubMeshData* sdata, apf::MeshEntity* vert);

// create an entity from its one-level downward adjacencies
apf::MeshEntity* createEntity(SubMeshData* sdata, apf::MeshEntity* entity);



extern "C" {
  apf::Mesh2* createSubMesh(apf::Mesh* m, int* el_list, int numel);
  void writeNewMesh(SubMeshData* sdata, const char* fname);
  apf::Mesh2* getNewMesh(SubMeshData* sdata);
} // extern C

#endif
