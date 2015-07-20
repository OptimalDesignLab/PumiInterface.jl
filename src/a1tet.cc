#include <iostream>  // includes cin and cout

#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include "funcs1.h"  // include user defind functions

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh("cube.dmg", "tet-mesh-1.smb");
  
  init(m);  // initilize the state of the function file


  // Problem 1
  std::cout << "Problem 1 output: " << std::endl;
  apf::MeshIterator* it = m->begin(3); // iterator over elements of mesh
  apf::MeshEntity* e;   // used during iteration to hold pointer to the current element
  double volume_e;  // holds volume of element e
  int i = 0;  // element counter
  while (e = m->iterate(it))   // loop over elements
  {
    volume_e = calcElementVolume_tet(e);  // get volume of currnet element
    std::cout << "Element " << i << " , volume = " << volume_e << std::endl;
    ++i;
  }

  std::cout << "\n" << std::endl;

  // Problem 2
  std::cout << "Problem 2 output: " << std::endl;
  // use MeshEntity* variable e from problem 1
  it = m->begin(0);   // iterator over vetices
  int numBound;   // holds number of edges bounded by a vertex
  i = 0;  // reset counter i to 0
  while (e = m->iterate(it))  // loop over vertices
  {
    numBound = getNumEdgesfromVertex(e);  // get number of edges
    std::cout << "Vertex " << i << " bounds " << numBound << " edges" << std::endl;
    ++i;
  }

  std::cout << "\n" << std::endl;

  //Problem 6
  std::cout << "Problem 6 output: " << std::endl;
  it = m->begin(1);  // iterator over edges
  apf::Vector3 r1;  // hold coordinates of first vertex of edge
  apf::Vector3 r2;  // hold coordinates of second vertex of edge
  apf::Vector3 r3;  // holds coordinates of new vertex
  i = 0;  // reset counter to zero

  while (e = m->iterate(it))
  {
    getRsfromEdge(e, r1, r2);  // put coordinates of vertices in r1, r2
    r3 = r1 + r2;
    r3 = r3/2;   // average coordinate
    createPoint(e, r3);  // user defind function to create a point 

    std::cout << "adding point to edge " << i << " at location " << r3 << std::endl;
    ++i;
  }


  std::cout << "\n" << std::endl;
  std::cout << "Program finished" << std::endl;

  apf::writeVtkFiles("outTet", m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
