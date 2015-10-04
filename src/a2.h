#ifndef A2_H
#define A2_H

#include <iostream>
#include <deque>
#include <queue>


extern "C"
{
  extern bool hasNode(apf::Mesh2* m, apf::MeshEntity* e);
  extern void addQueues(std::queue<apf::MeshEntity*> & q1, std::queue<apf::MeshEntity*> & q2);
  extern apf::MeshEntity* getStartEntity(apf::Mesh2* & m);
  extern void printType(apf::Mesh* m, apf::MeshEntity* e);
  extern void numberdofs(apf::Mesh2* m, apf::Numbering* nodeNums, int numN);
  extern void numberElements(apf::Mesh2*& m, apf::Numbering*& elNums, int numEl);
  extern void printElNumbers(apf::Mesh2*& m, apf::Numbering*& elNums);
  extern void reorder(apf::Mesh2* m, int ndof, const int nnodes, const int comp, apf::Numbering* dof_statusNumbering, apf::Numbering* nodeNums, apf::Numbering* elNums, apf::MeshEntity* els_reordered[]);
}

#endif
