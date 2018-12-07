#ifndef A2_H
#define A2_H

#include <iostream>
#include <deque>
#include <queue>


extern "C"
{
  extern bool hasNode(apf::Mesh2* m, apf::MeshEntity* e);
  extern void addQueues(std::queue<apf::MeshEntity*> & q1, std::queue<apf::MeshEntity*> & q2);
  extern apf::MeshEntity* getStartEntity(apf::Mesh2* & m, const double x, const double y);
  extern void printType(apf::Mesh* m, apf::MeshEntity* e);
  extern void numberdofs(apf::Mesh2* m, apf::Numbering* nodeNums, int numN);
  extern void numberElements(apf::Mesh2* m, apf::Numbering* elNums, int numEl);
  extern void printElNumbers(apf::Mesh2*& m, apf::Numbering*& elNums);

void reorder(apf::Mesh2* m_local, int ndof, const int comp, apf::Numbering* node_statusNumbering, apf::Numbering* nodeNums, apf::Numbering* elNums, const double start_coords[3]);

}

#endif
