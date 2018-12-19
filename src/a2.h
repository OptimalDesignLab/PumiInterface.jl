#ifndef A2_H
#define A2_H

#include <iostream>
#include <deque>
#include <queue>

// for passing data to numberEntity()
struct NumberXi
{
  NumberXi() = default;
  NumberXi(apf::Mesh2* _m, apf::Numbering* _n, int node_unused,
      int node_unseen, int node_queued) : m(_m), xiNums(_n),
      NODE_UNUSED(node_unused), NODE_UNSEEN(node_unseen),
      NODE_QUEUED(node_queued) {}
  apf::Mesh2* m;
  apf::Numbering* xiNums;
  const int NODE_UNUSED;
  const int NODE_UNSEEN;
  const int NODE_QUEUED;
};


extern "C"
{
  extern bool hasNode(apf::Mesh2* m, apf::MeshEntity* e);
  extern void addQueues(std::queue<apf::MeshEntity*> & q1, std::queue<apf::MeshEntity*> & q2);
  extern apf::MeshEntity* getStartEntity(apf::Mesh2* & m, const double x, const double y);
  extern void printType(apf::Mesh* m, apf::MeshEntity* e);
  extern void numberdofs(apf::Mesh2* m, apf::Numbering* nodeNums, int numN);

  int countXiDofs(apf::Mesh2* m);
  int numberXiDofs(apf::Mesh2* m, apf::Numbering* xiNums);

  void numberEntity(NumberXi& nxi, apf::MeshEntity* e, int& nodeLabel_i);
  extern void numberElements(apf::Mesh2* m, apf::Numbering* elNums, int numEl);
  extern void printElNumbers(apf::Mesh2*& m, apf::Numbering*& elNums);

  void reorder(apf::Mesh2* m_local, int ndof, const int comp, apf::Numbering* node_statusNumbering, apf::Numbering* nodeNums, apf::Numbering* elNums, const double start_coords[3]);

  int reorderXi(apf::Mesh2* m_local, apf::Numbering* xiNums,
               const double start_coords[3]);

}

#endif
