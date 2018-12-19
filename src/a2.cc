#include <climits>
#include <cassert>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <deque>
#include <queue>

#include "a2.h"

bool hasNode(apf::Mesh2* m_local, apf::MeshEntity* e)
{
  return m_local->getShape()->countNodesOn(m_local->getType(e)) > 0;
}

int nodeCount(apf::Mesh2* m_local, apf::MeshEntity* e)
{
  return m_local->getShape()->countNodesOn(m_local->getType(e));
}

// determine if a dof should be numbered or not
// Takes into account whether node_statusNumbering is NULL or not
bool shouldNumber(apf::Numbering* node_statusNumbering, apf::MeshEntity* e,
                  int node, int component)
{
  bool ret = true;
  // number if dof is free or loaded
  if (node_statusNumbering)
    ret = (apf::getNumber(node_statusNumbering, e, node, component) >= 2);

  return ret;

}


// returns the number of geometric dofs on each node of the given MeshEntity
int getNumXiDof(apf::Mesh* m, apf::MeshEntity* e)
{
  apf::ModelEntity* me = m->toModel(e);
  return m->getModelType(me);
}

// adds p2 to the end of q1
void addQueues(std::queue<apf::MeshEntity*> & q1, std::queue<apf::MeshEntity*> & q2)
{
  int sizeq2 = q2.size();
  for ( int i = 0; i < sizeq2; ++i)
  {
    q1.push(q2.front());
    q2.pop();
  }


}

// get starting entity for node reordering
// search for node vertex classified on geometric vertex that is closest to 
// the point (x,y) and has minimum connectivity
//
apf::MeshEntity* getStartEntity(apf::Mesh2* & m_local, const double start_coords[3])
{
  double x = start_coords[0];
  double y = start_coords[1];
  double z = start_coords[2];


  if (PCU_Comm_Self() == 0)
    std::cout << "requested starting coordinates = " << x << ", " << y << ", " << z << std::endl;
  apf::MeshEntity* e_min; // minimum degree meshentity
  apf::MeshEntity* e_i; // current meshentity
  apf::Vector3 coords; // coordinates of current point
  double dist;  // distance from coords to (x,y)
  double min_dist; // minimum distance to (x,y)
  apf::ModelEntity* me_i;
  int me_dimension;
  x = start_coords[0];
  y = start_coords[1];
  z = start_coords[2];

  apf::MeshIterator* it = m_local->begin(0); // iterator over verticies

  gmi_model* g = m_local->getModel();
/*
  std::cout << "g = " << g << std::endl;
  std::cout << "g.n[0] = " << g->n[0] << std::endl;
  std::cout << "g.n[1] = " << g->n[1] << std::endl;
  std::cout << "g.n[2] = " << g->n[2] << std::endl;
  std::cout << "g.n[3] = " << g->n[3] << std::endl;
*/
  bool is_null = g->n[0] == 0;
  // initilize
  e_i = m_local->deref(it);

  e_min = e_i;  // ensure we return a value if conditions are never
                // satisfied

  // calculate distance to (x,y)
  m_local-> getPoint(e_i, 0, coords);
  // no need to take square root if we are only interested in relative
  // distance
  min_dist = (coords[0] - x)*(coords[0] - x) + (coords[1] - y)*(coords[1] - y) + (coords[2] - z)*(coords[2] - z);



  while ( (e_i = m_local->iterate(it)) )
  {
    me_i = m_local->toModel(e_i);
    me_dimension = m_local->getModelType(me_i);
    if ( !me_dimension || is_null ) // if me_dimension == 0
    {

      // calculate distance to (x,y)
      m_local-> getPoint(e_i, 0, coords);
      // no need to take square root if we are only interested in relative
      // distance
      dist = (coords[0] - x)*(coords[0] - x) + (coords[1] - y)*(coords[1] - y) + (coords[2] - z)*(coords[2] - z);

      if (dist < min_dist)
      {
        e_min = e_i;
        min_dist = dist;
//        std::cout << "choosing this vertex" << std::endl;
      }
    }
  }

  m_local->getPoint(e_min, 0, coords);

  if (PCU_Comm_Self() == 0)
    std::cout << "starting entity coordinates = " << coords[0] << ", " << coords[1] << ", " << coords[2]  << std::endl;

  return e_min;

}


void printType(apf::Mesh* m_local, apf::MeshEntity* e)
{
  int type_enum = m_local->getType(e);
  // check type, print appropriate message
    
  switch(type_enum)
    {
        case apf::Mesh::TET :
          std::cout << " tetrahedron" << std::endl;
          break;
        case apf::Mesh::HEX :
          std::cout << " Hexahedron" << std::endl;
           break;
        case apf::Mesh::PRISM :
           std::cout << " prism" << std::endl;
           break;
        case apf::Mesh::PYRAMID :
          std::cout << " pyramid" << std::endl;
          break;
        case apf::Mesh::QUAD :
          std::cout << " quad" << std::endl;
          break;
        case apf::Mesh::TRIANGLE :
          std::cout << " triangle" << std::endl;
          break;
        case apf::Mesh::VERTEX :
          std::cout << "vertex" << std::endl;
          break;
        case apf::Mesh::EDGE :
          std::cout << "edge" << std::endl;
          break;
        default:
             std::cout << " has no matching type" << std::endl;
    }
}



// initially number all dofs with number greater than number of nodes, to show
// they have not received final number yet

void numberdofs(apf::Mesh2* m_local, apf::Numbering* nodeNums, int ndof, int comp)
{
//  apf::FieldShape* fieldshape = m_local->getShape();
  apf::MeshIterator* it;
  apf::MeshEntity* e;
  int numNodes_typei;
//  int numNodes_j;
  int k = ndof + 1;

  for (int i = 0; i < 4; ++i) // loop over entity types
  {
    it = m_local->begin(i);
    e = m_local->deref(it);
    numNodes_typei = nodeCount(m_local,e);
//    std::cout << "entity type " << i << " has " << numNodes_typei << " nodes" << std::endl;
    it = m_local->begin(i);

    if (numNodes_typei)  // if there are any dofs on this type of entity
    {
      
      while ( (e = m_local->iterate(it)) )
      {
        for ( int j = 0; j < numNodes_typei; ++j)
        {
          for ( int c = 0; c < comp; ++c)
          {
            apf::number(nodeNums, e, j, c, k); // label current node
//            std::cout << "initially numbered dof " << k << std::endl;
            k += 1;
          }
        }
      }
    }


  }

}

// counts the number of geometric degrees of freedom (entities on the boundary
// use CAD parametric coordinates)
int countXiDofs(apf::Mesh2* m)
{
  int ndof = 0;
  apf::FieldShape* fshape = m->getShape();
  apf::MeshEntity* e;
  apf::MeshIterator* it;

  for (int dim=0; dim <= m->getDimension(); ++dim)
  {
    if (!fshape->hasNodesIn(dim))
      continue;

    it = m->begin(dim);
    while ( (e = m->iterate(it))) 
    {
      apf::ModelEntity* me = m->toModel(e);
      int me_dim = m->getModelType(me);
      // the number of geometric dofs is always equal to the geometric dimension
      ndof += me_dim*fshape->countNodesOn(m->getType(e));
    }

    m->end(it);
  }  // end dim

  return ndof;
}  // end function

// do initial numbering of all dofs as ndof + 2
// returns ndof
int numberXiDofs(apf::Mesh2* m, apf::Numbering* xiNums)
{
  int ndof = countXiDofs(m);

  apf::FieldShape* fshape = m->getShape();
  apf::MeshEntity* e;
  apf::MeshIterator* it;

  for (int dim=0; dim <= m->getDimension(); ++dim)
  {
    if (!fshape->hasNodesIn(dim))
      continue;

    it = m->begin(dim);
    while ( (e = m->iterate(it)) )
    {
      for (int i=0; i < fshape->countNodesOn(m->getType(e)); ++i)
        for (int d1=0; d1 <= m->getDimension(); ++d1)
          apf::number(xiNums, e, i, d1, ndof + 2);

    }  // end while

    m->end(it);

  }  // end for dim

  return ndof;
}



void numberEntity(NumberXi& nxi, apf::MeshEntity* e, int& nodeLabel_i)
{
  int numNodes_i = nodeCount(nxi.m, e);
  int ncomp = getNumXiDof(nxi.m, e);
  for ( int i = 0; i < numNodes_i; ++i)
  {
    for ( int c = 0; c < ncomp; ++c) // loop over dof of the node
      apf::number(nxi.xiNums, e, i, c, --nodeLabel_i);

    for ( int c = ncomp; c < nxi.m->getDimension(); ++c)
      apf::number(nxi.xiNums, e, i, c, nxi.NODE_UNUSED);
  }
}  // function




// number elements with numbers greater than number of elements to show they 
// have not yet received their final number
void numberElements(apf::Mesh2* m_local, apf::Numbering* elNums, int numEl)
{

  if (elNums)
  {
    apf::MeshIterator* it = m_local->begin(2);
    apf::MeshEntity* e;
    int k = numEl + 1;

    while ( (e = m_local->iterate(it) ) )
    {
  //    std::cout << "labelling element " << k - numEl << " as " << k << std::endl;
      apf::number(elNums, e, 0, 0, k);
      ++k;
    }
  }

}

void printElNumbers(apf::Mesh2*& m_local, apf::Numbering*& elNums)
{
  apf::MeshIterator* it = m_local->begin(2);
  apf::MeshEntity* e;
  int i = 1;
  while ((e = m_local->iterate(it) ) )
  {
    int num = apf::getNumber(elNums, e, 0 , 0);
    std::cout << "element " << i << "1 number = " << num << std::endl;
    ++i;
  }
}

// reorder mesh nodes and elements
// multiple dof on the same node are labelled sequentially
// this is not actually a good idea because it will double the bandwidth of hte 
// sitffness matrix, but modifying the algorithm to treat them as separate
// dofs would be difficult
// return numberings of each, and array of element pointers 
// (so they can be iterated over in order)
// node_statusNumbering tells whether a dof is fixed (=3) or not.  If NULL,
// all dofs are numbered
// ndof is the number of actual degrees of freedom in the mesh
// comp is the number of dofs per node (the number of components in the dof numberings)
// nodeNums is a numbering over the dofs to be populated with  node numbers
// elNums is numbering over elements (faces) to be populated, if NULL, elements
// will not be numbered
//
// x, y, z coordinates in start_coords are used to deterimine the starting node
//   the mesh vertex classified on a model vertex closest to the point is chosen
void reorder(apf::Mesh2* m_local, int ndof, const int comp, 
             apf::Numbering* node_statusNumbering, apf::Numbering* nodeNums,
             apf::Numbering* elNums, const double start_coords[3])
{
// TODO: move node_statusNumbering checks out one loop level because 
//       it is node status now, not dof status

  // check that we won't overflow integers
  if (ndof > INT_MAX/4)
  {
    std::cerr << "Error: Cannot number this many dofs using C integers" << std::endl;
    return;
  }

  const int dim = m_local->getDimension();
  const int numEl = m_local->count(m_local->getDimension());  // counts the number of elements


  // initially number in range numnodes + 1 to 2*numnodes
  numberElements(m_local, elNums, numEl);
  numberdofs(m_local, nodeNums, ndof, comp);


  // create queues
  std::queue < apf::MeshEntity*> que1;  // main queue
  std::queue < apf::MeshEntity*> tmpQue;  // temporary queue for vertices

  apf::MeshEntity* e;
  e = getStartEntity(m_local, start_coords); // get starting node

  int nodeLabel_i = ndof + 1;  // one-based indexing
  int elementLabel_i = numEl;  // zero-based indexing
  int numNodes_i;

  // queue initial entity
  que1.push(e);
  
  while (que1.size() > 0 )
  {
    e = que1.front(); // get next element in queue
    que1.pop();  // remove that element from the que

    // label all nodes on entity e
    numNodes_i = nodeCount(m_local,e);
    for ( int i = 0; i < numNodes_i; ++i)
    {
      for ( int c = 0; c < comp; ++c) // loop over dof of the node
      {
        if (shouldNumber(node_statusNumbering, e, i, c))
          number(nodeNums, e, i, c, --nodeLabel_i);
        else // give all fixed dofs a label of 0 so they get written to vtk
           apf::number(nodeNums, e, i, c, 0);
      }
    }  // end for i

    // if e is a vertex, find adjacencies, look for unlabeled nodes
    if (m_local->getType(e) == apf::Mesh::VERTEX)
    {
      int numEdges_w = m_local->countUpward(e);
      for (int i = 0; i < numEdges_w; ++i)  // loop over edges
      {
        apf::MeshEntity* edge_i = m_local->getUpward(e, i);
        int numFaces_i = m_local->countUpward(edge_i);

        // queue faces adjacent to the edge if needed
        for (int j = 0; j < numFaces_i; ++j)  // loop over faces
        {
          apf::MeshEntity* face_j = m_local->getUpward(edge_i, j);

          if (elNums && dim == 2)
          {
            // label face (element) if it is not yet labelled
            int faceNum_j = apf::getNumber(elNums, face_j, 0, 0);
            if (faceNum_j > numEl) // if face not labelled
            {
              elementLabel_i -= 1;
              apf::number(elNums, face_j, 0, 0, elementLabel_i);
            }
          }

          // add face to queue if it hasn't been labelled yet
          if ( hasNode(m_local, face_j) )
          {
            int  nodeNum_j = apf::getNumber(nodeNums, face_j, 0, 0); 
            // node not labelled) and not in queue
            if ( (nodeNum_j > ndof) && (nodeNum_j <= 2*ndof)) 
            {
              // double node number to show is has been added to queue
              apf::number(nodeNums, face_j, 0, 0, nodeNum_j*2);
              que1.push(face_j); // add face to que
            }
          }

          // get regions adjacent to the face
          if (dim == 3)
          {
            int numRegions_i = m_local->countUpward(face_j);
            for (int k = 0; k < numRegions_i; ++k)
            {
              apf::MeshEntity* region_k = m_local->getUpward(face_j, k);

              // label element
              if (elNums)
              {
                int regionNum_k = apf::getNumber(elNums, region_k, 0, 0);
                if (regionNum_k > numEl)  // face not labelled
                {
                  elementLabel_i -= 1;
                  apf::number(elNums, region_k, 0, 0, elementLabel_i);
                }
              }

              // add region to queue if it hasn't been labelled yet
              if ( hasNode(m_local, region_k) )
              {
                int nodenum_k = apf::getNumber(nodeNums, region_k, 0, 0);
                if ((nodenum_k > ndof) && (nodenum_k <= 2*ndof))
                {
                  apf::number(nodeNums, region_k, 0, 0, nodenum_k*2);
                  que1.push(region_k);
                }
              }
            }  // end region loop
          }
        }  // end face loop

        // look at other vertex on edge
        apf::MeshEntity* otherVertex = apf::getEdgeVertOppositeVert(m_local, edge_i, e);
        int otherVertex_num = getNumber(nodeNums, otherVertex, 0,0);
        bool labeled = (otherVertex_num <= ndof);
        bool queued = (otherVertex_num > 2*ndof);  // vertex already queued

        if (hasNode(m_local, edge_i))
        {
          int edgeNode_num = getNumber(nodeNums, edge_i,0,0);
          bool edgeNotLabeled = ( (edgeNode_num > ndof) && ( edgeNode_num <= 2*ndof) ); // edge not labeled nor in que


          if ((labeled || queued) && edgeNotLabeled)
          {
            // label all nodes on edge
            int numEdgeNodes = nodeCount(m_local, edge_i);
            for (int j = 0; j < numEdgeNodes; ++j)
            {
              for (int c = 0; c < comp; ++c) // loop over dofs of node
              {
               if (shouldNumber(node_statusNumbering, edge_i, j, c))
                  apf::number(nodeNums, edge_i, j, c, --nodeLabel_i);
               else
                  apf::number(nodeNums, edge_i, j, c, 0);
              }
            }
          } else { // add edge to queue first, othervertex second (via list)   
          
              if (edgeNotLabeled) 
              {
                apf::number(nodeNums, edge_i, 0, 0, edgeNode_num*2);  
                que1.push(edge_i);
              }

              if ( (!labeled) && (!queued) )
              {
                // double node number to signal this node is queued
                apf::number(nodeNums, otherVertex, 0, 0, 2*otherVertex_num);
                tmpQue.push(otherVertex);
              }
          }          

        } else  // if edge does not have node, deal with otherVertex
        {
          if ( (!labeled) && (!queued))  // if otherVertex not labelled
          {
            // double node number to signal this node queued
            apf::number(nodeNums,otherVertex, 0,0, 2*otherVertex_num);
            tmpQue.push(otherVertex);
          }
        }

      }  // end edge loop

      // copy tmpQue into que1, empty tmpQue
      addQueues(que1, tmpQue);


      }  // end if (vertex)
    } // end while loop over que    


  if (elNums)
  {
    if (elementLabel_i != 0)
    {
      std::cerr << "Warning: element numbering not sane" << std::endl;
      std::cerr << "final elementLabel_i = " << elementLabel_i << std::endl;
    } else
    {
      std::cout << "element reordering is sane" << std::endl;
    }
  }

  if (nodeLabel_i != 1)
  {
    std::cerr << "Warning: node numbering not sane" << std::endl;
    std::cerr << "final nodeLabel_i = " << nodeLabel_i << std::endl;
  } else
  {
//    std::cout << "node reordering is sane" << std::endl;
  }

  assert(nodeLabel_i == 1);
}


// similar to reorder(), but numbers geometric dofs
// xiNums must have the same FieldShape as the mesh
int reorderXi(apf::Mesh2* m_local, apf::Numbering* xiNums,
               const double start_coords[3])
{

  const int dim = m_local->getDimension();

  // initially number in range numnodes + 1 to 2*numnodes
  int ndof = numberXiDofs(m_local, xiNums);
  const int NODE_UNUSED = ndof + 1;
  const int NODE_UNSEEN = ndof + 2;  // not numbered, not in queue
  const int NODE_QUEUED = ndof + 3;  // in queue (not numbered)
  NumberXi nxi(m_local, xiNums, NODE_UNUSED, NODE_UNSEEN, NODE_QUEUED);

  // create queues
  std::queue < apf::MeshEntity*> que1;  // main queue
  std::queue < apf::MeshEntity*> tmpQue;  // temporary queue for vertices

  apf::MeshEntity* e;
  e = getStartEntity(m_local, start_coords); // get starting node

  int nodeLabel_i = ndof + 1;  // one-based indexing

  // queue initial entity
  que1.push(e);
  
  while (que1.size() > 0 )
  {
    e = que1.front(); // get next element in queue
    que1.pop();  // remove that element from the que

    // label all nodes on entity e
    numberEntity(nxi, e, nodeLabel_i);


    // if e is a vertex, find adjacencies, look for unlabeled nodes
    if (m_local->getType(e) == apf::Mesh::VERTEX)
    {
      int numEdges_w = m_local->countUpward(e);
      for (int i = 0; i < numEdges_w; ++i)  // loop over edges
      {
        apf::MeshEntity* edge_i = m_local->getUpward(e, i);
        int numFaces_i = m_local->countUpward(edge_i);

        // queue faces adjacent to the edge if needed
        for (int j = 0; j < numFaces_i; ++j)  // loop over faces
        {
          apf::MeshEntity* face_j = m_local->getUpward(edge_i, j);

          // add face to queue if it hasn't been labelled yet
          if ( hasNode(m_local, face_j) )
          {
            int nodeNum_j = apf::getNumber(xiNums, face_j, 0, 0); 
            // node not labelled) and not in queue
            if (nodeNum_j == NODE_UNSEEN) 
            {
              // double node number to show is has been added to queue
              apf::number(xiNums, face_j, 0, 0, NODE_QUEUED);
              que1.push(face_j); // add face to que
            }
          }

          // get regions adjacent to the face
          if (dim == 3)
          {
            int numRegions_i = m_local->countUpward(face_j);
            for (int k = 0; k < numRegions_i; ++k)
            {
              apf::MeshEntity* region_k = m_local->getUpward(face_j, k);

              // add region to queue if it hasn't been labelled yet
              if ( hasNode(m_local, region_k) )
              {
                int nodenum_k = apf::getNumber(xiNums, region_k, 0, 0);
                if (nodenum_k == NODE_UNSEEN)
                {
                  apf::number(xiNums, region_k, 0, 0, NODE_QUEUED);
                  que1.push(region_k);
                }
              }
            }  // end region loop
          }
        }  // end face loop

        // look at other vertex on edge
        apf::MeshEntity* otherVertex = apf::getEdgeVertOppositeVert(m_local, edge_i, e);
        int otherVertex_num = getNumber(xiNums, otherVertex, 0,0);
        bool labeled = (otherVertex_num <= NODE_UNUSED);
        bool queued = (otherVertex_num == NODE_QUEUED);  // vertex already queued

        if (hasNode(m_local, edge_i))
        {
          int edgeNode_num = getNumber(xiNums, edge_i,0,0);
          bool edgeNotLabeled = ( edgeNode_num == NODE_UNSEEN ); // edge not labeled nor in que


          if ((labeled || queued) && edgeNotLabeled)
          {
            numberEntity(nxi, edge_i, nodeLabel_i);
          } else { // add edge to queue first, othervertex second (via list)   
          
              if (edgeNotLabeled)
              {
                apf::number(xiNums, edge_i, 0, 0, NODE_QUEUED);
                que1.push(edge_i);
              }

              if ( (!labeled) && (!queued) )
              {
                // double node number to signal this node is queued
                apf::number(xiNums, otherVertex, 0, 0, NODE_QUEUED);
                tmpQue.push(otherVertex);
              }
          }          

        } else  // if edge does not have node, deal with otherVertex
        {
          if ( (!labeled) && (!queued))  // if otherVertex not labelled
          {
            // double node number to signal this node queued
            apf::number(xiNums, otherVertex, 0, 0, NODE_QUEUED);
            tmpQue.push(otherVertex);
          }
        }

      }  // end edge loop

      // copy tmpQue into que1, empty tmpQue
      addQueues(que1, tmpQue);

      }  // end if (vertex)
    } // end while loop over que    


  if (nodeLabel_i != 1)
  {
    std::cerr << "Warning: node numbering not sane" << std::endl;
    std::cerr << "final nodeLabel_i = " << nodeLabel_i << std::endl;
  } else
  {
//    std::cout << "node reordering is sane" << std::endl;
  }

  assert(nodeLabel_i == 1);

  return ndof;
}

