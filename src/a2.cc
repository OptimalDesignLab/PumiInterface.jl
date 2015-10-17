#include <climits>
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
apf::MeshEntity* getStartEntity(apf::Mesh2* & m_local, const double x, const double y)
{
  apf::MeshEntity* e_min; // minimum degree meshentity
  apf::MeshEntity* e_i; // current meshentity
  apf::Vector3 coords; // coordinates of current point
  double dist;  // distance from coords to (x,y)
  double min_dist; // minimum distance to (x,y)
  apf::ModelEntity* me_i;
  int me_dimension;

  apf::MeshIterator* it = m_local->begin(0); // iterator over verticies


  // initilize
  e_i = m_local->deref(it);

  e_min = e_i;  // ensure we return a value if conditions are never
                // satisfied

  // calculate distance to (x,y)
  m_local-> getPoint(e_i, 0, coords);
  // no need to take square root if we are only interested in relative
  // distance
  min_dist = (coords[0] - x)*(coords[0] - x) + (coords[1] - y)*(coords[1] - y);



  while ( (e_i = m_local->iterate(it)) )
  {
    me_i = m_local->toModel(e_i);
    me_dimension = m_local->getModelType(me_i);
    if ( !me_dimension ) // if me_dimension == 0
    {

      // calculate distance to (x,y)
      m_local-> getPoint(e_i, 0, coords);
      // no need to take square root if we are only interested in relative
      // distance
      dist = (coords[0] - x)*(coords[0] - x) + (coords[1] - y)*(coords[1] - y);

      if (dist < min_dist)
      {
        e_min = e_i;
        min_dist = dist;
        std::cout << "choosing this vertex" << std::endl;
      }
    }
  }

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
                  std::cout << " edge" << std::endl;
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

// number elements with numbers greater than number of elements to show they 
// have not yet received their final number
void numberElements(apf::Mesh2* m_local, apf::Numbering* elNums, int numEl)
{
  apf::MeshIterator* it = m_local->begin(2);
  apf::MeshEntity* e;
  int k = numEl + 1;

  std::cout << "elNums = " << elNums << std::endl; 
  while ( (e = m_local->iterate(it) ) )
  {
    std::cout << "labelling element " << k - numEl << " as " << k << std::endl;
    std::cout << " e = " << e << std::endl;
    apf::number(elNums, e, 0, 0, k);
    ++k;
  }

}

void printElNumbers(apf::Mesh2*& m_local, apf::Numbering*& elNums)
{
  apf::MeshIterator* it = m_local->begin(2);
  apf::MeshEntity* e;
  int i = 1;
  while ((e = m_local->iterate(it) ) )
  {
//    int num = apf::getNumber(elNums, e, 0 , 0);
//    std::cout << "element " << i << "1 number = " << num << std::endl;
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
// node_statusNumbering tells whether a dof is fixed (=3) or not
// ndof is the number of actual degrees of freedom in the mesh
// comp is the number of dofs per node (the number of components in the dof numberings)
// nodeNums is a numbering over the dofs to be populated with global node numbers
// elNums is numbering over elements (faces) to be populated
// els_reordered is array of pointers to the faces in the new order
// x, y are the coordinates of the point used to deterimine the starting node
//   the mesh vertex classified on a model vertex closest to (x,y) is chosen
void reorder(apf::Mesh2* m_local, int ndof, const int comp, apf::Numbering* node_statusNumbering, apf::Numbering* nodeNums, apf::Numbering* elNums, const double x, const double y)
{
// TODO: move node_statusNumbering checks out one loop level because 
//       it is node status now, not dof status

  std::cout << "Entered reorder" << std::endl;
  apf::FieldShape* fieldshape = m_local->getShape();
  
  std::cout << " field has nodes in dimension: 0 => " << fieldshape->hasNodesIn(0) << std::endl;
  std::cout << " 1 => " << fieldshape->hasNodesIn(1) << std::endl;
  std::cout << " 2 => " << fieldshape->hasNodesIn(2) << std::endl;
  std::cout << " 3 => " << fieldshape->hasNodesIn(3) << std::endl;
 
  // use temporary numbering to count number of nodes
//  apf::Numbering* tmp = apf::numberOwnedNodes(m, "tmpnumber");
//  const int numN = apf::countNodes(tmp);

  // check that we won't overflow integers
  if (ndof > INT_MAX/4)
  {
    std::cerr << "Error: Cannot number this many dofs using C integers" << std::endl;
    return;
  }

  const int numEl = m_local->count(m_local->getDimension());  // counts the number of elements


  // create empty numberings of the proper shape
//  apf::Numbering* elNums = apf::createNumbering(m, "elementNumbers", apf::getConstant(2), 1);
  numberElements(m_local, elNums, numEl);
  std::cout << "finished initial numbering of elements" << std::endl;
//  apf::Numbering* nodeNums = createNumbering(m, "nodeNumbers", m_local->getShape(), 1);
  numberdofs(m_local, nodeNums, ndof, comp);
  std::cout << "finished initial numbering of nodes" << std::endl;



  //print initial numberings
  apf::writeVtkFiles("number_orig", m_local);


  printElNumbers( m_local, elNums);

  // create queues
  std::queue < apf::MeshEntity*> que1;
  std::queue < apf::MeshEntity*> tmpQue;

//  apf::MeshIterator* it = m_local->begin(m_local->getDimension());
//  apf::MeshIterator* nodeIt = m_local->begin(0);

  apf::MeshEntity* e;

  // get starting entity
//  e = m_local->iterate(nodeIt);

  int nodeLabel_i = ndof + 1;  // one-based indexing
  int elementLabel_i = numEl;  // zero-based indexing
  int numNodes_i;

  std::cout << "starting nodelabel+1 = " << nodeLabel_i << std::endl;
  std::cout << "starting elementlabel_i = " << elementLabel_i << std::endl;

  e = getStartEntity(m_local, x, y); // get starting node
  // queue initial entity
  que1.push(e);
  
  while (que1.size() > 0 )
  {
    e = que1.front(); // get next element in queue
    que1.pop();  // remove that element from the que

    numNodes_i = nodeCount(m_local,e); // get number of nodes on this entity

//    std::cout << std::endl;
//    std::cout << "at beginning of while loop" << std::endl;
//    std::cout << "type of current entity = ";
    printType(m_local, e);

    int type_enum_e = m_local->getType(e);

    if ( type_enum_e == apf::Mesh::VERTEX)
    {
//      int num_e = apf::getNumber(nodeNums, e, 0, 0);
 //     std::cout << "current entity number = " << num_e << std::endl;
    } 

    for ( int i = 0; i < numNodes_i; ++i) // label all nodes on this entity
    {
      for ( int c = 0; c < comp; ++c) // loop over dof of the node
      {
//        int nodenum_i = apf::getNumber(nodeNums, e, i, c);
        int dof_status = apf::getNumber(node_statusNumbering, e, i, c);
        if (dof_status >= 2) // if node is free for loaded
        {
//          int nodenum_i = apf::getNumber(nodeNums, e, i, c);
          nodeLabel_i -= 1;  // decrement nodelabel
          number(nodeNums, e, i, c, nodeLabel_i);
//          std::cout << "relabeling node " << nodenum_i << " to " << nodeLabel_i << std::endl;
        } else {
//           std::cout << "found non free node " << nodenum_i << " , status = " << dof_status << std::endl;
           apf::number(nodeNums, e, i, c, 0);
          }

      }
    }

    // if e is a vertex, find adjacencies, look for unlabeled nodes
    if (m_local->getType(e) == apf::Mesh::VERTEX)
    {
      // count number of edges that contain vertex e
      int numEdges_w = m_local->countUpward(e);
      for (int i = 0; i < numEdges_w; ++i)  // loop over edges
      {
//        std::cout << "  on edge " << i << std::endl;
        apf::MeshEntity* edge_i = m_local->getUpward(e, i); // get the edge
//        std::cout << " edge has type" << std::endl;
//        printType(m, edge_i);
        int numFaces_i = m_local->countUpward(edge_i); // get number of faces on the edge

        for (int j = 0; j < numFaces_i; ++j)  // loop over faces
        {
//          std::cout << "    on face " << j << std::endl;
          apf::MeshEntity* face_j = m_local->getUpward(edge_i, j);  // get face
//          std::cout << " face has type";
//          printType(m, face_j);
          int faceNum_j = apf::getNumber(elNums, face_j, 0, 0); // get face number
          int numNodes_j = nodeCount(m_local, face_j); // get number of nodes on face


//          std::cout << "    face number = " << faceNum_j << std::endl;
          // label face if it is not yet labelled
          if (faceNum_j > numEl) // if face not labelled
          {
            elementLabel_i -= 1;  // decrement element label
//            std::cout << "    renumbering face " << faceNum_j << " to " << elementLabel_i << std::endl;
            apf::number(elNums, face_j, 0, 0, elementLabel_i); // number element
          }

          // check if face has nodes that need labelling
          int nodeNum_j;
          for (int k = 0; k < numNodes_j; ++k)
          {
            // get number of node on face
            nodeNum_j = apf::getNumber(nodeNums, face_j, k, 0); 
            if ( (nodeNum_j > ndof) && (nodeNum_j <= 2*ndof)) // node not labelled) and not in queue
            {
//              faceNum_j = apf::getNumber( elNums, face_j, 0, 0); // get new face number
//              std::cout << "adding face " << faceNum_j << " to que" << std::endl; 
//            // double node number to show is has been added to queue
              apf::number(nodeNums, face_j, 0, 0, nodeNum_j*2);
              que1.push(face_j); // add face to que
              break; // only add face to que once
            }
          }


        }  // end face loop

        // look at other vertex on edge
        apf::MeshEntity* otherVertex = apf::getEdgeVertOppositeVert(m_local, edge_i, e);
        int otherVertex_num = getNumber(nodeNums, otherVertex, 0,0);
        bool labeled = (otherVertex_num <= ndof);
        bool queued = (otherVertex_num > 2*ndof);  // vertex already queued

//        std::cout << " other vertex number = " << otherVertex_num << std::endl;

        if (hasNode(m_local, edge_i))
        {
//          std::cout << "edge has node " << std::endl;
          int edgeNode_num = getNumber(nodeNums, edge_i,0,0); // get number of first node on edge
          bool edgeNotLabeled = ( (edgeNode_num > ndof) && ( edgeNode_num <= 2*ndof) ); // edge not labeled nor in que
//          std::cout << "  other vertexlabeled? = " << labeled << " , queued ? = " << queued << " edge nodes not labeled? = " << edgeNotLabeled << std::endl;
          if ((labeled || queued) && edgeNotLabeled)
          {
            // label all nodes on edge
            int numEdgeNodes = nodeCount(m_local, e);
            for (int j = 0; j < numEdgeNodes; ++j)
            {
//              int nodeNum_j = apf::getNumber(nodeNums, edge_i, j, c);
              for (int c = 0; c < comp; ++c) // loop over dofs of node
              {
//                int nodeNum_j = apf::getNumber(nodeNums, edge_i, j, c);
                int dof_status = apf::getNumber(node_statusNumbering, edge_i, j, c);
                if (dof_status >= 2) // if node is free or loaded
                {
                  nodeLabel_i -= 1;
//                  std::cout << " relabelling node " << nodeNum_j << " to " << nodeLabel_i << std::endl;
                  apf::number(nodeNums, edge_i, j, c, nodeLabel_i);
                } else {
//                  std::cout << " found non free node " << nodeNum_j << " , status = " << dof_status << std::endl;
                  apf::number(nodeNums, edge_i, j, c, 0);
                }
              }
            }
          } else { // add edge to queue first, othervertex second (via list)             

              if (edgeNotLabeled) 
              {
//                std::cout << " adding edge to que" << std::endl;
                que1.push(edge_i);
              }

              if ( (!labeled) && (!queued) )
              {
                // double node number to signal this node is queued
                apf::number(nodeNums, otherVertex, 0, 0, 2*otherVertex_num);
//                std::cout << " adding vertex " << otherVertex_num << " to tmpque, doubling number to " << 2*otherVertex_num << std::endl;
                tmpQue.push(otherVertex);
              }
          }          

        } else {
//            if ( (!labeled) )
            if ( (!labeled) && (!queued))  // if otherVertex not labelled
            {
//              std::cout << " edge does not have node" << std::endl;
              // double node number to signal this node queued
              apf::number(nodeNums,otherVertex, 0,0, 2*otherVertex_num);
//              std::cout << " adding vertex " << otherVertex_num << " to tmpque, doubling number to " << 2*otherVertex_num << std::endl;
  
              tmpQue.push(otherVertex);
            }
        }


      }  // end edge loop

      // copy tmpQue into que1, empty tmpQue
      addQueues(que1, tmpQue);


      }  // end if (vertex)
    } // end while loop over que    


  if (elementLabel_i != 0)
  {
    std::cerr << "Warning: element numbering not sane" << std::endl;
    std::cerr << "final elementLabel_i = " << elementLabel_i << std::endl;
  }

  if (nodeLabel_i != 1)
  {
    std::cerr << "Warning: node numbering not sane" << std::endl;
    std::cerr << "final nodeLabel_i = " << nodeLabel_i << std::endl;
  }


  apf::writeVtkFiles("number", m_local);
//  m_local->destroyNative();
//  apf::destroyMesh(m);
//  PCU_Comm_Free();
// MPI_Finalize();
}

