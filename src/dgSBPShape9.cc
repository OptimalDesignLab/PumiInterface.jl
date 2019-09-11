//#include "dgSBPShape1.h"
#include "apfMesh.h"
#include "apfVector.h"
#include "apfMatrix.h"
#include "dgSBPShape9.h"
#include <cstdio>
#include <cassert>

namespace apf {

// DG version of Pumi's built-in Lagrange types
class DGLagrange : public FieldShape
{
private:  // this should generally be after public:, but this thing has so many
          // internally defined classes it would be impossible to find
  int degree; // degree of fieldshape
  int dim;    // 2D or 3D
  std::string name;
  apf::FieldShape* fshape;  // the CG fieldshape, used for shape functions

public:
  DGLagrange(int _degree, int _dim) : degree(_degree), dim(_dim)
  {
    char namebuf[128];
    int n = sprintf(namebuf, "DGLagrange%d", _degree);
    assert(n <= 128);
    name = namebuf;

    fshape = getLagrange(_degree);
    registerSelf(name.c_str());
  }
  const char* getName() const  {return name.c_str(); }


  // for the EntityShape classes, we use the EntityShape from the CG field
  // when there are nodes on the given entity, and one of the below classes
  // otherwise (which error if the shape function values are request, to
  // indicate something has gone wrong).
  class Vertex : public EntityShape
  {

  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const&, NewArray<double>& values) const
    {
      fail("unimplimented getValues called in DGLagrange Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called DGLagrange Vertex");
    }

    int countNodes() const {return 0;}
  };


  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DGLagrange Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DGLagrange Edge");
    }

    int countNodes() const {return 0;}


  };

  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DGLagrange Triangle");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DGLagrange Triangle");
    }

    int countNodes() const {return 0;}
  };

  class Tetrahedron : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DGLagrange Tetrahedron");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DGLagrange Tetrahedron");
    }

    int countNodes() const {return 0;}
  };


  EntityShape* getEntityShape(int type)
  {
    // thee EntityShape's are the same for 
    static Vertex vertex;
    static Edge edge;
    static Triangle triangle;
    //      static Quad quad;
    static Tetrahedron tet;
    //      static Prism prism;
    //      static Pyramid pyramid;
    //     static Hexahedron hex;
    static EntityShape* shapes[Mesh::TYPES] =
    {
      &vertex,
      &edge,
      &triangle,
      NULL, // quad
      &tet,
      NULL, // hex
      NULL,  //prism
      NULL};  //pyramid

    if (type == Mesh::TRIANGLE && dim == 2)
    {
      return fshape->getEntityShape(Mesh::TRIANGLE);
    } else if (type == Mesh::TET && dim == 3) {
      return fshape->getEntityShape(Mesh::TET);
    } else {
      return shapes[type];
    }
  }
  bool hasNodesIn(int dimension)
  {
    if (dimension == dim)
      return true;
    else
      return false;
  }
  int countNodesOn(int type)
  {
    if (type == Mesh::TRIANGLE || type == Mesh::TET)
      return fshape->getEntityShape(type)->countNodes();
    else
      return 0;
  }
  int getOrder() {return degree;}
  void getNodeXi(int type, int node, Vector3& xi)
  {
    // return the xi coordinates of the specified node of the specified type
    // *in the coordinate system of that type*
    // which makes this function not useful
    //
    //TODO: need to get the Xi coords for the CG lagrange elements
    if (type != Mesh::TRIANGLE)
      fail("getNodeXi() implemented for Triangles only");

    if (degree == 1)
    {
      switch (node)
      {
        case 0:
          {xi = Vector3(0.0, 0.0, 0.0); break; }
        case 1:
          {xi = Vector3(1.0, 0.0, 0.0); break; }
        case 2:
          {xi = Vector3(0.0, 1.0, 0.0); break; }
        default:
          {xi = Vector3(0, 0, 0); break; }
      }
    } else if (degree == 2)
    {
      // the triangular Lagrange quadratic is always Serendipity
      // It appears the SerendepityQuadratic class is only for quads 
      switch (node)
      {
        case 0:
          {xi = Vector3(0.0, 0.0, 0.0); break; }
        case 1:
          {xi = Vector3(1.0, 0.0, 0.0); break; }
        case 2:
          {xi = Vector3(0.0, 1.0, 0.0); break; }
        case 3:
          {xi = Vector3(0.5, 0.0, 0.0); break; }
        case 4:
          {xi = Vector3(0.5, 0.5, 0.0); break; }
        case 5:
          {xi = Vector3(0.0, 0.5, 0.0); break; }
        default:
          {xi = Vector3(0, 0, 0); break; }
      }
    } else if (degree == 3)
    {
      switch (node)
      {
        case 0:  // vertices
          {xi = Vector3(0.0,   0.0,   0.0); break; }
        case 1:
          {xi = Vector3(1.0,   0.0,   0.0); break; }
        case 2:
          {xi = Vector3(0.0,   1.0,   0.0); break; }
        case 3:  // 1st edge
          {xi = Vector3(1.0/3, 0.0,   0.0); break; }
        case 4:
          {xi = Vector3(2.0/3, 0.0,   0.0); break; }
        case 5:  // 2nd edge
          {xi = Vector3(2.0/3, 1./3,  0.0); break; }
        case 6:
          {xi = Vector3(1.0/3, 2.0/3, 0.0); break; }
        case 7:  // 3rd edge
          {xi = Vector3(0,     2.0/3, 0.0); break; }
        case 8:
          {xi = Vector3(0,     1.0/3, 0.0); break; }
        case 9:  // interior
          {xi = Vector3(1.0/3, 1.0/3, 0.0); break; }
        default:
          {xi = Vector3(0, 0, 0); break; }
      }
    } else
    {
      fail("degree > 3 not supported");
    } 

  } // function
};  // class DGLagrange


FieldShape* getDG9SBPShape(int order, int dim)
{
  if (dim != 2)
    fail("DGLagrange only supported in 2D");

  static DGLagrange linear1(1, 2);
  static DGLagrange quadratic1(2, 2);
  static DGLagrange cubic1(3, 2);
  switch (order) {
  case 1:
    return &linear1;
  case 2:
    return &quadratic1;
    case 3:
    return &cubic1;
  default:
    {std::cout << "order " << order << " is not supported by apfSBPShape.cc" << std::endl;
    std::cout << std::endl;
    return NULL;}
  }
}

} // end namespace apf
