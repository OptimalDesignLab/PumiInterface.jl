//#include "dgSBPShape1.h"
#include "apfMesh.h"
#include "apfVector.h"
#include "apfMatrix.h"
#include "dgSBPShape4.h"

namespace apf {

// sbp-gamma for DG
class DG4SBPLinear : public FieldShape
{
public:
  DG4SBPLinear() { registerSelf(apf::DG4SBPLinear::getName()); }
  const char* getName() const  {return "DG4SBPLinear"; }
  class Vertex : public EntityShape
    // use shape function value, derivative functions inherited from base EntityShape (which return fail('unimplimented')	  

    // need to register name with PUMI?
  {

  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const&, NewArray<double>& values) const
    {
      fail("unimplimented getValues called in DG4SBPLinear Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called DG4SBPLinear Vertex");
    }

    int countNodes() const {return 0;}
  };
  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG4SBPLinear Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DG4SBPLinear Edge");
    }

    int countNodes() const {return 0;}


  };
  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG4SBPLinear Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG4SBPLinear Triangle");
    }

    int countNodes() const {return 6;}

    void alignSharedNodes(Mesh* m, MeshEntity* elem, MeshEntity* shared, int order[])
      // elem is the triangle 
      // shared is the entity (edge or vertex) being shared
      // order[] contains the mapping such that order[i], where i is the local node number, give
      // the position of that node in the canonical ordering
    {
      // nothing to do here for linear element because they have no shared nodes on edges

    }
  };

  class Tetrahedron : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG4SBPLinear Tetrahedron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG4SBPLinear Tetrahedron");
    }

    int countNodes() const {return 0;}

    void alignSharedNodes(Mesh* m, MeshEntity* elem, MeshEntity* shared, int order[])
      // elem is the triangle 
      // shared is the entity (edge or vertex) being shared
      // order[] contains the mapping such that order[i], where i is the local node number, give
      // the position of that node in the canonical ordering
    {
      // nothing to do here for linear element because they have no shared nodes on edges

    }
  };


  EntityShape* getEntityShape(int type)
  {
    static Vertex vertex;
    static Edge edge;
    static Triangle triangle;
    //      static Quad quad;
    static Tetrahedron tet;
    //      static Prism prism;
    //      static Pyramid pyramid;
    //     static Hexahedron hex;
    static EntityShape* shapes[Mesh::TYPES] =
    {&vertex,
      &edge,
      &triangle,
      NULL, // quad
      &tet,
      NULL, // hex
      NULL,  //prism
      NULL};  //pyramid
    return shapes[type];
  }
  bool hasNodesIn(int dimension)
  {
    if (dimension == 2)
      return true;
    else
      return false;
  }
  int countNodesOn(int type)
  {
    if (type == Mesh::TRIANGLE)
      return 6;
    else
      return 0;
  }
  int getOrder() {return 1;}
  void getNodeXi(int type, int node, Vector3& xi)
  {
    // return the xi coordinates of the specified node of the specified type
    // *in the coordinate system of that type*
    // which makes this function not useful
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

  }
};  // class DG4SBPLinear



class DG4SBPQuadratic : public FieldShape
{
public:
  DG4SBPQuadratic() { registerSelf(apf::DG4SBPQuadratic::getName()); }
  const char* getName() const { return "DG4SBPQuadratic"; }
  class Vertex : public EntityShape
    // use shape function value, derivative functions inherited from base EntityShape (which return fail('unimplimented')	  

    // need to register name with PUMI?
  {

  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const&, NewArray<double>& values) const
    {
      //          values.allocate(1);
      //          values[0] = 1.0;
      fail("unimplimented getValues called in DG4SBPQuadratic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getLocalGradients called in DG4SBPQuadratic Vertex");
    }

    int countNodes() const {return 0;}
  };
  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG4SBPQuadratic Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in SBPQuadratiac Edge");
    }

    int countNodes() const {return 0;}


  };
  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG4SBPQuadratic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG4SBPQuadratic Triangle");
    }

    int countNodes() const {return 10;}

    void alignSharedNodes(Mesh* m, MeshEntity* elem, MeshEntity* shared, int order[])
      // elem is the triangle 
      // shared is the entity (edge or vertex) being shared
      // order[] contains the mapping such that order[i], where i is the local node number, give
      // the position of that node in the canonical ordering
    {

    }
  };

  class Tetrahedron : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG4SBPQuadratic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG4SBPQuadratic Tetrahedron");
    }

    int countNodes() const {return 0;}

    void alignSharedNodes(Mesh* m, MeshEntity* elem, MeshEntity* shared, int order[])
      // elem is the triangle 
      // shared is the entity (edge or vertex) being shared
      // order[] contains the mapping such that order[i], where i is the local node number, give
      // the position of that node in the canonical ordering
    {
      // nothing to do here for linear element because they have no shared nodes on edges

    }
  };	

  EntityShape* getEntityShape(int type)
  {
    static Vertex vertex;
    static Edge edge;
    static Triangle triangle;
    //      static Quad quad;
    static Tetrahedron tet;
    //      static Prism prism;
    //      static Pyramid pyramid;
    //     static Hexahedron hex;
    static EntityShape* shapes[Mesh::TYPES] =
    {&vertex,
      &edge,
      &triangle,
      NULL, // quad
      &tet,
      NULL, // hex
      NULL,  //prism
      NULL};  //pyramid
    return shapes[type];
  }

  bool hasNodesIn(int dimension)
  {
    if (dimension == 2)
      return true;
    else
      return false;
  }

  int countNodesOn(int type)
  {
    if (type == Mesh::TRIANGLE)
    {
      return 10;
    } else
    {
      return 0;
    }
  }

  int getOrder() {return 2;}

  void getNodeXi(int type, int node, Vector3& xi)
  {
    // return the xi coordinates of the specified node of the specified type
    // *in the coordinate system of that type*
    // which makes this function not useful
    if (type != Mesh::TRIANGLE)
    {
      xi = Vector3(0, 0, 0);
      return;
    }
    switch (node)
    {
    case 0:
      {xi = Vector3(0.0, 0.0, 0.0); break; }
    case 1:
      {xi = Vector3(1.0, 0.0, 0.0); break; }
    case 2:
      {xi = Vector3(0.0, 1.0, 0.0); break; }
    case 3:
      {xi = Vector3(0.27639320225002106, 0.0, 0.0); break; }
    case 4:
      {xi = Vector3(0.7236067977499789, 0.0, 0.0); break; }
    case 5:
      {xi = Vector3(0.7236067977499789, 0.27639320225002106, 0.0); break; }
    case 6:
      {xi = Vector3(0.27639320225002106, 0.7236067977499789, 0.0); break; }
    case 7:
      {xi = Vector3(0.0, 0.7236067977499789, 0.0); break; }
    case 8:
      {xi = Vector3(0.0, 0.27639320225002106, 0.0); break; }
    case 9:
      {xi = Vector3(0.3333333333333333, 0.3333333333333333, 0.0); break; }
    default:
      {xi = Vector3(0, 0, 0); break; }
    }  // end function getNodeXi
};  // class DG4SBPQuadratic

FieldShape* getDG4SBPShape(int order)
{
  static DG4SBPLinear linear1;
  static DG4SBPQuadratic quadratic1;
  // add an if statement here or something to support other orders
  switch (order) {
  case 1:
    return &linear1;
  case 2:
    return &quadratic1;
  // case 3:
    // return &cubic1;
  // case 4:
    // return &quartic1;
  default:
    std::cout << "order " << order << " is not supported by apfSBPShape.cc" << std::endl;
    std::cout << std::endl;
    return NULL;
  }
}

} // end namespace apf
