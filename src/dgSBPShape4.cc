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
    {
      &vertex,
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
  {
  public:

    // use shape function value, derivative functions inherited from base EntityShape (which return fail('unimplimented')	  

    // need to register name with PUMI?
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
    } 
  }  // end function getNodeXi
};  // class DG4SBPQuadratic

class DG4SBPCubic : public FieldShape
{
public:
  DG4SBPCubic() { registerSelf(apf::DG4SBPCubic::getName()); }
  const char* getName() const { return "DG4SBPCubic"; }
  class Vertex : public EntityShape
    // use shape function value, derivative functions inherited from base EntityShape (which return fail('unimplimented')	  

    // need to register name with PUMI?
  {

  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const&, NewArray<double>& values) const
    {
      values.allocate(1);
      values[0] = 1.0;
      //            fail("unimplimented getValues called in DG4SBPCubic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called in DG4SBPCubic Vertex");
    }

    int countNodes() const {return 0;}
  };
  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG4SBPCubic Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DG4SBPCubic Edge");
    }

    int countNodes() const {return 0;}


  };
  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG4SBPCubic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG4SBPCubic Triangle");
    }

    int countNodes() const {return 15;}

    void alignSharedNodes(Mesh* m, MeshEntity* elem, MeshEntity* shared, int order[])
    {
      // elem is the triangle 
      // shared is the entity (edge or vertex) being shared
      // order[] contains the mapping such that order[i], where i is the local node number, gives
      // the position of that node in the canonical ordering

      // which is the index of shared in getDownward(elm)
      // rotate is the number of rotations required  ( a complete circle is n rotations, where n is the 
      // number of sides of elem
      // flip determines whether to flip the nodes
    }
  };

  class Tetrahedron : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG4SBPCubic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG4SBPCubic Tetrahdron");
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
      return 15;
    } else
    {
      return 0;
    }

  }

  int getOrder() {return 3;}
  void getNodeXi(int type, int node, Vector3& xi)
  {
    // return the xi coordinates of the specified node of the specified type
    // *in the coordinate system of that type*
    // which makes this function not useful
    if (type != Mesh::TRIANGLE)
    {
      xi = Vector3(0,0,0);
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
      {xi = Vector3(0.5, 0.0, 0.0); break; }
    case 4:
      {xi = Vector3(0.5, 0.5, 0.0); break; }
    case 5:
      {xi = Vector3(0.0, 0.5, 0.0); break; }
    case 6:
      {xi = Vector3(0.20734517566359092, 0.5853096486728182, 0.0); break; }
    case 7:
      {xi = Vector3(0.20734517566359092, 0.20734517566359092, 0.0); break; }
    case 8:
      {xi = Vector3(0.5853096486728182, 0.20734517566359092, 0.0); break; }
    case 9:
      {xi = Vector3(0.17267316464601146, 0.0, 0.0); break; }
    case 10:
      {xi = Vector3(0.8273268353539885, 0.0, 0.0); break; }
    case 11:
      {xi = Vector3(0.8273268353539885, 0.17267316464601146, 0.0); break; }
    case 12:
      {xi = Vector3(0.17267316464601146, 0.8273268353539885, 0.0); break; }
    case 13:
      {xi = Vector3(0.0, 0.8273268353539885, 0.0); break; }
    case 14:
      {xi = Vector3(0.0, 0.17267316464601146, 0.0); break; }
    default:
      {xi = Vector3(0, 0, 0); break; }
    } 
  }
};  // class DG4SBPCubic



class DG4SBPQuartic : public FieldShape
{
public:
  DG4SBPQuartic() { registerSelf(apf::DG4SBPQuartic::getName()); }

  const char* getName() const { return "DG4SBPQuartic"; }


  class Vertex : public EntityShape
  {
    // need to register name with PUMI?
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const&, NewArray<double>& values) const
    {
      fail("unimplimented getValues called in DG4SBPQuartic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called in DG4SBPQuartic Vertex");
    }

    int countNodes() const {return 0;}
  };

  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG4SBPQuartic Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DG4SBPQuartic Edge");
    }

    int countNodes() const {return 0;}


  };
  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG4SBPQuartic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG4SBPQuartic Triangle");
    }

    int countNodes() const {return 24;}

    void alignSharedNodes(Mesh* m, MeshEntity* elem, MeshEntity* shared, int order[])
    {
      // elem is the triangle 
      // shared is the entity (edge or vertex) being shared
      // order[] contains the mapping such that order[i], where i is the local node number, gives
      // the position of that node in the canonical ordering

      // which is the index of shared in getDownward(elm)
      // rotate is the number of rotations required  ( a complete circle is n rotations, where n is the 
      // number of sides of elem
      // flip determines whether to flip the nodes
    }
  };

  class Tetrahedron : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG4SBPQuartic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG4SBPQuartic Tetrahdron");
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
    {
      &vertex,
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
      return 24;
    } else
    {
      return 0;
    }
  }

  int getOrder() {return 4;}

  void getNodeXi(int type, int node, Vector3& xi)
  {
    // return the xi coordinates of the specified node of the specified type
    // *in the coordinate system of that type*
    // which makes this function not useful
    if (type != Mesh::TRIANGLE)
    {
      xi = Vector3(0,0,0);
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
      {xi = Vector3(0.423068519326303, 0.15386296134739397, 0.0); break; }
    case 4:
      {xi = Vector3(0.423068519326303, 0.423068519326303, 0.0); break; }
    case 5:
      {xi = Vector3(0.15386296134739397, 0.423068519326303, 0.0); break; }
    case 6:
      {xi = Vector3(0.11747233803526758, 0.0, 0.0); break; }
    case 7:
      {xi = Vector3(0.8825276619647324, 0.0, 0.0); break; }
    case 8:
      {xi = Vector3(0.8825276619647324, 0.11747233803526758, 0.0); break; }
    case 9:
      {xi = Vector3(0.11747233803526758, 0.8825276619647324, 0.0); break; }
    case 10:
      {xi = Vector3(0.0, 0.8825276619647324, 0.0); break; }
    case 11:
      {xi = Vector3(0.0, 0.11747233803526758, 0.0); break; }
    case 12:
      {xi = Vector3(0.35738424175967753, 0.0, 0.0); break; }
    case 13:
      {xi = Vector3(0.6426157582403225, 0.0, 0.0); break; }
    case 14:
      {xi = Vector3(0.6426157582403225, 0.35738424175967753, 0.0); break; }
    case 15:
      {xi = Vector3(0.35738424175967753, 0.6426157582403225, 0.0); break; }
    case 16:
      {xi = Vector3(0.0, 0.6426157582403225, 0.0); break; }
    case 17:
      {xi = Vector3(0.0, 0.35738424175967753, 0.0); break; }
    case 18:
      {xi = Vector3(0.09915272947160811, 0.1693974424688585, 0.0); break; }
    case 19:
      {xi = Vector3(0.7314498280595334, 0.1693974424688585, 0.0); break; }
    case 20:
      {xi = Vector3(0.7314498280595334, 0.09915272947160811, 0.0); break; }
    case 21:
      {xi = Vector3(0.09915272947160811, 0.7314498280595334, 0.0); break; }
    case 22:
      {xi = Vector3(0.1693974424688585, 0.7314498280595334, 0.0); break; }
    case 23:
      {xi = Vector3(0.1693974424688585, 0.09915272947160811, 0.0); break; }
    default:
      {xi = Vector3(0, 0, 0); break; }
    }
  } 
};  // class DG4SBPQuartic

FieldShape* getDG4SBPShape(int order)
{
  static DG4SBPLinear linear1;
  static DG4SBPQuadratic quadratic1;
  static DG4SBPCubic cubic1;
  static DG4SBPQuartic quartic1;
  // add an if statement here or something to support other orders
  switch (order) {
  case 1:
    return &linear1;
  case 2:
    return &quadratic1;
    case 3:
    return &cubic1;
    case 4:
    return &quartic1;
  default:
    std::cout << "order " << order << " is not supported by apfSBPShape.cc" << std::endl;
    std::cout << std::endl;
    return NULL;
  }
}

} // end namespace apf
