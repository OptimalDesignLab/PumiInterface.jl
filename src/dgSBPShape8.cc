//#include "dgSBPShape1.h"
#include "apfMesh.h"
#include "apfVector.h"
#include "apfMatrix.h"
#include "dgSBPShape8.h"

namespace apf {

// sbp-DiagonalE with min Frobenius norm for DG
class DG8SBPLinear : public FieldShape
{
public:
  DG8SBPLinear() { registerSelf(apf::DG8SBPLinear::getName()); }
  const char* getName() const  {return "DG8SBPLinear"; }


  class Vertex : public EntityShape
    // use shape function value, derivative functions inherited from base EntityShape (which return fail('unimplimented')	  

    // need to register name with PUMI?
  {

  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const&, NewArray<double>& values) const
    {
      fail("unimplimented getValues called in DG8SBPLinear Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called DG8SBPLinear Vertex");
    }

    int countNodes() const {return 0;}
  };


  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG8SBPLinear Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DG8SBPLinear Edge");
    }

    int countNodes() const {return 0;}


  };


  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG8SBPLinear Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG8SBPLinear Triangle");
    }

    int countNodes() const {return 7;}

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
      fail("unimplimented getValues() called in DG8SBPLinear Tetrahedron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG8SBPLinear Tetrahedron");
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
      return 7;
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
      case 6:
        {xi = Vector3(0.3333333333333333, 0.3333333333333333, 0.0); break; }
      default:
        {xi = Vector3(0, 0, 0); break; }
    }

  }
};  // class DG8SBPLinear



class DG8SBPQuadratic : public FieldShape
{
public:
  DG8SBPQuadratic() { registerSelf(apf::DG8SBPQuadratic::getName()); }
  const char* getName() const { return "DG8SBPQuadratic"; }

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
      fail("unimplimented getValues called in DG8SBPQuadratic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getLocalGradients called in DG8SBPQuadratic Vertex");
    }

    int countNodes() const {return 0;}
  };

  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG8SBPQuadratic Edge");
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
      fail("unimplimented getValues() called in DG8SBPQuadratic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG8SBPQuadratic Triangle");
    }

    int countNodes() const {return 12;}

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
      fail("unimplimented getValues() called in DG8SBPQuadratic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG8SBPQuadratic Tetrahedron");
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
      return 12;
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
        {xi = Vector3(0.21285435711180825, 0.5742912857763836, 0.0); break; }
      case 4:
        {xi = Vector3(0.21285435711180825, 0.21285435711180825, 0.0); break; }
      case 5:
        {xi = Vector3(0.5742912857763836, 0.21285435711180825, 0.0); break; }
      case 6:
        {xi = Vector3(0.27639320225002106, 0.0, 0.0); break; }
      case 7:
        {xi = Vector3(0.7236067977499789, 0.0, 0.0); break; }
      case 8:
        {xi = Vector3(0.7236067977499789, 0.27639320225002106, 0.0); break; }
      case 9:
        {xi = Vector3(0.27639320225002106, 0.7236067977499789, 0.0); break; }
      case 10:
        {xi = Vector3(0.0, 0.7236067977499789, 0.0); break; }
      case 11:
        {xi = Vector3(0.0, 0.27639320225002106, 0.0); break; }
      default:
        {xi = Vector3(0, 0, 0); break; }
    } 
  }  // end function getNodeXi
};  // class DG8SBPQuadratic

class DG8SBPCubic : public FieldShape
{
public:
  DG8SBPCubic() { registerSelf(apf::DG8SBPCubic::getName()); }
  const char* getName() const { return "DG8SBPCubic"; }
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
      //            fail("unimplimented getValues called in DG8SBPCubic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called in DG8SBPCubic Vertex");
    }

    int countNodes() const {return 0;}
  };
  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG8SBPCubic Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DG8SBPCubic Edge");
    }

    int countNodes() const {return 0;}


  };
  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG8SBPCubic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG8SBPCubic Triangle");
    }

    int countNodes() const {return 18;}

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
      fail("unimplimented getValues() called in DG8SBPCubic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG8SBPCubic Tetrahdron");
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
      return 18;
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
        {xi = Vector3(0.4243860251718814, 0.1512279496562372, 0.0); break; }
      case 7:
        {xi = Vector3(0.4243860251718814, 0.4243860251718814, 0.0); break; }
      case 8:
        {xi = Vector3(0.1512279496562372, 0.4243860251718814, 0.0); break; }
      case 9:
        {xi = Vector3(0.14200508409677795, 0.7159898318064442, 0.0); break; }
      case 10:
        {xi = Vector3(0.14200508409677795, 0.14200508409677795, 0.0); break; }
      case 11:
        {xi = Vector3(0.7159898318064442, 0.14200508409677795, 0.0); break; }
      case 12:
        {xi = Vector3(0.17267316464601146, 0.0, 0.0); break; }
      case 13:
        {xi = Vector3(0.8273268353539885, 0.0, 0.0); break; }
      case 14:
        {xi = Vector3(0.8273268353539885, 0.17267316464601146, 0.0); break; }
      case 15:
        {xi = Vector3(0.17267316464601146, 0.8273268353539885, 0.0); break; }
      case 16:
        {xi = Vector3(0.0, 0.8273268353539885, 0.0); break; }
      case 17:
        {xi = Vector3(0.0, 0.17267316464601146, 0.0); break; }
      default:
        {xi = Vector3(0, 0, 0); break; }
    } 
  }
};  // class DG8SBPCubic



class DG8SBPQuartic : public FieldShape
{
public:
  DG8SBPQuartic() { registerSelf(apf::DG8SBPQuartic::getName()); }

  const char* getName() const { return "DG8SBPQuartic"; }


  class Vertex : public EntityShape
  {
    // need to register name with PUMI?
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const&, NewArray<double>& values) const
    {
      fail("unimplimented getValues called in DG8SBPQuartic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called in DG8SBPQuartic Vertex");
    }

    int countNodes() const {return 0;}
  };

  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG8SBPQuartic Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DG8SBPQuartic Edge");
    }

    int countNodes() const {return 0;}


  };
  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG8SBPQuartic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG8SBPQuartic Triangle");
    }

    int countNodes() const {return 27;}

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
      fail("unimplimented getValues() called in DG8SBPQuartic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG8SBPQuartic Tetrahdron");
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
      return 27;
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
        {xi = Vector3(0.10367750814280517, 0.7926449837143896, 0.0); break; }
      case 4:
        {xi = Vector3(0.10367750814280517, 0.10367750814280517, 0.0); break; }
      case 5:
        {xi = Vector3(0.7926449837143896, 0.10367750814280517, 0.0); break; }
      case 6:
        {xi = Vector3(0.2653313804842097, 0.46933723903158064, 0.0); break; }
      case 7:
        {xi = Vector3(0.2653313804842097, 0.2653313804842097, 0.0); break; }
      case 8:
        {xi = Vector3(0.46933723903158064, 0.2653313804842097, 0.0); break; }
      case 9:
        {xi = Vector3(0.35738424175967753, 0.0, 0.0); break; }
      case 10:
        {xi = Vector3(0.6426157582403225, 0.0, 0.0); break; }
      case 11:
        {xi = Vector3(0.6426157582403225, 0.35738424175967753, 0.0); break; }
      case 12:
        {xi = Vector3(0.35738424175967753, 0.6426157582403225, 0.0); break; }
      case 13:
        {xi = Vector3(0.0, 0.6426157582403225, 0.0); break; }
      case 14:
        {xi = Vector3(0.0, 0.35738424175967753, 0.0); break; }
      case 15:
        {xi = Vector3(0.11747233803526758, 0.0, 0.0); break; }
      case 16:
        {xi = Vector3(0.8825276619647324, 0.0, 0.0); break; }
      case 17:
        {xi = Vector3(0.8825276619647324, 0.11747233803526758, 0.0); break; }
      case 18:
        {xi = Vector3(0.11747233803526758, 0.8825276619647324, 0.0); break; }
      case 19:
        {xi = Vector3(0.0, 0.8825276619647324, 0.0); break; }
      case 20:
        {xi = Vector3(0.0, 0.11747233803526758, 0.0); break; }
      case 21:
        {xi = Vector3(0.5870855671333673, 0.0882739606015811, 0.0); break; }
      case 22:
        {xi = Vector3(0.3246404722650515, 0.0882739606015811, 0.0); break; }
      case 23:
        {xi = Vector3(0.3246404722650515, 0.5870855671333673, 0.0); break; }
      case 24:
        {xi = Vector3(0.5870855671333673, 0.3246404722650515, 0.0); break; }
      case 25:
        {xi = Vector3(0.0882739606015811, 0.3246404722650515, 0.0); break; }
      case 26:
        {xi = Vector3(0.0882739606015811, 0.5870855671333673, 0.0); break; }
      default:
        {xi = Vector3(0, 0, 0); break; }
    }
  } 
};  // class DG8SBPQuartic

FieldShape* getDG8SBPShape(int order)
{
  static DG8SBPLinear linear1;
  static DG8SBPQuadratic quadratic1;
  static DG8SBPCubic cubic1;
  static DG8SBPQuartic quartic1;
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
