//#include "dgSBPShape1.h"
#include "apfMesh.h"
#include "apfVector.h"
#include "apfMatrix.h"
#include "dgSBPShape6.h"

namespace apf {

// sbp-Omega 2p for DG
class DG6SBPLinear : public FieldShape
{
public:
  DG6SBPLinear() { registerSelf(apf::DG6SBPLinear::getName()); }
  const char* getName() const  {return "DG6SBPLinear"; }


  class Vertex : public EntityShape
    // use shape function value, derivative functions inherited from base EntityShape (which return fail('unimplimented')	  

    // need to register name with PUMI?
  {

  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const&, NewArray<double>& values) const
    {
      fail("unimplimented getValues called in DG6SBPLinear Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called DG6SBPLinear Vertex");
    }

    int countNodes() const {return 0;}
  };


  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG6SBPLinear Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DG6SBPLinear Edge");
    }

    int countNodes() const {return 0;}


  };


  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG6SBPLinear Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG6SBPLinear Triangle");
    }

    int countNodes() const {return 3;}

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
      fail("unimplimented getValues() called in DG6SBPLinear Tetrahedron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG6SBPLinear Tetrahedron");
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
      return 3;
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
        {xi = Vector3(0.16666666666666666, 0.6666666666666667, 0.0); break; }
      case 1:
        {xi = Vector3(0.16666666666666666, 0.16666666666666666, 0.0); break; }
      case 2:
        {xi = Vector3(0.6666666666666667, 0.16666666666666666, 0.0); break; }
      default:
        {xi = Vector3(0, 0, 0); break; }
    }

  }
};  // class DG6SBPLinear



class DG6SBPQuadratic : public FieldShape
{
public:
  DG6SBPQuadratic() { registerSelf(apf::DG6SBPQuadratic::getName()); }
  const char* getName() const { return "DG6SBPQuadratic"; }

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
      fail("unimplimented getValues called in DG6SBPQuadratic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getLocalGradients called in DG6SBPQuadratic Vertex");
    }

    int countNodes() const {return 0;}
  };

  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG6SBPQuadratic Edge");
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
      fail("unimplimented getValues() called in DG6SBPQuadratic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG6SBPQuadratic Triangle");
    }

    int countNodes() const {return 6;}

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
      fail("unimplimented getValues() called in DG6SBPQuadratic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG6SBPQuadratic Tetrahedron");
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
      return 6;
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
        {xi = Vector3(0.4459484909159649, 0.10810301816807022, 0.0); break; }
      case 1:
        {xi = Vector3(0.4459484909159649, 0.4459484909159649, 0.0); break; }
      case 2:
        {xi = Vector3(0.10810301816807022, 0.4459484909159649, 0.0); break; }
      case 3:
        {xi = Vector3(0.09157621350977074, 0.8168475729804585, 0.0); break; }
      case 4:
        {xi = Vector3(0.09157621350977074, 0.09157621350977074, 0.0); break; }
      case 5:
        {xi = Vector3(0.8168475729804585, 0.09157621350977074, 0.0); break; }
      default:
        {xi = Vector3(0, 0, 0); break; }
    } 
  }  // end function getNodeXi
};  // class DG6SBPQuadratic

class DG6SBPCubic : public FieldShape
{
public:
  DG6SBPCubic() { registerSelf(apf::DG6SBPCubic::getName()); }
  const char* getName() const { return "DG6SBPCubic"; }
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
      //            fail("unimplimented getValues called in DG6SBPCubic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called in DG6SBPCubic Vertex");
    }

    int countNodes() const {return 0;}
  };
  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG6SBPCubic Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DG6SBPCubic Edge");
    }

    int countNodes() const {return 0;}


  };
  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG6SBPCubic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG6SBPCubic Triangle");
    }

    int countNodes() const {return 12;}

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
      fail("unimplimented getValues() called in DG6SBPCubic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG6SBPCubic Tetrahdron");
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
        {xi = Vector3(0.06308901449150225, 0.8738219710169954, 0.0); break; }
      case 1:
        {xi = Vector3(0.06308901449150225, 0.06308901449150225, 0.0); break; }
      case 2:
        {xi = Vector3(0.8738219710169954, 0.06308901449150225, 0.0); break; }
      case 3:
        {xi = Vector3(0.24928674517091048, 0.501426509658179, 0.0); break; }
      case 4:
        {xi = Vector3(0.24928674517091048, 0.24928674517091048, 0.0); break; }
      case 5:
        {xi = Vector3(0.501426509658179, 0.24928674517091048, 0.0); break; }
      case 6:
        {xi = Vector3(0.31035245103378434, 0.6365024991213988, 0.0); break; }
      case 7:
        {xi = Vector3(0.053145049844816994, 0.6365024991213988, 0.0); break; }
      case 8:
        {xi = Vector3(0.053145049844816994, 0.31035245103378434, 0.0); break; }
      case 9:
        {xi = Vector3(0.31035245103378434, 0.053145049844816994, 0.0); break; }
      case 10:
        {xi = Vector3(0.6365024991213988, 0.053145049844816994, 0.0); break; }
      case 11:
        {xi = Vector3(0.6365024991213988, 0.31035245103378434, 0.0); break; }
      default:
        {xi = Vector3(0, 0, 0); break; }
    } 
  }
};  // class DG6SBPCubic



class DG6SBPQuartic : public FieldShape
{
public:
  DG6SBPQuartic() { registerSelf(apf::DG6SBPQuartic::getName()); }

  const char* getName() const { return "DG6SBPQuartic"; }


  class Vertex : public EntityShape
  {
    // need to register name with PUMI?
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const&, NewArray<double>& values) const
    {
      fail("unimplimented getValues called in DG6SBPQuartic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called in DG6SBPQuartic Vertex");
    }

    int countNodes() const {return 0;}
  };

  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG6SBPQuartic Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DG6SBPQuartic Edge");
    }

    int countNodes() const {return 0;}


  };
  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG6SBPQuartic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG6SBPQuartic Triangle");
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
      fail("unimplimented getValues() called in DG6SBPQuartic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG6SBPQuartic Tetrahdron");
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
      return 18;
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
        {xi = Vector3(0.04599732020598892, 0.9080053595880222, 0.0); break; }
      case 1:
        {xi = Vector3(0.04599732020598892, 0.04599732020598892, 0.0); break; }
      case 2:
        {xi = Vector3(0.9080053595880222, 0.04599732020598892, 0.0); break; }
      case 3:
        {xi = Vector3(0.47737011568388227, 0.045259768632235464, 0.0); break; }
      case 4:
        {xi = Vector3(0.47737011568388227, 0.47737011568388227, 0.0); break; }
      case 5:
        {xi = Vector3(0.045259768632235464, 0.47737011568388227, 0.0); break; }
      case 6:
        {xi = Vector3(0.1768381152202795, 0.646323769559441, 0.0); break; }
      case 7:
        {xi = Vector3(0.1768381152202795, 0.1768381152202795, 0.0); break; }
      case 8:
        {xi = Vector3(0.646323769559441, 0.1768381152202795, 0.0); break; }
      case 9:
        {xi = Vector3(0.40187325969515103, 0.19625348060969794, 0.0); break; }
      case 10:
        {xi = Vector3(0.40187325969515103, 0.40187325969515103, 0.0); break; }
      case 11:
        {xi = Vector3(0.19625348060969794, 0.40187325969515103, 0.0); break; }
      case 12:
        {xi = Vector3(0.030525835755582265, 0.7388316469821417, 0.0); break; }
      case 13:
        {xi = Vector3(0.23064251726227614, 0.7388316469821417, 0.0); break; }
      case 14:
        {xi = Vector3(0.23064251726227614, 0.030525835755582265, 0.0); break; }
      case 15:
        {xi = Vector3(0.030525835755582265, 0.23064251726227614, 0.0); break; }
      case 16:
        {xi = Vector3(0.7388316469821417, 0.23064251726227614, 0.0); break; }
      case 17:
        {xi = Vector3(0.7388316469821417, 0.030525835755582265, 0.0); break; }
      default:
        {xi = Vector3(0, 0, 0); break; }
    }
  } 
};  // class DG6SBPQuartic

FieldShape* getDG6SBPShape(int order)
{
  static DG6SBPLinear linear1;
  static DG6SBPQuadratic quadratic1;
  static DG6SBPCubic cubic1;
  static DG6SBPQuartic quartic1;
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
