//#include "dgSBPShape1.h"
#include "apfMesh.h"
#include "apfVector.h"
#include "apfMatrix.h"
#include "dgSBPShape5.h"

namespace apf {

// sbp-gamma for DG
class DG5SBPLinear : public FieldShape
{
public:
  DG5SBPLinear() { registerSelf(apf::DG5SBPLinear::getName()); }
  const char* getName() const  {return "DG5SBPLinear"; }


  class Vertex : public EntityShape
    // use shape function value, derivative functions inherited from base EntityShape (which return fail('unimplimented')	  

    // need to register name with PUMI?
  {

  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const&, NewArray<double>& values) const
    {
      fail("unimplimented getValues called in DG5SBPLinear Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called DG5SBPLinear Vertex");
    }

    int countNodes() const {return 0;}
  };


  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG5SBPLinear Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DG5SBPLinear Edge");
    }

    int countNodes() const {return 0;}


  };


  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG5SBPLinear Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG5SBPLinear Triangle");
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
      fail("unimplimented getValues() called in DG5SBPLinear Tetrahedron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG5SBPLinear Tetrahedron");
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
        {xi = Vector3(0.21132486540518713, 0.0, 0.0); break; }
      case 1:
        {xi = Vector3(0.7886751345948129, 0.0, 0.0); break; }
      case 2:
        {xi = Vector3(0.7886751345948129, 0.21132486540518713, 0.0); break; }
      case 3:
        {xi = Vector3(0.21132486540518713, 0.7886751345948129, 0.0); break; }
      case 4:
        {xi = Vector3(0.0, 0.7886751345948129, 0.0); break; }
      case 5:
        {xi = Vector3(0.0, 0.21132486540518713, 0.0); break; }
      case 6:
        {xi = Vector3(0.3333333333333333, 0.3333333333333333, 0.0); break; }
      default:
        {xi = Vector3(0, 0, 0); break; }
    }

  }
};  // class DG5SBPLinear



class DG5SBPQuadratic : public FieldShape
{
public:
  DG5SBPQuadratic() { registerSelf(apf::DG5SBPQuadratic::getName()); }
  const char* getName() const { return "DG5SBPQuadratic"; }

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
      fail("unimplimented getValues called in DG5SBPQuadratic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getLocalGradients called in DG5SBPQuadratic Vertex");
    }

    int countNodes() const {return 0;}
  };

  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG5SBPQuadratic Edge");
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
      fail("unimplimented getValues() called in DG5SBPQuadratic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG5SBPQuadratic Triangle");
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
      fail("unimplimented getValues() called in DG5SBPQuadratic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG5SBPQuadratic Tetrahedron");
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
        {xi = Vector3(0.5, 0.0, 0.0); break; }
      case 1:
        {xi = Vector3(0.5, 0.5, 0.0); break; }
      case 2:
        {xi = Vector3(0.0, 0.5, 0.0); break; }
      case 3:
        {xi = Vector3(0.20468064157076207, 0.5906387168584759, 0.0); break; }
      case 4:
        {xi = Vector3(0.20468064157076207, 0.20468064157076207, 0.0); break; }
      case 5:
        {xi = Vector3(0.5906387168584759, 0.20468064157076207, 0.0); break; }
      case 6:
        {xi = Vector3(0.1127016653792583, 0.0, 0.0); break; }
      case 7:
        {xi = Vector3(0.8872983346207417, 0.0, 0.0); break; }
      case 8:
        {xi = Vector3(0.8872983346207417, 0.1127016653792583, 0.0); break; }
      case 9:
        {xi = Vector3(0.1127016653792583, 0.8872983346207417, 0.0); break; }
      case 10:
        {xi = Vector3(0.0, 0.8872983346207417, 0.0); break; }
      case 11:
        {xi = Vector3(0.0, 0.1127016653792583, 0.0); break; }
      default:
        {xi = Vector3(0, 0, 0); break; }
    } 
  }  // end function getNodeXi
};  // class DG5SBPQuadratic

class DG5SBPCubic : public FieldShape
{
public:
  DG5SBPCubic() { registerSelf(apf::DG5SBPCubic::getName()); }
  const char* getName() const { return "DG5SBPCubic"; }
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
      //            fail("unimplimented getValues called in DG5SBPCubic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called in DG5SBPCubic Vertex");
    }

    int countNodes() const {return 0;}
  };
  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG5SBPCubic Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DG5SBPCubic Edge");
    }

    int countNodes() const {return 0;}


  };
  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG5SBPCubic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG5SBPCubic Triangle");
    }

    int countNodes() const {return 21;}

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
      fail("unimplimented getValues() called in DG5SBPCubic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG5SBPCubic Tetrahdron");
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
      return 21;
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
        {xi = Vector3(0.4239374074448131, 0.15212518511037376, 0.0); break; }
      case 1:
        {xi = Vector3(0.4239374074448131, 0.4239374074448131, 0.0); break; }
      case 2:
        {xi = Vector3(0.15212518511037376, 0.4239374074448131, 0.0); break; }
      case 3:
        {xi = Vector3(0.33000947820757187, 0.0, 0.0); break; }
      case 4:
        {xi = Vector3(0.6699905217924281, 0.0, 0.0); break; }
      case 5:
        {xi = Vector3(0.6699905217924281, 0.33000947820757187, 0.0); break; }
      case 6:
        {xi = Vector3(0.33000947820757187, 0.6699905217924281, 0.0); break; }
      case 7:
        {xi = Vector3(0.0, 0.6699905217924281, 0.0); break; }
      case 8:
        {xi = Vector3(0.0, 0.33000947820757187, 0.0); break; }
      case 9:
        {xi = Vector3(0.06943184420297377, 0.0, 0.0); break; }
      case 10:
        {xi = Vector3(0.9305681557970262, 0.0, 0.0); break; }
      case 11:
        {xi = Vector3(0.9305681557970262, 0.06943184420297377, 0.0); break; }
      case 12:
        {xi = Vector3(0.06943184420297377, 0.9305681557970262, 0.0); break; }
      case 13:
        {xi = Vector3(0.0, 0.9305681557970262, 0.0); break; }
      case 14:
        {xi = Vector3(0.0, 0.06943184420297377, 0.0); break; }
      case 15:
        {xi = Vector3(0.10948299428608156, 0.7328086641350939, 0.0); break; }
      case 16:
        {xi = Vector3(0.15770834157882446, 0.7328086641350939, 0.0); break; }
      case 17:
        {xi = Vector3(0.15770834157882446, 0.10948299428608156, 0.0); break; }
      case 18:
        {xi = Vector3(0.10948299428608156, 0.15770834157882446, 0.0); break; }
      case 19:
        {xi = Vector3(0.7328086641350939, 0.15770834157882446, 0.0); break; }
      case 20:
        {xi = Vector3(0.7328086641350939, 0.10948299428608156, 0.0); break; }
      default:
        {xi = Vector3(0, 0, 0); break; }
    } 
  }
};  // class DG5SBPCubic



class DG5SBPQuartic : public FieldShape
{
public:
  DG5SBPQuartic() { registerSelf(apf::DG5SBPQuartic::getName()); }

  const char* getName() const { return "DG5SBPQuartic"; }


  class Vertex : public EntityShape
  {
    // need to register name with PUMI?
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const&, NewArray<double>& values) const
    {
      fail("unimplimented getValues called in DG5SBPQuartic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called in DG5SBPQuartic Vertex");
    }

    int countNodes() const {return 0;}
  };

  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG5SBPQuartic Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocaGradients() called in DG5SBPQuartic Edge");
    }

    int countNodes() const {return 0;}


  };
  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      fail("unimplimented getValues() called in DG5SBPQuartic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG5SBPQuartic Triangle");
    }

    int countNodes() const {return 28;}

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
      fail("unimplimented getValues() called in DG5SBPQuartic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      fail("unimplimented getLocalGradients() called in DG5SBPQuartic Tetrahdron");
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
      return 28;
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
        {xi = Vector3(0.5, 0.0, 0.0); break; }
      case 1:
        {xi = Vector3(0.5, 0.5, 0.0); break; }
      case 2:
        {xi = Vector3(0.0, 0.5, 0.0); break; }
      case 3:
        {xi = Vector3(0.2307653449471585, 0.0, 0.0); break; }
      case 4:
        {xi = Vector3(0.7692346550528415, 0.0, 0.0); break; }
      case 5:
        {xi = Vector3(0.7692346550528415, 0.2307653449471585, 0.0); break; }
      case 6:
        {xi = Vector3(0.2307653449471585, 0.7692346550528415, 0.0); break; }
      case 7:
        {xi = Vector3(0.0, 0.7692346550528415, 0.0); break; }
      case 8:
        {xi = Vector3(0.0, 0.2307653449471585, 0.0); break; }
      case 9:
        {xi = Vector3(0.046910077030668074, 0.0, 0.0); break; }
      case 10:
        {xi = Vector3(0.9530899229693319, 0.0, 0.0); break; }
      case 11:
        {xi = Vector3(0.9530899229693319, 0.046910077030668074, 0.0); break; }
      case 12:
        {xi = Vector3(0.046910077030668074, 0.9530899229693319, 0.0); break; }
      case 13:
        {xi = Vector3(0.0, 0.9530899229693319, 0.0); break; }
      case 14:
        {xi = Vector3(0.0, 0.046910077030668074, 0.0); break; }
      case 15:
        {xi = Vector3(0.5815956036636433, 0.307905177124701, 0.0); break; }
      case 16:
        {xi = Vector3(0.11049921921165569, 0.307905177124701, 0.0); break; }
      case 17:
        {xi = Vector3(0.11049921921165569, 0.5815956036636433, 0.0); break; }
      case 18:
        {xi = Vector3(0.5815956036636433, 0.11049921921165569, 0.0); break; }
      case 19:
        {xi = Vector3(0.307905177124701, 0.11049921921165569, 0.0); break; }
      case 20:
        {xi = Vector3(0.307905177124701, 0.5815956036636433, 0.0); break; }
      case 21:
        {xi = Vector3(0.06242060369267859, 0.8146202065854535, 0.0); break; }
      case 22:
        {xi = Vector3(0.1229591897218679, 0.8146202065854535, 0.0); break; }
      case 23:
        {xi = Vector3(0.1229591897218679, 0.06242060369267859, 0.0); break; }
      case 24:
        {xi = Vector3(0.06242060369267859, 0.1229591897218679, 0.0); break; }
      case 25:
        {xi = Vector3(0.8146202065854535, 0.1229591897218679, 0.0); break; }
      case 26:
        {xi = Vector3(0.8146202065854535, 0.06242060369267859, 0.0); break; }
      case 27:
        {xi = Vector3(0.3333333333333333, 0.3333333333333333, 0.0); break; }
      default:
        {xi = Vector3(0, 0, 0); break; }
    }
  } 
};  // class DG5SBPQuartic

FieldShape* getDG5SBPShape(int order)
{
  static DG5SBPLinear linear1;
  static DG5SBPQuadratic quadratic1;
  static DG5SBPCubic cubic1;
  static DG5SBPQuartic quartic1;
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
