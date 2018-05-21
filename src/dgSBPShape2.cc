//#include "dgSBPShape1.h"
#include "apfMesh.h"
#include "apfVector.h"
#include "apfMatrix.h"
#include "dgSBPShape2.h"

namespace apf {

// sbp-gamma for DG
class DG2SBPLinear : public FieldShape
{
  public:
    DG2SBPLinear() { registerSelf(apf::DG2SBPLinear::getName()); }
    const char* getName() const  {return "DG2SBPLinear"; }
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
            fail("unimplimented getValues called in DG2SBPLinear Vertex");    
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
           fail("unimplimented getValues called DG2SBPLinear Vertex");
        }
		
        int countNodes() const {return 0;}
    };
    class Edge : public EntityShape
    {
      public:
	  
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          /*
          values.allocate(2);
          values[0] = (1.0-xi[0])/2.0;
          values[1] = (1.0+xi[0])/2.0;
          */
          fail("unimplimented getValues() called in DG2SBPLinear Edge");
        }

        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          /*
          grads.allocate(2);
          grads[0] = Vector3(-0.5,0,0);
          grads[1] = Vector3( 0.5,0,0);
          */
          fail("unimplimented getLocaGradients() called in DG2SBPLinear Edge");
        }
		
        int countNodes() const {return 0;}
		

    };
    class Triangle : public EntityShape
    {
      public:
	  
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          /*
          values.allocate(3);
          values[0] = 1-xi[0]-xi[1];
          values[1] = xi[0];
          values[2] = xi[1];
          */
          fail("unimplimented getValues() called in DG2SBPLinear Triangle");
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          /*
          grads.allocate(3);
          grads[0] = Vector3(-1,-1,0);
          grads[1] = Vector3( 1, 0,0);
          grads[2] = Vector3( 0, 1,0);
          */
          fail("unimplimented getLocalGradients() called in DG2SBPLinear Triangle");
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
          /*
          values.allocate(3);
          values[0] = 1-xi[0]-xi[1];
          values[1] = xi[0];
          values[2] = xi[1];
          */
          fail("unimplimented getValues() called in DG2SBPLinear Tetrahedron");
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          /*
          grads.allocate(3);
          grads[0] = Vector3(-1,-1,0);
          grads[1] = Vector3( 1, 0,0);
          grads[2] = Vector3( 0, 1,0);
          */
          fail("unimplimented getLocalGradients() called in DG2SBPLinear Tetrahedron");
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
      {xi = Vector3(0.0, 0.0, 0.0); break; }
    case 1:
      {xi = Vector3(1.0, 0.0, 0.0); break; }
    case 2:
      {xi = Vector3(0.0, 1.0, 0.0); break; }
    default:
      {xi = Vector3(0, 0, 0); break; }
    }

  }
};  // class DG2SBPLinear



class DG2SBPQuadratic : public FieldShape
{
public:
  DG2SBPQuadratic() { registerSelf(apf::DG2SBPQuadratic::getName()); }
  const char* getName() const { return "DG2SBPQuadratic"; }
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
      fail("unimplimented getValues called in DG2SBPQuadratic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getLocalGradients called in DG2SBPQuadratic Vertex");
    }

    int countNodes() const {return 0;}
  };
  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      /*
         values.allocate(2);
         values[0] = (1.0-xi[0])/2.0;
         values[1] = (1.0+xi[0])/2.0;
         */
      fail("unimplimented getValues() called in DG2SBPQuadratic Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      /*
         grads.allocate(2);
         grads[0] = Vector3(-0.5,0,0);
         grads[1] = Vector3( 0.5,0,0);
         */
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
      /*
         values.allocate(3);
         values[0] = 1-xi[0]-xi[1];
         values[1] = xi[0];
         values[2] = xi[1];
         */
      fail("unimplimented getValues() called in DG2SBPQuadratic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      /*
         grads.allocate(3);
         grads[0] = Vector3(-1,-1,0);
         grads[1] = Vector3( 1, 0,0);
         grads[2] = Vector3( 0, 1,0);
         */
      fail("unimplimented getLocalGradients() called in DG2SBPQuadratic Triangle");
    }

    int countNodes() const {return 7;}

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
      /*
         values.allocate(3);
         values[0] = 1-xi[0]-xi[1];
         values[1] = xi[0];
         values[2] = xi[1];
         */
      fail("unimplimented getValues() called in DG2SBPQuadratic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      /*
         grads.allocate(3);
         grads[0] = Vector3(-1,-1,0);
         grads[1] = Vector3( 1, 0,0);
         grads[2] = Vector3( 0, 1,0);
         */
      fail("unimplimented getLocalGradients() called in DG2SBPQuadratic Tetrahedron");
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
      return 7;
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

  }  // end function getNodeXi
};  // class DG2SBPQuadratic




class DG2SBPCubic : public FieldShape
{
public:
  DG2SBPCubic() { registerSelf(apf::DG2SBPCubic::getName()); }
  const char* getName() const { return "DG2SBPCubic"; }
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
      //            fail("unimplimented getValues called in DG2SBPCubic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called in DG2SBPCubic Vertex");
    }

    int countNodes() const {return 0;}
  };
  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      /*
         double xiv = xi[0];
         const double n1 = 0.0;
         const double n2 = 0.2934695559090417;
         const double n3 = 0.706530440905599;
         const double n4 = 1.0;
         values.allocate(4);
         values[0] = (xiv - n2)*(xiv - n3)*(xiv - n4)/(-0.20734517566347355);
         values[1] = (xiv - n1)*(xiv - n3)*(xiv - n4)/(0.08564618241975609);
         values[2] = (xiv - n1)*(xiv - n2)*(xiv - n4)/(-0.08564618241982433);
         values[3] = (xiv - n1)*(xiv - n2)*(xiv - n3)/(0.2073451756638735);
         */

      fail("unimplimented getValues() called in DG2SBPCubic Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      /*
         grads.allocate(2);
         grads[0] = Vector3(-0.5,0,0);
         grads[1] = Vector3( 0.5,0,0);
         */
      fail("unimplimented getLocaGradients() called in DG2SBPCubic Edge");
    }

    int countNodes() const {return 0;}


  };
  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      /*          
                  values.allocate(12);
                  values[0] = 1-xi[0]-xi[1];
                  values[1] = xi[0];
                  values[2] = xi[1];
                  values[3] = 0.0;

                  values[4] = 0.0;
                  values[5] = 0.0;
                  values[6] = 0.0;
                  values[7] = 0.0;
                  values[8] = 0.0;
                  values[9] = 0.0;
                  values[10] = 0.0;
                  values[11] = 0.0;
                  */

      fail("unimplimented getValues() called in DG2SBPCubic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      /*          
                  grads.allocate(12);
                  grads[0] = Vector3(-1,-1,0);
                  grads[1] = Vector3( 1, 0,0);
                  grads[2] = Vector3( 0, 1,0);
                  grads[3] = Vector3( 0, 0, 0);
                  grads[4] = Vector3( 0, 0, 0);
                  grads[5] = Vector3( 0, 0, 0);
                  grads[6] = Vector3( 0, 0, 0);
                  grads[7] = Vector3( 0, 0, 0);
                  grads[8] = Vector3( 0, 0, 0);
                  grads[9] = Vector3( 0, 0, 0);
                  grads[10] = Vector3( 0, 0, 0);
                  grads[11] = Vector3( 0, 0, 0);
                  */          
      fail("unimplimented getLocalGradients() called in DG2SBPCubic Triangle");
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
      /*
         values.allocate(3);
         values[0] = 1-xi[0]-xi[1];
         values[1] = xi[0];
         values[2] = xi[1];
         */
      fail("unimplimented getValues() called in DG2SBPCubic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      /*
         grads.allocate(3);
         grads[0] = Vector3(-1,-1,0);
         grads[1] = Vector3( 1, 0,0);
         grads[2] = Vector3( 0, 1,0);
         */
      fail("unimplimented getLocalGradients() called in DG2SBPCubic Tetrahdron");
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
      {xi = Vector3(0.0, 0.0, 0.0); break; }
    case 1:
      {xi = Vector3(1.0, 0.0, 0.0); break; }
    case 2:
      {xi = Vector3(0.0, 1.0, 0.0); break; }
    case 3:
      {xi = Vector3(0.20734517566359092, 0.5853096486728182, 0.0); break; }
    case 4:
      {xi = Vector3(0.20734517566359092, 0.20734517566359092, 0.0); break; }
    case 5:
      {xi = Vector3(0.5853096486728182, 0.20734517566359092, 0.0); break; }
    case 6:
      {xi = Vector3(0.7065304440909599, 0.0, 0.0); break; }
    case 7:
      {xi = Vector3(0.29346955590904017, 0.0, 0.0); break; }
    case 8:
      {xi = Vector3(0.29346955590904017, 0.7065304440909599, 0.0); break; }
    case 9:
      {xi = Vector3(0.7065304440909599, 0.29346955590904017, 0.0); break; }
    case 10:
      {xi = Vector3(0.0, 0.29346955590904017, 0.0); break; }
    case 11:
      {xi = Vector3(0.0, 0.7065304440909599, 0.0); break; }
    default:
      {xi = Vector3(0, 0, 0); break; }
    } 
  }
};  // class DG2SBPCubic



class DG2SBPQuartic : public FieldShape
{
public:
  DG2SBPQuartic() { registerSelf(apf::DG2SBPQuartic::getName()); }
  const char* getName() const { return "DG2SBPQuartic"; }


  class Vertex : public EntityShape

    // need to register name with PUMI?
  {

  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const&, NewArray<double>& values) const
    {
      //          values.allocate(1);
      //          values[0] = 1.0;
      fail("unimplimented getValues called in DG2SBPQuartic Vertex");    
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>&) const
    {
      fail("unimplimented getValues called in DG2SBPQuartic Vertex");
    }

    int countNodes() const {return 0;}
  };
  class Edge : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      /*
         values.allocate(2);
         values[0] = (1.0-xi[0])/2.0;
         values[1] = (1.0+xi[0])/2.0;
         */
      fail("unimplimented getValues() called in DG2SBPQuartic Edge");
    }

    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      /*
         grads.allocate(2);
         grads[0] = Vector3(-0.5,0,0);
         grads[1] = Vector3( 0.5,0,0);
         */
      fail("unimplimented getLocaGradients() called in DG2SBPQuartic Edge");
    }

    int countNodes() const {return 0;}


  };
  class Triangle : public EntityShape
  {
  public:

    void getValues(Mesh*, MeshEntity*,
                   Vector3 const& xi, NewArray<double>& values) const
    {
      /*
         values.allocate(3);
         values[0] = 1-xi[0]-xi[1];
         values[1] = xi[0];
         values[2] = xi[1];
         */
      fail("unimplimented getValues() called in DG2SBPQuartic Triangle");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      /*
         grads.allocate(3);
         grads[0] = Vector3(-1,-1,0);
         grads[1] = Vector3( 1, 0,0);
         grads[2] = Vector3( 0, 1,0);
         */
      fail("unimplimented getLocalGradients() called in DG2SBPQuartic Triangle");
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
      /*
         values.allocate(3);
         values[0] = 1-xi[0]-xi[1];
         values[1] = xi[0];
         values[2] = xi[1];
         */
      fail("unimplimented getValues() called in DG2SBPQuartic Tetrahdron");
    }
    void getLocalGradients(Mesh*, MeshEntity*,
                           Vector3 const&, NewArray<Vector3>& grads) const
    {
      /*
         grads.allocate(3);
         grads[0] = Vector3(-1,-1,0);
         grads[1] = Vector3( 1, 0,0);
         grads[2] = Vector3( 0, 1,0);
         */
      fail("unimplimented getLocalGradients() called in DG2SBPQuartic Tetrahdron");
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
      {xi = Vector3(0.5, 0.0, 0.0); break; }
    case 4:
      {xi = Vector3(0.5, 0.5, 0.0); break; }
    case 5:
      {xi = Vector3(0.0, 0.5, 0.0); break; }
    case 6:
      {xi = Vector3(0.13079159382974495, 0.7384168123405102, 0.0); break; }
    case 7:
      {xi = Vector3(0.13079159382974495, 0.13079159382974495, 0.0); break; }
    case 8:
      {xi = Vector3(0.7384168123405102, 0.13079159382974495, 0.0); break; }
    case 9:
      {xi = Vector3(0.4247639617258106, 0.1504720765483788, 0.0); break; }
    case 10:
      {xi = Vector3(0.4247639617258106, 0.4247639617258106, 0.0); break; }
    case 11:
      {xi = Vector3(0.1504720765483788, 0.4247639617258106, 0.0); break; }
    case 12:
      {xi = Vector3(0.7886751345948129, 0.0, 0.0); break; }
    case 13:
      {xi = Vector3(0.2113248654051872, 0.0, 0.0); break; }
    case 14:
      {xi = Vector3(0.2113248654051872, 0.7886751345948129, 0.0); break; }
    case 15:
      {xi = Vector3(0.7886751345948129, 0.2113248654051872, 0.0); break; }
    case 16:
      {xi = Vector3(0.0, 0.2113248654051872, 0.0); break; }
    case 17:
      {xi = Vector3(0.0, 0.7886751345948129, 0.0); break; }
    default:
      {xi = Vector3(0, 0, 0); break; }
    }
  }
};  // class DG2SBPQuartic




FieldShape* getDG2SBPShape(int order)
{
  static DG2SBPLinear linear1;
  static DG2SBPQuadratic quadratic1;
  static DG2SBPCubic cubic1;
  static DG2SBPQuartic quartic1;
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
    return NULL;
  }
}

} // end namespace apf
