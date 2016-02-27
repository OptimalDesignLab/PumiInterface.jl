#include "apfShape.h"
#include "apfMesh.h"
#include "apfVector.h"
#include "apfMatrix.h"
#include "apfSBPShape.h"

namespace apf {

class DG1SBPLinear : public FieldShape
{
  public:
    const char* getName() const  {return "DG1SBPLinear"; }
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
            fail("unimplimented getValues called in DG1SBPLinear Vertex");    
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
           fail("unimplimented getValues called DG1SBPLinear Vertex");
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
          fail("unimplimented getValues() called in DG1SBPLinear Edge");
        }

        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          /*
          grads.allocate(2);
          grads[0] = Vector3(-0.5,0,0);
          grads[1] = Vector3( 0.5,0,0);
          */
          fail("unimplimented getLocaGradients() called in DG1SBPLinear Edge");
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
          fail("unimplimented getValues() called in DG1SBPLinear Triangle");
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
          fail("unimplimented getLocalGradients() called in DG1SBPLinear Triangle");
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
          fail("unimplimented getValues() called in DG1SBPLinear Tetrahedron");
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
          fail("unimplimented getLocalGradients() called in DG1SBPLinear Tetrahedron");
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
          { xi = Vector3(1/6,1/6,0); break; }
        case 1:
          { xi = Vector3(2/3,1/6,0); break; }
        case 2:
          { xi = Vector3(1/6, 2/3,0); break; }
        default: 
          { xi = Vector3(0, 0,0); break; }
      }

    }
};  // class DG1SBPLinear



class DG1SBPQuadratic : public FieldShape
{
  public:
    const char* getName() const { return "DG1SBPQuadratic"; }
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
            fail("unimplimented getValues called in DG1SBPQuadratic Vertex");    
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
           fail("unimplimented getLocalGradients called in DG1SBPQuadratic Vertex");
        }
		
        int countNodes() const {return 1;}
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
          fail("unimplimented getValues() called in DG1SBPQuadratic Edge");
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
		
        int countNodes() const {return 3;}
		

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
          fail("unimplimented getValues() called in DG1SBPQuadratic Triangle");
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
          fail("unimplimented getLocalGradients() called in DG1SBPQuadratic Triangle");
        }
		
        int countNodes() const {return 7;}
		
		void alignSharedNodes(Mesh* m, MeshEntity* elem, MeshEntity* shared, int order[])
		// elem is the triangle 
		// shared is the entity (edge or vertex) being shared
		// order[] contains the mapping such that order[i], where i is the local node number, give
		// the position of that node in the canonical ordering
		{
		  if (m->getType(shared) == Mesh::EDGE)
		  {
		    order[0] = 0;  // no orientation change for quadratic edges
		  }
		  
		  // no need to consider shared vertices
		
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
          fail("unimplimented getValues() called in DG1SBPQuadratic Tetrahdron");
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
          fail("unimplimented getLocalGradients() called in DG1SBPQuadratic Tetrahedron");
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
      if (dimension <= 2)
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if (type == Mesh::VERTEX)
	  {
        return 1;
	  } else if ( type == Mesh::EDGE) 
	  {
	    return 1;
	  } else if ( type == Mesh::TRIANGLE)
	  {
	    return 1;
	  }  else
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
	  if (type == Mesh::VERTEX)
	  {
        xi = Vector3(0,0,0);
	  } else if (type == Mesh::EDGE)
	  {  
	    xi = Vector3(0,0,0);
	  }  else if (type == Mesh::TRIANGLE)
	  {
	    xi = Vector3(1.0/3.0, 1.0/3.0, 0);
	  } else
	  {
	    xi = Vector3(0,0,0);
	  }
	  
    }
};  // class DG1SBPQuadratic




class DG1SBPCubic : public FieldShape
{
  public:
    const char* getName() const { return "DG1SBPCubic"; }
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
//            fail("unimplimented getValues called in DG1SBPCubic Vertex");    
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
           fail("unimplimented getValues called in DG1SBPCubic Vertex");
        }
		
        int countNodes() const {return 1;}
    };
    class Edge : public EntityShape
    {
      public:
	  
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
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
          
  //        fail("unimplimented getValues() called in DG1SBPCubic Edge");
        }

        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          /*
          grads.allocate(2);
          grads[0] = Vector3(-0.5,0,0);
          grads[1] = Vector3( 0.5,0,0);
          */
          fail("unimplimented getLocaGradients() called in DG1SBPCubic Edge");
        }
		
        int countNodes() const {return 4;}
		

    };
    class Triangle : public EntityShape
    {
      public:
	  
        void getValues(Mesh*, MeshEntity*,
            Vector3 const& xi, NewArray<double>& values) const
        {
          
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
          
//          fail("unimplimented getValues() called in DG1SBPCubic Triangle");
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          
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
          
//          fail("unimplimented getLocalGradients() called in DG1SBPCubic Triangle");
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
		int which, rotate;
		bool flip;
		
		getAlignment(m, elem, shared, which, flip, rotate); // populate, which, flip, rotate
		
		  if (m->getType(shared) == Mesh::EDGE)
		  {
		    if (flip)
			{ 
			  order[0] = 1;
			  order[1] = 0;
			} else   // not flipping
			{  
			  order[0] = 0;
			  order[1] = 1;
			}
		  }
		  
		  // no need to consider shared vertices
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
          fail("unimplimented getValues() called in DG1SBPCubic Tetrahdron");
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
          fail("unimplimented getLocalGradients() called in DG1SBPCubic Tetrahdron");
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
      if (dimension <= 2)
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if (type == Mesh::VERTEX)
	  {
        return 1;
	  } else if ( type == Mesh::EDGE) 
	  {
	    return 2;
	  } else if ( type == Mesh::TRIANGLE)
	  {
	    return 3;
	  }  else
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
	  if (type == Mesh::VERTEX)
	  {
        xi = Vector3(0,0,0);
	  } else if (type == Mesh::EDGE)
	  { 
        if (node == 0)
	    {
		  xi = Vector3(-0.4130608881811966, 0, 0);
		} else if (node == 1)
		{
		  xi = Vector3(0.41306088818191977, 0, 0);
		} else  // default case
		{
		  xi = Vector3(0,0, 0);
		}
	  }  else if (type == Mesh::TRIANGLE)
	  {
	    switch(node) {
		  case 0:
		    xi = Vector3(0.20734517566359092, 0.20734517566359092, 0);
			break;
		  case 1:
		    xi = Vector3 (0.5853096486728182, 0.20734517566359092, 0);
			break;
		  case 2: 
		    xi = Vector3( 0.20734517566359092, 0.5853096486728182, 0);
			break;
		  default:
		    xi = Vector3(0, 0, 0);
		  }
	  } else
	  {
	    xi = Vector3(0,0,0);
	  }
	  
    }
};  // class DG1SBPCubic



class DG1SBPQuartic : public FieldShape
{
  public:
    const char* getName() const { return "DG1SBPQuartic"; }
	
	
    class Vertex : public EntityShape
	
	// need to register name with PUMI?
    {
	
      public:
	  
        void getValues(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<double>& values) const
        {
//          values.allocate(1);
//          values[0] = 1.0;
            fail("unimplimented getValues called in DG1SBPQuartic Vertex");    
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
           fail("unimplimented getValues called in DG1SBPQuartic Vertex");
        }
		
        int countNodes() const {return 1;}
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
          fail("unimplimented getValues() called in DG1SBPQuartic Edge");
        }

        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          /*
          grads.allocate(2);
          grads[0] = Vector3(-0.5,0,0);
          grads[1] = Vector3( 0.5,0,0);
          */
          fail("unimplimented getLocaGradients() called in DG1SBPQuartic Edge");
        }
		
        int countNodes() const {return 5;}
		

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
          fail("unimplimented getValues() called in DG1SBPQuartic Triangle");
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
          fail("unimplimented getLocalGradients() called in DG1SBPQuartic Triangle");
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
		int which, rotate;
		bool flip;
		
		getAlignment(m, elem, shared, which, flip, rotate); // populate, which, flip, rotate
		
		  if (m->getType(shared) == Mesh::EDGE)
		  {
		    if (flip)
			{ 
			  order[0] = 2;
			  order[1] = 1;
			  order[2] = 0;
			} else   // not flipping
			{  
			  order[0] = 0;
			  order[1] = 1;
			  order[2] = 2;
			}
		  }
		  
		  // no need to consider shared vertices
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
          fail("unimplimented getValues() called in DG1SBPQuartic Tetrahdron");
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
          fail("unimplimented getLocalGradients() called in DG1SBPQuartic Tetrahdron");
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
      if (dimension <= 2)
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if (type == Mesh::VERTEX)
	  {
        return 1;
	  } else if ( type == Mesh::EDGE) 
	  {
	    return 3;
	  } else if ( type == Mesh::TRIANGLE)
	  {
	    return 6;
	  }  else
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
	  if (type == Mesh::VERTEX)
	  {
        xi = Vector3(0,0,0);
	  } else if (type == Mesh::EDGE)
	  { 
        if (node == 0)
	    {
		  xi = Vector3(-0.5773502691896256, 0, 0);
		} else if (node == 1)
		{
		  xi = Vector3(0.0, 0, 0);
		} else if (node == 2)
		{
		  xi = Vector3(0.5773502691896257, 0, 0);
		} else  // default case
		{
		  xi = Vector3(0,0, 0);
		}
	  }  else if (type == Mesh::TRIANGLE)
	  {
	    switch(node) {
		  case 0:
		    xi = Vector3(0.13079159382974495, 0.13079159382974495, 0);
			break;
		  case 1:
		    xi = Vector3 (0.4247639617258106, 0.1504720765483788, 0);
			break;
		  case 2: 
		    xi = Vector3( 0.7384168123405102, 0.13079159382974495, 0);
			break;
		  case 3:
		    xi = Vector3( 0.4247639617258106, 0.4247639617258106, 0);
			break;
		  case 4: 
		    xi = Vector3(0.13079159382974495, 0.7384168123405102, 0);
			break;
		  case 5:
		    xi = Vector3(0.1504720765483788, 0.4247639617258106, 0);
                    break;
		  default:
		    xi = Vector3(0, 0, 0);
		  }
	  } else
	  {
	    xi = Vector3(0,0,0);
	  }
	  
    }
};  // class DG1SBPQuartic




FieldShape* getDG1SBPShape(int order)
{
  static DG1SBPLinear linear1;
  static DG1SBPQuadratic quadratic1;
  static DG1SBPCubic cubic1;
  static DG1SBPQuartic quartic1;
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
