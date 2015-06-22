#include "apfShape.h"
#include "apfMesh.h"
#include "apfVector.h"
#include "apfMatrix.h"
#include "apfSBPShape.h"

namespace apf {

class SBPLinear : public FieldShape
{
  public:
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "SBPLinear"; }
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
            fail("unimplimented getValues called");    
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
           fail("unimplimented getValues called");
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
          fail("unimplimented getValues() called");
        }

        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          /*
          grads.allocate(2);
          grads[0] = Vector3(-0.5,0,0);
          grads[1] = Vector3( 0.5,0,0);
          */
          fail("unimplimented getLocaGradients() called");
        }
		
        int countNodes() const {return 2;}
		

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
          fail("unimplimented getValues() called");
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
          fail("unimplimented getLocalGradients() called");
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
	

	EntityShape* getEntityShape(int type)
    {
      static Vertex vertex;
      static Edge edge;
      static Triangle triangle;
//      static Quad quad;
//      static Tetrahedron tet;
//      static Prism prism;
//      static Pyramid pyramid;
 //     static Hexahedron hex;
      static EntityShape* shapes[Mesh::TYPES] =
      {&vertex,
       &edge,
       &triangle,
       NULL, // quad
       NULL,
       NULL, // hex
       NULL,  //prism
       NULL};  //pyramid
      return shapes[type];
    }
    bool hasNodesIn(int dimension)
    {
      if (dimension == 0)
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if (type == Mesh::VERTEX)
        return 1;
      else
        return 0;
    }
    int getOrder() {return 1;}
    void getNodeXi(int type, int node, Vector3& xi)
    {
	  // return the xi coordinates of the specified node of the specified type
	  // *in the coordinate system of that type*
	  // which makes this function not useful
      xi = Vector3(0,0,0);
    }
};  // class SBPLinear



class SBPQuadratic : public FieldShape
{
  public:
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "SBPQuadratic"; }
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
            fail("unimplimented getValues called");    
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
           fail("unimplimented getValues called");
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
          fail("unimplimented getValues() called");
        }

        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          /*
          grads.allocate(2);
          grads[0] = Vector3(-0.5,0,0);
          grads[1] = Vector3( 0.5,0,0);
          */
          fail("unimplimented getLocaGradients() called");
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
          fail("unimplimented getValues() called");
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
          fail("unimplimented getLocalGradients() called");
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
	

	EntityShape* getEntityShape(int type)
    {
      static Vertex vertex;
      static Edge edge;
      static Triangle triangle;
//      static Quad quad;
//      static Tetrahedron tet;
//      static Prism prism;
//      static Pyramid pyramid;
 //     static Hexahedron hex;
      static EntityShape* shapes[Mesh::TYPES] =
      {&vertex,
       &edge,
       &triangle,
       NULL, // quad
       NULL,
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
};  // class SBPQuadratic




class SBPCubic : public FieldShape
{
  public:
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "SBPCubic"; }
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
            fail("unimplimented getValues called");    
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
           fail("unimplimented getValues called");
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
          fail("unimplimented getValues() called");
        }

        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          /*
          grads.allocate(2);
          grads[0] = Vector3(-0.5,0,0);
          grads[1] = Vector3( 0.5,0,0);
          */
          fail("unimplimented getLocaGradients() called");
        }
		
        int countNodes() const {return 4;}
		

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
          fail("unimplimented getValues() called");
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
          fail("unimplimented getLocalGradients() called");
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
	

	EntityShape* getEntityShape(int type)
    {
      static Vertex vertex;
      static Edge edge;
      static Triangle triangle;
//      static Quad quad;
//      static Tetrahedron tet;
//      static Prism prism;
//      static Pyramid pyramid;
 //     static Hexahedron hex;
      static EntityShape* shapes[Mesh::TYPES] =
      {&vertex,
       &edge,
       &triangle,
       NULL, // quad
       NULL,
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
		  xi = Vector3(0.29346955590904017, 0, 0);
		} else if (node == 1)
		{
		  xi = Vector3(0.7065304440909599, 0, 0);
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
};  // class SBPCubic



class SBPQuartic : public FieldShape
{
  public:
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "SBPQuartic"; }
	
	
    class Vertex : public EntityShape
	
	// need to register name with PUMI?
    {
	
      public:
	  
        void getValues(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<double>& values) const
        {
//          values.allocate(1);
//          values[0] = 1.0;
            fail("unimplimented getValues called");    
        }
        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>&) const
        {
           fail("unimplimented getValues called");
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
          fail("unimplimented getValues() called");
        }

        void getLocalGradients(Mesh*, MeshEntity*,
            Vector3 const&, NewArray<Vector3>& grads) const
        {
          /*
          grads.allocate(2);
          grads[0] = Vector3(-0.5,0,0);
          grads[1] = Vector3( 0.5,0,0);
          */
          fail("unimplimented getLocaGradients() called");
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
          fail("unimplimented getValues() called");
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
          fail("unimplimented getLocalGradients() called");
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
	

	EntityShape* getEntityShape(int type)
    {
      static Vertex vertex;
      static Edge edge;
      static Triangle triangle;
//      static Quad quad;
//      static Tetrahedron tet;
//      static Prism prism;
//      static Pyramid pyramid;
 //     static Hexahedron hex;
      static EntityShape* shapes[Mesh::TYPES] =
      {&vertex,
       &edge,
       &triangle,
       NULL, // quad
       NULL,
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
		  xi = Vector3(0.2113248654051872, 0, 0);
		} else if (node == 1)
		{
		  xi = Vector3(0.5, 0, 0);
		} else if (node == 2)
		{
		  xi = Vector3(0.7886751345948129, 0, 0);
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
		  default:
		    xi = Vector3(0, 0, 0);
		  }
	  } else
	  {
	    xi = Vector3(0,0,0);
	  }
	  
    }
};  // class SBPQuartic




FieldShape* getSBPShape(int order)
{
  static SBPLinear linear1;
  static SBPQuadratic quadratic1;
  static SBPCubic cubic1;
  static SBPQuartic quartic1;
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
