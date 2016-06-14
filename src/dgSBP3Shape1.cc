#include "apfShape.h"
#include "apfMesh.h"
#include "apfVector.h"
#include "apfMatrix.h"
#include "dgSBPShape1.h"

namespace apf {

class DG1SBP3Linear : public FieldShape
{
  public:
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "dg1SBP3Linear"; }

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
    
        int countNodes() const {return 4;}
    
        void alignSharedNodes(Mesh* m, MeshEntity* elem, MeshEntity* shared, int order[])
        // elem is the triangle 
        // shared is the entity (edge or vertex) being shared
        // order[] contains the mapping such that order[i], where i is the local node number, give
        // the position of that node in the canonical ordering
        {
          // nothing to do here for linear element because they have no shared nodes on edges
        
        }
    };  // class Tetrahedron
	

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
      if (dimension == 3)
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if (type == Mesh::TET)
        return 4;
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
          { xi = Vector3(1.0/8,1.0/8,1.0/8); break; }
        case 1:
          { xi = Vector3(5.0/8,1.0/8,1.0/8); break; }
        case 2:
          { xi = Vector3(1.0/8, 5.0/8,1.0/8); break; }
        case 3:
          { xi = Vector3(1.0/8, 1.0/8, 5.0/8); break; }
        default: 
          { xi = Vector3(0, 0,0); break; }
      }

    }
};  // class SBPLinear



class DG1SBP3Quadratic : public FieldShape
{
  public:
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "SBP3Quadratic"; }
	
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
		
        int countNodes() const {return 6;}
// no nodes on face, so no need to align them		
/*		
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
*/		
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
    
        int countNodes() const {return 11;}
    
    void alignSharedNodes(Mesh* m, MeshEntity* elem, MeshEntity* shared, int order[])
    // elem is the triangle 
    // shared is the entity (edge or vertex) being shared
    // order[] contains the mapping such that order[i], where i is the local node number, give
    // the position of that node in the canonical ordering
    {
      // vertex does not have orientation, only one node on edge
	  if (m->getType(shared) == Mesh::EDGE)
	  {
		  order[0] = 0;
	  }
    
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
      if (dimension <= 1 || dimension == 3)
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
	    return 0;
	  }  else if ( type == Mesh::TET)
	  {
		  return 1;
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
	  // which makes this function not useful, because the user could define the origins differently
	  if (type == Mesh::VERTEX)
	  {
        xi = Vector3(0,0,0);
	  } else if (type == Mesh::EDGE)
	  {  
	    xi = Vector3(0,0,0);
	  }  else if (type == Mesh::TET)
	  {
	    xi = Vector3(0.25, 0.25, 0.25);
	  } else
	  {
	    xi = Vector3(0,0,0);
	  }
	  
    }
};  // class SBPQuadratic




class DG1SBP3Cubic : public FieldShape
{
  public:
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "SBP3Cubic"; }
	
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
		
        int countNodes() const {return 10;}
		
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
			
			// edges cannot be rotated
		  }
		  
		  if (m->getType(shared) == Mesh::TRIANGLE)
		  {
			  order[0] = 0;  // only one node classified on the face
		  }
			  
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
      if (dimension <= 3)
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
	    return 1;
	  }  else if ( type == Mesh::TET)
	  {
		  return 4;
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
	  if (type == Mesh::VERTEX)
	  {
        xi = Vector3(0,0,0);
	  } else if (type == Mesh::EDGE)
	  { 
        if (node == 0)
	    {
		  xi = Vector3(0.30480589839889616,0.0,0.0);
		} else if (node == 1)
		{
		  xi = Vector3(0.6951941016011038, 0, 0);
		} else  // default case
		{
		  xi = Vector3(0,0, 0);
		}
	  }  else if (type == Mesh::TRIANGLE)
	  {
		  if (node == 0)
		  {
			  xi = Vector3(1.0/3.0, 1.0/3.0, 0);
		  } else  // default case
		  {
			  xi = Vector3(0.0, 0.0, 0.0);
		  }
	  } else if (type == Mesh::TET)
	  {
	     switch(node) {
		  case 0:
		    xi = Vector3(0.1524029491994481,0.1524029491994481,0.1524029491994481);
			break;
		  case 1:
		    xi = Vector3 (0.5427911524016557,0.1524029491994481,0.1524029491994481);
			break;
		  case 2: 
		    xi = Vector3( 0.1524029491994481,0.5427911524016557,0.1524029491994481);
			break;
		  case 3:
		    xi = Vector3(0.1524029491994481,0.1524029491994481,0.5427911524016557);
                    break;
		  default:
		    xi = Vector3(0, 0, 0);
		  }
	  }
	  
    }
};  // class SBPCubic



class DG1SBP3Quartic : public FieldShape
{
  public:
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "SBP3Quartic"; }
	
	
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
		
        int countNodes() const {return 15;}
		
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
    
        int countNodes() const {return 45;}
    
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
			
			// edges cannot be flipped
		  }
		  
		  if (m->getType(shared) == Mesh::TRIANGLE)
		  {
			  if (flip)  // reverse order of all nodes
			  {
				  order[0] = 2;
				  order[1] = 1;
				  order[2] = 0;
			  }
			  
			  if (rotate)
			  {				  
				  for (int i = 0; i < 3; ++i)
				  {
					order[i] += rotate;  // this works because therea re only 3 nodes on  a face
					  
					if (order[i] > 2)
					{
					  order[i] -= 3;
					}
				  }
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
      if (dimension <= 3)
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
	    return 3;
	  }  else if ( type == Mesh::TET)
	  {
		return 11;
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
	  if (type == Mesh::VERTEX)
	  {
        xi = Vector3(0,0,0);
	  } else if (type == Mesh::EDGE)
	  { 
        if (node == 0)
	    {
		  xi = Vector3(0.25737274681480826, 0, 0);
		} else if (node == 1)
		{
		  xi = Vector3(0.5, 0, 0);
		} else if (node == 2)
		{
		  xi = Vector3(0.7426272531851917, 0, 0);
		} else  // default case
		{
		  xi = Vector3(0,0, 0);
		}
	  }  else if (type == Mesh::TRIANGLE)
	  {
	    switch(node) {
		  case 0:
		    xi = Vector3(0.22504424155412348,0.22504424155412348,0.0);
			break;
		  case 1:
		    xi = Vector3 (0.549911516891753,0.22504424155412348,0.0);
			break;
		  case 2: 
		    xi = Vector3( 0.22504424155412348,0.549911516891753,0.0);
			break;
		  default:
		    xi = Vector3(0, 0, 0);
		  }
	  } else if (type == Mesh::TET)
	  {
		switch(node) {
	      case 0:
		    xi = Vector3(0.09472900091823398,0.09472900091823398,0.09472900091823398);
		    break;
		  case 1:
		    xi = Vector3(0.715812997245298,0.09472900091823398,0.09472900091823398);
			break;
		  case 2:
		    xi = Vector3( 0.09472900091823398,0.715812997245298,0.09472900091823398);
			break;
		  case 3:
		    xi = Vector3( 0.09472900091823398,0.09472900091823398,0.715812997245298);
			break;
		  case 4:
		    xi = Vector3( 0.39128583990222227,0.10871416009777772,0.10871416009777772);
			break;
		  case 5:
		    xi = Vector3(0.39128583990222227,0.39128583990222227,0.10871416009777772);
			break;
		  case 6:
		    xi = Vector3(0.10871416009777772,0.39128583990222227,0.10871416009777772);
			break;
		  case 7:
		    xi = Vector3(0.39128583990222227,0.10871416009777772,0.39128583990222227);
			break;
		  case 8:
		    xi = Vector3(0.10871416009777772,0.10871416009777772,0.39128583990222227);
			break;
		  case 9:
		    xi = Vector3(0.10871416009777772,0.39128583990222227,0.39128583990222227);
			break;
		  case 10:
		    xi = Vector3( 0.25, 0.25, 0.25);
			break;
		  default:
		    xi = Vector3(0, 0, 0);

	    }  // end switch block
	  }  // end if type == tet
	  
    }  // end getNodeXi
};  // class SBPQuartic




FieldShape* getDG1SBP3Shape(int order)
{
  static DG1SBP3Linear linear1;
  static DG1SBP3Quadratic quadratic1;
  static DG1SBP3Cubic cubic1;
  static DG1SBP3Quartic quartic1;
  switch (order) {
    case 1:
	  return &linear1;
//	case 2:
//	  return &quadratic1;
//	case 3:
//	  return &cubic1;
//	case 4:
//	  return &quartic1;
	default:
	  std::cout << "order " << order << " is not supported by dgSBP3Shape1.cc" << std::endl;
	  return NULL;
    }
}

} // end namespace apf
