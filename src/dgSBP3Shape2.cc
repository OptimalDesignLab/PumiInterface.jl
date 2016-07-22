#include "apfShape.h"
#include "apfMesh.h"
#include "apfVector.h"
#include "apfMatrix.h"
#include <iostream>
#include "dgSBP3Shape2.h"

namespace apf {

// SBP-Gamma for DG
class DG2SBP3Linear : public FieldShape
{
  public:

    DG2SBP3Linear() { registerSelf(apf::DG2SBP3Linear::getName()); }
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "DG2SBP3Linear"; }

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
          { xi = Vector3(0, 0, 0); break; }
        case 1:
          { xi = Vector3(1.0, 0, 0); break; }
        case 2:
          { xi = Vector3(0, 1.0, 0); break; }
        case 3:
          { xi = Vector3(0, 0, 1.0); break; }
        default: 
          { xi = Vector3(0, 0,0); break; }
      }

    }
};  // class SBPLinear



class DG2SBP3Quadratic : public FieldShape
{
  public:
    DG2SBP3Quadratic() { registerSelf(apf::DG2SBP3Quadratic::getName()); }
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "DG2SBP3Quadratic"; }
	
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

      // no need to do this for DG
    
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
      if (dimension == 3)
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if (type == Mesh::TET)
        return 11;
      else
        return 0;
    }

    int getOrder() {return 2;}
    void getNodeXi(int type, int node, Vector3& xi)
    {
	  // return the xi coordinates of the specified node of the specified type
	  // *in the coordinate system of that type*
	  // which makes this function not useful, because the user could define the origins differently
        switch (node)
        {
        case 0:
          { xi = Vector3(0.0, 0.0, 0.0); break; }
        case 1:
          { xi = Vector3(1.0, 0.0, 0.0); break; }
        case 2:
          { xi = Vector3(0.0, 1.0, 0.0); break; }
        case 3:
          { xi = Vector3(0.0, 0.0, 1.0); break; }
        case 4:
          { xi = Vector3(0.5, 0.0, 0.0); break; }
        case 5:
          { xi = Vector3(0.5, 0.5, 0.0); break; }
        case 6:
          { xi = Vector3(0.0, 0.5, 0.0); break; }
        case 7:
          { xi = Vector3(0.0, 0.0, 0.5); break; }
        case 8:
          { xi = Vector3(0.5, 0.0, 0.5); break; }
        case 9:
          { xi = Vector3(0.0, 0.5, 0.5); break; }
        case 10:
          { xi = Vector3(0.25, 0.25, 0.25); break; }
        default: 
          { xi = Vector3(0, 0,0); break; }
        }

    }
};  // class SBPQuadratic




class DG2SBP3Cubic : public FieldShape
{
  public:
    DG2SBP3Cubic() { registerSelf(apf::DG2SBP3Cubic::getName()); }
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "DG2SBP3Cubic"; }
	
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
      if (dimension == 3)
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if (type == Mesh::TET)
        return 24;
      else
        return 0;
    }
    int getOrder() {return 3;}
    void getNodeXi(int type, int node, Vector3& xi)
    {
	  // return the xi coordinates of the specified node of the specified type
	  // *in the coordinate system of that type*
        switch (node)
        {
        case 0:
          { xi = Vector3(0.0, 0.0, 0.0); break; }
        case 1:
          { xi = Vector3(1.0, 0.0, 0.0); break; }
        case 2:
          { xi = Vector3(0.0, 1.0, 0.0); break; }
        case 3:
          { xi = Vector3(0.0, 0.0, 1.0); break; }
        case 4:
          { xi = Vector3(0.3333333333333333, 0.3333333333333333, 0.0); break; }
        case 5:
          { xi = Vector3(0.3333333333333333, 0.0, 0.3333333333333333); break; }
        case 6:
          { xi = Vector3(0.3333333333333333, 0.3333333333333333, 0.3333333333333333); break; }
        case 7:
          { xi = Vector3(0.0, 0.3333333333333333, 0.3333333333333333); break; }
        case 8:
          { xi = Vector3(0.1524029491994481, 0.1524029491994481, 0.5427911524016557); break; }
        case 9:
          { xi = Vector3(0.1524029491994481, 0.5427911524016557, 0.1524029491994481); break; }
        case 10:
          { xi = Vector3(0.1524029491994481, 0.1524029491994481, 0.1524029491994481); break; }
        case 11:
          { xi = Vector3(0.5427911524016557, 0.1524029491994481, 0.1524029491994481); break; }
        case 12:
          { xi = Vector3(0.6951941016011038, 0.0, 0.0); break; }
        case 13:
          { xi = Vector3(0.30480589839889616, 0.0, 0.0); break; }
        case 14:
          { xi = Vector3(0.30480589839889616, 0.6951941016011038, 0.0); break; }
        case 15:
          { xi = Vector3(0.6951941016011038, 0.30480589839889616, 0.0); break; }
        case 16:
          { xi = Vector3(0.0, 0.30480589839889616, 0.0); break; }
        case 17:
          { xi = Vector3(0.0, 0.6951941016011038, 0.0); break; }
        case 18:
          { xi = Vector3(0.0, 0.0, 0.6951941016011038); break; }
        case 19:
          { xi = Vector3(0.0, 0.0, 0.30480589839889616); break; }
        case 20:
          { xi = Vector3(0.30480589839889616, 0.0, 0.6951941016011038); break; }
        case 21:
          { xi = Vector3(0.6951941016011038, 0.0, 0.30480589839889616); break; }
        case 22:
          { xi = Vector3(0.0, 0.30480589839889616, 0.6951941016011038); break; }
        case 23:
          { xi = Vector3(0.0, 0.6951941016011038, 0.30480589839889616); break; }
        default: 
          { xi = Vector3(0, 0,0); break; }
        }
    }
};  // class SBPCubic



class DG2SBP3Quartic : public FieldShape
{
  public:
    DG2SBP3Quartic() { registerSelf(apf::DG2SBP3Quartic::getName()); }
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "DG2SBP3Quartic"; }
	
	
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
      if (dimension == 3)
        return true;
      else
        return false;
    }
    int countNodesOn(int type)
    {
      if (type == Mesh::TET)
        return 45;
      else
        return 0;
    }
    int getOrder() {return 4;}
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
          {xi = Vector3(0.0, 0.0, 1.0); break; }
        case 4:
          {xi = Vector3(0.09472900091823398, 0.09472900091823398, 0.715812997245298); break; }
        case 5:
          {xi = Vector3(0.09472900091823398, 0.715812997245298, 0.09472900091823398); break; }
        case 6:
          {xi = Vector3(0.09472900091823398, 0.09472900091823398, 0.09472900091823398); break; }
        case 7:
          {xi = Vector3(0.715812997245298, 0.09472900091823398, 0.09472900091823398); break; }
        case 8:
          {xi = Vector3(0.5, 0.0, 0.0); break; }
        case 9:
          {xi = Vector3(0.5, 0.5, 0.0); break; }
        case 10:
          {xi = Vector3(0.0, 0.5, 0.0); break; }
        case 11:
          {xi = Vector3(0.0, 0.0, 0.5); break; }
        case 12:
          {xi = Vector3(0.5, 0.0, 0.5); break; }
        case 13:
          {xi = Vector3(0.0, 0.5, 0.5); break; }
        case 14:
          {xi = Vector3(0.10871416009777772, 0.39128583990222227, 0.39128583990222227); break; }
        case 15:
          {xi = Vector3(0.39128583990222227, 0.10871416009777772, 0.39128583990222227); break; }
        case 16:
          {xi = Vector3(0.39128583990222227, 0.39128583990222227, 0.10871416009777772); break; }
        case 17:
          {xi = Vector3(0.10871416009777772, 0.10871416009777772, 0.39128583990222227); break; }
        case 18:
          {xi = Vector3(0.10871416009777772, 0.39128583990222227, 0.10871416009777772); break; }
        case 19:
          {xi = Vector3(0.39128583990222227, 0.10871416009777772, 0.10871416009777772); break; }
        case 20:
          {xi = Vector3(0.7426272531851917, 0.0, 0.0); break; }
        case 21:
          {xi = Vector3(0.25737274681480826, 0.0, 0.0); break; }
        case 22:
          {xi = Vector3(0.25737274681480826, 0.7426272531851917, 0.0); break; }
        case 23:
          {xi = Vector3(0.7426272531851917, 0.25737274681480826, 0.0); break; }
        case 24:
          {xi = Vector3(0.0, 0.25737274681480826, 0.0); break; }
        case 25:
          {xi = Vector3(0.0, 0.7426272531851917, 0.0); break; }
        case 26:
          {xi = Vector3(0.0, 0.0, 0.7426272531851917); break; }
        case 27:
          {xi = Vector3(0.0, 0.0, 0.25737274681480826); break; }
        case 28:
          {xi = Vector3(0.25737274681480826, 0.0, 0.7426272531851917); break; }
        case 29:
          {xi = Vector3(0.7426272531851917, 0.0, 0.25737274681480826); break; }
        case 30:
          {xi = Vector3(0.0, 0.25737274681480826, 0.7426272531851917); break; }
        case 31:
          {xi = Vector3(0.0, 0.7426272531851917, 0.25737274681480826); break; }
        case 32:
          {xi = Vector3(0.22504424155412348, 0.549911516891753, 0.0); break; }
        case 33:
          {xi = Vector3(0.22504424155412348, 0.22504424155412348, 0.0); break; }
        case 34:
          {xi = Vector3(0.549911516891753, 0.22504424155412348, 0.0); break; }
        case 35:
          {xi = Vector3(0.549911516891753, 0.0, 0.22504424155412348); break; }
        case 36:
          {xi = Vector3(0.22504424155412348, 0.0, 0.22504424155412348); break; }
        case 37:
          {xi = Vector3(0.22504424155412348, 0.0, 0.549911516891753); break; }
        case 38:
          {xi = Vector3(0.22504424155412348, 0.549911516891753, 0.22504424155412348); break; }
        case 39:
          {xi = Vector3(0.549911516891753, 0.22504424155412348, 0.22504424155412348); break; }
        case 40:
          {xi = Vector3(0.22504424155412348, 0.22504424155412348, 0.549911516891753); break; }
        case 41:
          {xi = Vector3(0.0, 0.22504424155412348, 0.549911516891753); break; }
        case 42:
          {xi = Vector3(0.0, 0.22504424155412348, 0.22504424155412348); break; }
        case 43:
          {xi = Vector3(0.0, 0.549911516891753, 0.22504424155412348); break; }
        case 44:
          {xi = Vector3(0.25, 0.25, 0.25); break; }
        default:
          {xi = Vector3(0, 0, 0); break; }
      }
  
    }  // end getNodeXi
};  // class SBPQuartic




FieldShape* getDG2SBP3Shape(int order)
{
  static DG2SBP3Linear linear1;
  static DG2SBP3Quadratic quadratic1;
  static DG2SBP3Cubic cubic1;
  static DG2SBP3Quartic quartic1;
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
	  std::cout << "order " << order << " is not supported by dgSBP3Shape1.cc" << std::endl;
	  return NULL;
    }
}

} // end namespace apf
