#include "apfShape.h"
#include "apfMesh.h"
#include "apfVector.h"
#include "apfMatrix.h"
#include <iostream>
#include "dgSBP3Shape4.h"

namespace apf {
// 3d Diagonal E operators

class DG4SBP3Linear : public FieldShape
{
  public:

    DG4SBP3Linear() { registerSelf(apf::DG4SBP3Linear::getName()); }
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "DG4SBP3Linear"; }

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
    
        int countNodes() const {return 13;}
    
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
        return 13;
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
          {xi = Vector3(0.16666666666666663, 0.6666666666666667, 0.0); break; }
        case 1:
          {xi = Vector3(0.16666666666666663, 0.16666666666666663, 0.0); break; }
        case 2:
          {xi = Vector3(0.6666666666666667, 0.16666666666666663, 0.0); break; }
        case 3:
          {xi = Vector3(0.6666666666666667, 0.0, 0.16666666666666663); break; }
        case 4:
          {xi = Vector3(0.16666666666666663, 0.0, 0.16666666666666663); break; }
        case 5:
          {xi = Vector3(0.16666666666666663, 0.0, 0.6666666666666667); break; }
        case 6:
          {xi = Vector3(0.16666666666666663, 0.6666666666666667, 0.16666666666666663); break; }
        case 7:
          {xi = Vector3(0.6666666666666667, 0.16666666666666663, 0.16666666666666663); break; }
        case 8:
          {xi = Vector3(0.16666666666666663, 0.16666666666666663, 0.6666666666666667); break; }
        case 9:
          {xi = Vector3(0.0, 0.16666666666666663, 0.6666666666666667); break; }
        case 10:
          {xi = Vector3(0.0, 0.16666666666666663, 0.16666666666666663); break; }
        case 11:
          {xi = Vector3(0.0, 0.6666666666666667, 0.16666666666666663); break; }
        case 12:
          {xi = Vector3(0.25, 0.25, 0.25); break; }
        default:
          {xi = Vector3(0, 0, 0); break; }
      }

    }
};  // class SBPLinear



class DG4SBP3Quadratic : public FieldShape
{
  public:
    DG4SBP3Quadratic() { registerSelf(apf::DG4SBP3Quadratic::getName()); }
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "DG4SBP3Quadratic"; }
	
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
    
        int countNodes() const {return 36;}
    
      void alignSharedNodes(Mesh* m, MeshEntity* elem, MeshEntity* shared, int order[])
      // elem is the triangle 
      // shared is the entity (edge or vertex) being shared
      // order[] contains the mapping such that order[i], where i is the local node number, give
      // the position of that node in the canonical ordering
      {

        // no need to do this for DG
      
      }
    };	// class Tetrahdronf

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
      if (type == Mesh::VERTEX)
	  {
            return 0;
	  } else if ( type == Mesh::EDGE) 
	  {
	    return 0;
	  } else if ( type == Mesh::TRIANGLE)
	  {
	    return 0;
	  }  else if ( type == Mesh::TET)
	  {
	     return 36;
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
        switch (node)
        {
          case 0:
            {xi = Vector3(0.4459484909159649, 0.10810301816807022, 0.0); break; }
          case 1:
            {xi = Vector3(0.4459484909159649, 0.4459484909159649, 0.0); break; }
          case 2:
            {xi = Vector3(0.10810301816807022, 0.4459484909159649, 0.0); break; }
          case 3:
            {xi = Vector3(0.10810301816807022, 0.0, 0.4459484909159649); break; }
          case 4:
            {xi = Vector3(0.4459484909159649, 0.0, 0.4459484909159649); break; }
          case 5:
            {xi = Vector3(0.4459484909159649, 0.0, 0.10810301816807022); break; }
          case 6:
            {xi = Vector3(0.4459484909159649, 0.10810301816807022, 0.4459484909159649); break; }
          case 7:
            {xi = Vector3(0.10810301816807022, 0.4459484909159649, 0.4459484909159649); break; }
          case 8:
            {xi = Vector3(0.4459484909159649, 0.4459484909159649, 0.10810301816807022); break; }
          case 9:
            {xi = Vector3(0.0, 0.4459484909159649, 0.10810301816807022); break; }
          case 10:
            {xi = Vector3(0.0, 0.4459484909159649, 0.4459484909159649); break; }
          case 11:
            {xi = Vector3(0.0, 0.10810301816807022, 0.4459484909159649); break; }
          case 12:
            {xi = Vector3(0.09157621350977074, 0.8168475729804585, 0.0); break; }
          case 13:
            {xi = Vector3(0.09157621350977074, 0.09157621350977074, 0.0); break; }
          case 14:
            {xi = Vector3(0.8168475729804585, 0.09157621350977074, 0.0); break; }
          case 15:
            {xi = Vector3(0.8168475729804585, 0.0, 0.09157621350977074); break; }
          case 16:
            {xi = Vector3(0.09157621350977074, 0.0, 0.09157621350977074); break; }
          case 17:
            {xi = Vector3(0.09157621350977074, 0.0, 0.8168475729804585); break; }
          case 18:
            {xi = Vector3(0.09157621350977074, 0.8168475729804585, 0.09157621350977074); break; }
          case 19:
            {xi = Vector3(0.8168475729804585, 0.09157621350977074, 0.09157621350977074); break; }
          case 20:
            {xi = Vector3(0.09157621350977074, 0.09157621350977074, 0.8168475729804585); break; }
          case 21:
            {xi = Vector3(0.0, 0.09157621350977074, 0.8168475729804585); break; }
          case 22:
            {xi = Vector3(0.0, 0.09157621350977074, 0.09157621350977074); break; }
          case 23:
            {xi = Vector3(0.0, 0.8168475729804585, 0.09157621350977074); break; }
          case 24:
            {xi = Vector3(0.20199409829748438, 0.502459307676599, 0.09355249572843227); break; }
          case 25:
            {xi = Vector3(0.20199409829748438, 0.20199409829748438, 0.09355249572843227); break; }
          case 26:
            {xi = Vector3(0.502459307676599, 0.20199409829748438, 0.09355249572843227); break; }
          case 27:
            {xi = Vector3(0.502459307676599, 0.09355249572843227, 0.20199409829748438); break; }
          case 28:
            {xi = Vector3(0.20199409829748438, 0.09355249572843227, 0.20199409829748438); break; }
          case 29:
            {xi = Vector3(0.20199409829748438, 0.09355249572843227, 0.502459307676599); break; }
          case 30:
            {xi = Vector3(0.20199409829748438, 0.502459307676599, 0.20199409829748438); break; }
          case 31:
            {xi = Vector3(0.502459307676599, 0.20199409829748438, 0.20199409829748438); break; }
          case 32:
            {xi = Vector3(0.20199409829748438, 0.20199409829748438, 0.502459307676599); break; }
          case 33:
            {xi = Vector3(0.09355249572843227, 0.20199409829748438, 0.502459307676599); break; }
          case 34:
            {xi = Vector3(0.09355249572843227, 0.20199409829748438, 0.20199409829748438); break; }
          case 35:
            {xi = Vector3(0.09355249572843227, 0.502459307676599, 0.20199409829748438); break; }
          default:
            {xi = Vector3(0, 0, 0); break; }
        }

    }
};  // class SBPQuadratic

// not implemented
class DG4SBP3Cubic : public FieldShape
{
  public:
    DG4SBP3Cubic() { registerSelf(apf::DG4SBP3Cubic::getName()); }
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "DG4SBP3Cubic"; }
	
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
    
        int countNodes() const {return 69;}
    
        void alignSharedNodes(Mesh* m, MeshEntity* elem, MeshEntity* shared, int order[])
        // elem is the triangle 
        // shared is the entity (edge or vertex) being shared
        // order[] contains the mapping such that order[i], where i is the local node number, give
        // the position of that node in the canonical ordering
        {

          // no need to do this for DG
        
        }
    };  // class Tetrahderon	

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
      if (type == Mesh::VERTEX)
	  {
            return 0;
	  } else if ( type == Mesh::EDGE) 
	  {
	    return 0;
	  } else if ( type == Mesh::TRIANGLE)
	  {
	    return 0;
	  }  else if ( type == Mesh::TET)
	  {
	     return 69;
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
	  // which makes this function not useful, because the user could define the origins differently
        switch (node)
        {
          case 0:
            {xi = Vector3(0.07253084755798583, 0.07253084755798583, 0.7824074573260424); break; }
          case 1:
            {xi = Vector3(0.07253084755798583, 0.7824074573260424, 0.07253084755798589); break; }
          case 2:
            {xi = Vector3(0.07253084755798583, 0.07253084755798583, 0.07253084755798583); break; }
          case 3:
            {xi = Vector3(0.7824074573260424, 0.07253084755798583, 0.07253084755798583); break; }
          case 4:
            {xi = Vector3(0.3011452278959582, 0.3011452278959582, 0.09656431631212525); break; }
          case 5:
            {xi = Vector3(0.3011452278959582, 0.09656431631212525, 0.3011452278959582); break; }
          case 6:
            {xi = Vector3(0.3011452278959582, 0.3011452278959582, 0.3011452278959582); break; }
          case 7:
            {xi = Vector3(0.09656431631212525, 0.3011452278959582, 0.3011452278959582); break; }
          case 8:
            {xi = Vector3(0.06308901449150228, 0.8738219710169954, 0.0); break; }
          case 9:
            {xi = Vector3(0.06308901449150228, 0.06308901449150228, 0.0); break; }
          case 10:
            {xi = Vector3(0.8738219710169954, 0.06308901449150228, 0.0); break; }
          case 11:
            {xi = Vector3(0.8738219710169954, 0.0, 0.06308901449150228); break; }
          case 12:
            {xi = Vector3(0.06308901449150228, 0.0, 0.06308901449150228); break; }
          case 13:
            {xi = Vector3(0.06308901449150228, 0.0, 0.8738219710169954); break; }
          case 14:
            {xi = Vector3(0.06308901449150228, 0.8738219710169954, 0.06308901449150228); break; }
          case 15:
            {xi = Vector3(0.8738219710169954, 0.06308901449150228, 0.06308901449150228); break; }
          case 16:
            {xi = Vector3(0.06308901449150228, 0.06308901449150228, 0.8738219710169954); break; }
          case 17:
            {xi = Vector3(0.0, 0.06308901449150228, 0.8738219710169954); break; }
          case 18:
            {xi = Vector3(0.0, 0.06308901449150228, 0.06308901449150228); break; }
          case 19:
            {xi = Vector3(0.0, 0.8738219710169954, 0.06308901449150228); break; }
          case 20:
            {xi = Vector3(0.24928674517091048, 0.501426509658179, 0.0); break; }
          case 21:
            {xi = Vector3(0.24928674517091048, 0.24928674517091043, 0.0); break; }
          case 22:
            {xi = Vector3(0.501426509658179, 0.24928674517091043, 0.0); break; }
          case 23:
            {xi = Vector3(0.501426509658179, 0.0, 0.24928674517091048); break; }
          case 24:
            {xi = Vector3(0.24928674517091043, 0.0, 0.24928674517091048); break; }
          case 25:
            {xi = Vector3(0.24928674517091043, 0.0, 0.501426509658179); break; }
          case 26:
            {xi = Vector3(0.24928674517091048, 0.501426509658179, 0.24928674517091048); break; }
          case 27:
            {xi = Vector3(0.501426509658179, 0.24928674517091043, 0.24928674517091048); break; }
          case 28:
            {xi = Vector3(0.24928674517091048, 0.24928674517091043, 0.501426509658179); break; }
          case 29:
            {xi = Vector3(0.0, 0.24928674517091048, 0.501426509658179); break; }
          case 30:
            {xi = Vector3(0.0, 0.24928674517091048, 0.24928674517091043); break; }
          case 31:
            {xi = Vector3(0.0, 0.501426509658179, 0.24928674517091043); break; }
          case 32:
            {xi = Vector3(0.08334025586798127, 0.5574175088221697, 0.2759019794418678); break; }
          case 33:
            {xi = Vector3(0.08334025586798127, 0.08334025586798122, 0.2759019794418678); break; }
          case 34:
            {xi = Vector3(0.5574175088221697, 0.08334025586798122, 0.2759019794418678); break; }
          case 35:
            {xi = Vector3(0.5574175088221697, 0.2759019794418678, 0.08334025586798127); break; }
          case 36:
            {xi = Vector3(0.08334025586798122, 0.2759019794418678, 0.08334025586798127); break; }
          case 37:
            {xi = Vector3(0.08334025586798122, 0.2759019794418678, 0.5574175088221697); break; }
          case 38:
            {xi = Vector3(0.08334025586798127, 0.5574175088221697, 0.08334025586798127); break; }
          case 39:
            {xi = Vector3(0.5574175088221697, 0.08334025586798122, 0.08334025586798127); break; }
          case 40:
            {xi = Vector3(0.08334025586798127, 0.08334025586798122, 0.5574175088221697); break; }
          case 41:
            {xi = Vector3(0.2759019794418678, 0.08334025586798127, 0.5574175088221697); break; }
          case 42:
            {xi = Vector3(0.2759019794418678, 0.08334025586798127, 0.08334025586798122); break; }
          case 43:
            {xi = Vector3(0.2759019794418678, 0.5574175088221697, 0.08334025586798122); break; }
          case 44:
            {xi = Vector3(0.3103524510337843, 0.6365024991213988, 0.0); break; }
          case 45:
            {xi = Vector3(0.05314504984481694, 0.6365024991213988, 0.0); break; }
          case 46:
            {xi = Vector3(0.05314504984481694, 0.3103524510337843, 0.0); break; }
          case 47:
            {xi = Vector3(0.3103524510337843, 0.05314504984481694, 0.0); break; }
          case 48:
            {xi = Vector3(0.6365024991213988, 0.05314504984481694, 0.0); break; }
          case 49:
            {xi = Vector3(0.6365024991213988, 0.3103524510337843, 0.0); break; }
          case 50:
            {xi = Vector3(0.6365024991213988, 0.0, 0.3103524510337843); break; }
          case 51:
            {xi = Vector3(0.6365024991213988, 0.0, 0.05314504984481694); break; }
          case 52:
            {xi = Vector3(0.3103524510337843, 0.0, 0.05314504984481694); break; }
          case 53:
            {xi = Vector3(0.05314504984481694, 0.0, 0.3103524510337843); break; }
          case 54:
            {xi = Vector3(0.05314504984481694, 0.0, 0.6365024991213988); break; }
          case 55:
            {xi = Vector3(0.3103524510337843, 0.0, 0.6365024991213988); break; }
          case 56:
            {xi = Vector3(0.05314504984481694, 0.6365024991213988, 0.3103524510337843); break; }
          case 57:
            {xi = Vector3(0.3103524510337843, 0.6365024991213988, 0.05314504984481694); break; }
          case 58:
            {xi = Vector3(0.6365024991213988, 0.3103524510337843, 0.05314504984481694); break; }
          case 59:
            {xi = Vector3(0.6365024991213988, 0.05314504984481694, 0.3103524510337843); break; }
          case 60:
            {xi = Vector3(0.3103524510337843, 0.05314504984481694, 0.6365024991213988); break; }
          case 61:
            {xi = Vector3(0.05314504984481694, 0.3103524510337843, 0.6365024991213988); break; }
          case 62:
            {xi = Vector3(0.0, 0.3103524510337843, 0.6365024991213988); break; }
          case 63:
            {xi = Vector3(0.0, 0.05314504984481694, 0.6365024991213988); break; }
          case 64:
            {xi = Vector3(0.0, 0.05314504984481694, 0.3103524510337843); break; }
          case 65:
            {xi = Vector3(0.0, 0.3103524510337843, 0.05314504984481694); break; }
          case 66:
            {xi = Vector3(0.0, 0.6365024991213988, 0.05314504984481694); break; }
          case 67:
            {xi = Vector3(0.0, 0.6365024991213988, 0.3103524510337843); break; }
          case 68:
            {xi = Vector3(0.25, 0.25, 0.25); break; }
          default:
            {xi = Vector3(0, 0, 0); break; }
        }  // end case statement
    }  // end function getNodeXi
};  // class SBPCubic

// not implemented
class DG4SBP3Quartic : public FieldShape
{
  public:
    DG4SBP3Quartic() { registerSelf(apf::DG4SBP3Quartic::getName()); }
//    SBPLinear() { registerSelf(apf::Linear::getName()); }  // use inherited/default constructor?
    const char* getName() const { return "DG4SBP3Quartic"; }
	
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
    
        int countNodes() const {return 99;}
    
        void alignSharedNodes(Mesh* m, MeshEntity* elem, MeshEntity* shared, int order[])
        // elem is the triangle 
        // shared is the entity (edge or vertex) being shared
        // order[] contains the mapping such that order[i], where i is the local node number, give
        // the position of that node in the canonical ordering
        {

          // no need to do this for DG
        
        }
    };  // class Tetrahderon	

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
      if (type == Mesh::VERTEX)
	  {
            return 0;
	  } else if ( type == Mesh::EDGE) 
	  {
	    return 0;
	  } else if ( type == Mesh::TRIANGLE)
	  {
	    return 0;
	  }  else if ( type == Mesh::TET)
	  {
	     return 99;
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
	  // which makes this function not useful, because the user could define the origins differently
        switch (node)
        {

          case 0:
            {xi = Vector3(0.33333333333333337, 0.33333333333333337, 0.0); break; }
          case 1:
            {xi = Vector3(0.33333333333333337, 0.0, 0.33333333333333337); break; }
          case 2:
            {xi = Vector3(0.33333333333333337, 0.33333333333333337, 0.33333333333333337); break; }
          case 3:
            {xi = Vector3(0.0, 0.33333333333333337, 0.33333333333333337); break; }
          case 4:
            {xi = Vector3(0.0612198866702342, 0.0612198866702342, 0.8163403399892974); break; }
          case 5:
            {xi = Vector3(0.0612198866702342, 0.8163403399892974, 0.0612198866702342); break; }
          case 6:
            {xi = Vector3(0.0612198866702342, 0.0612198866702342, 0.0612198866702342); break; }
          case 7:
            {xi = Vector3(0.8163403399892974, 0.0612198866702342, 0.0612198866702342); break; }
          case 8:
            {xi = Vector3(0.41379484327386623, 0.08620515672613382, 0.08620515672613371); break; }
          case 9:
            {xi = Vector3(0.08620515672613382, 0.41379484327386623, 0.08620515672613371); break; }
          case 10:
            {xi = Vector3(0.08620515672613377, 0.08620515672613377, 0.4137948432738663); break; }
          case 11:
            {xi = Vector3(0.41379484327386623, 0.41379484327386623, 0.08620515672613371); break; }
          case 12:
            {xi = Vector3(0.41379484327386623, 0.08620515672613377, 0.4137948432738663); break; }
          case 13:
            {xi = Vector3(0.08620515672613377, 0.41379484327386623, 0.4137948432738663); break; }
          case 14:
            {xi = Vector3(0.1705693077517602, 0.6588613844964796, 0.0); break; }
          case 15:
            {xi = Vector3(0.1705693077517602, 0.17056930775176027, 0.0); break; }
          case 16:
            {xi = Vector3(0.6588613844964796, 0.17056930775176027, 0.0); break; }
          case 17:
            {xi = Vector3(0.6588613844964796, 0.0, 0.1705693077517602); break; }
          case 18:
            {xi = Vector3(0.17056930775176027, 0.0, 0.1705693077517602); break; }
          case 19:
            {xi = Vector3(0.17056930775176027, 0.0, 0.6588613844964796); break; }
          case 20:
            {xi = Vector3(0.1705693077517602, 0.6588613844964796, 0.1705693077517602); break; }
          case 21:
            {xi = Vector3(0.6588613844964796, 0.17056930775176027, 0.1705693077517602); break; }
          case 22:
            {xi = Vector3(0.1705693077517602, 0.17056930775176027, 0.6588613844964796); break; }
          case 23:
            {xi = Vector3(0.0, 0.1705693077517602, 0.6588613844964796); break; }
          case 24:
            {xi = Vector3(0.0, 0.1705693077517602, 0.17056930775176027); break; }
          case 25:
            {xi = Vector3(0.0, 0.6588613844964796, 0.17056930775176027); break; }
          case 26:
            {xi = Vector3(0.4592925882927231, 0.08141482341455375, 0.0); break; }
          case 27:
            {xi = Vector3(0.4592925882927231, 0.4592925882927231, 0.0); break; }
          case 28:
            {xi = Vector3(0.08141482341455375, 0.4592925882927231, 0.0); break; }
          case 29:
            {xi = Vector3(0.08141482341455375, 0.0, 0.4592925882927231); break; }
          case 30:
            {xi = Vector3(0.4592925882927231, 0.0, 0.4592925882927231); break; }
          case 31:
            {xi = Vector3(0.4592925882927231, 0.0, 0.08141482341455375); break; }
          case 32:
            {xi = Vector3(0.4592925882927231, 0.08141482341455375, 0.4592925882927231); break; }
          case 33:
            {xi = Vector3(0.08141482341455375, 0.4592925882927231, 0.4592925882927231); break; }
          case 34:
            {xi = Vector3(0.4592925882927231, 0.4592925882927231, 0.08141482341455375); break; }
          case 35:
            {xi = Vector3(0.0, 0.4592925882927231, 0.08141482341455375); break; }
          case 36:
            {xi = Vector3(0.0, 0.4592925882927231, 0.4592925882927231); break; }
          case 37:
            {xi = Vector3(0.0, 0.08141482341455375, 0.4592925882927231); break; }
          case 38:
            {xi = Vector3(0.05054722831703096, 0.8989055433659381, 0.0); break; }
          case 39:
            {xi = Vector3(0.05054722831703096, 0.05054722831703096, 0.0); break; }
          case 40:
            {xi = Vector3(0.8989055433659381, 0.05054722831703096, 0.0); break; }
          case 41:
            {xi = Vector3(0.8989055433659381, 0.0, 0.05054722831703096); break; }
          case 42:
            {xi = Vector3(0.05054722831703096, 0.0, 0.05054722831703096); break; }
          case 43:
            {xi = Vector3(0.05054722831703096, 0.0, 0.8989055433659381); break; }
          case 44:
            {xi = Vector3(0.05054722831703096, 0.8989055433659381, 0.05054722831703096); break; }
          case 45:
            {xi = Vector3(0.8989055433659381, 0.05054722831703096, 0.05054722831703096); break; }
          case 46:
            {xi = Vector3(0.05054722831703096, 0.05054722831703096, 0.8989055433659381); break; }
          case 47:
            {xi = Vector3(0.0, 0.05054722831703096, 0.8989055433659381); break; }
          case 48:
            {xi = Vector3(0.0, 0.05054722831703096, 0.05054722831703096); break; }
          case 49:
            {xi = Vector3(0.0, 0.8989055433659381, 0.05054722831703096); break; }
          case 50:
            {xi = Vector3(0.05444165199353779, 0.6867639094442219, 0.20435278656870248); break; }
          case 51:
            {xi = Vector3(0.05444165199353779, 0.05444165199353779, 0.20435278656870248); break; }
          case 52:
            {xi = Vector3(0.6867639094442219, 0.05444165199353779, 0.20435278656870248); break; }
          case 53:
            {xi = Vector3(0.6867639094442219, 0.20435278656870248, 0.05444165199353779); break; }
          case 54:
            {xi = Vector3(0.05444165199353779, 0.20435278656870248, 0.05444165199353779); break; }
          case 55:
            {xi = Vector3(0.05444165199353779, 0.20435278656870248, 0.6867639094442219); break; }
          case 56:
            {xi = Vector3(0.05444165199353779, 0.6867639094442219, 0.05444165199353779); break; }
          case 57:
            {xi = Vector3(0.6867639094442219, 0.05444165199353779, 0.05444165199353779); break; }
          case 58:
            {xi = Vector3(0.05444165199353779, 0.05444165199353779, 0.6867639094442219); break; }
          case 59:
            {xi = Vector3(0.20435278656870248, 0.05444165199353779, 0.6867639094442219); break; }
          case 60:
            {xi = Vector3(0.20435278656870248, 0.05444165199353779, 0.05444165199353779); break; }
          case 61:
            {xi = Vector3(0.20435278656870248, 0.6867639094442219, 0.05444165199353779); break; }
          case 62:
            {xi = Vector3(0.22605711552510166, 0.4769172533468544, 0.07096851560294226); break; }
          case 63:
            {xi = Vector3(0.22605711552510166, 0.22605711552510166, 0.07096851560294226); break; }
          case 64:
            {xi = Vector3(0.4769172533468544, 0.22605711552510166, 0.07096851560294226); break; }
          case 65:
            {xi = Vector3(0.4769172533468544, 0.07096851560294226, 0.22605711552510166); break; }
          case 66:
            {xi = Vector3(0.22605711552510166, 0.07096851560294226, 0.22605711552510166); break; }
          case 67:
            {xi = Vector3(0.22605711552510166, 0.07096851560294226, 0.4769172533468544); break; }
          case 68:
            {xi = Vector3(0.22605711552510166, 0.4769172533468544, 0.22605711552510166); break; }
          case 69:
            {xi = Vector3(0.4769172533468544, 0.22605711552510166, 0.22605711552510166); break; }
          case 70:
            {xi = Vector3(0.22605711552510166, 0.22605711552510166, 0.4769172533468544); break; }
          case 71:
            {xi = Vector3(0.07096851560294226, 0.22605711552510166, 0.4769172533468544); break; }
          case 72:
            {xi = Vector3(0.07096851560294226, 0.22605711552510166, 0.22605711552510166); break; }
          case 73:
            {xi = Vector3(0.07096851560294226, 0.4769172533468544, 0.22605711552510166); break; }
          case 74:
            {xi = Vector3(0.008394777409957643, 0.7284923929554042, 0.0); break; }
          case 75:
            {xi = Vector3(0.2631128296346381, 0.7284923929554042, 0.0); break; }
          case 76:
            {xi = Vector3(0.2631128296346381, 0.008394777409957643, 0.0); break; }
          case 77:
            {xi = Vector3(0.008394777409957643, 0.2631128296346381, 0.0); break; }
          case 78:
            {xi = Vector3(0.7284923929554042, 0.2631128296346381, 0.0); break; }
          case 79:
            {xi = Vector3(0.7284923929554042, 0.008394777409957643, 0.0); break; }
          case 80:
            {xi = Vector3(0.7284923929554042, 0.0, 0.008394777409957643); break; }
          case 81:
            {xi = Vector3(0.7284923929554042, 0.0, 0.2631128296346381); break; }
          case 82:
            {xi = Vector3(0.008394777409957643, 0.0, 0.2631128296346381); break; }
          case 83:
            {xi = Vector3(0.2631128296346381, 0.0, 0.008394777409957643); break; }
          case 84:
            {xi = Vector3(0.2631128296346381, 0.0, 0.7284923929554042); break; }
          case 85:
            {xi = Vector3(0.008394777409957643, 0.0, 0.7284923929554042); break; }
          case 86:
            {xi = Vector3(0.2631128296346381, 0.7284923929554042, 0.008394777409957643); break; }
          case 87:
            {xi = Vector3(0.008394777409957643, 0.7284923929554042, 0.2631128296346381); break; }
          case 88:
            {xi = Vector3(0.7284923929554042, 0.008394777409957643, 0.2631128296346381); break; }
          case 89:
            {xi = Vector3(0.7284923929554042, 0.2631128296346381, 0.008394777409957643); break; }
          case 90:
            {xi = Vector3(0.008394777409957643, 0.2631128296346381, 0.7284923929554042); break; }
          case 91:
            {xi = Vector3(0.2631128296346381, 0.008394777409957643, 0.7284923929554042); break; }
          case 92:
            {xi = Vector3(0.0, 0.008394777409957643, 0.7284923929554042); break; }
          case 93:
            {xi = Vector3(0.0, 0.2631128296346381, 0.7284923929554042); break; }
          case 94:
            {xi = Vector3(0.0, 0.2631128296346381, 0.008394777409957643); break; }
          case 95:
            {xi = Vector3(0.0, 0.008394777409957643, 0.2631128296346381); break; }
          case 96:
            {xi = Vector3(0.0, 0.7284923929554042, 0.2631128296346381); break; }
          case 97:
            {xi = Vector3(0.0, 0.7284923929554042, 0.008394777409957643); break; }
          case 98:
            {xi = Vector3(0.25, 0.25, 0.25); break; }
          default:
            {xi = Vector3(0, 0, 0); break; }
        }  // end case statement

    }  // end function getNodeXi
};  // class SBPQuartic





FieldShape* getDG4SBP3Shape(int order)
{
  static DG4SBP3Linear linear1;
  static DG4SBP3Quadratic quadratic1;
  static DG4SBP3Cubic cubic1;
  static DG4SBP3Quartic quartic1;
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
      std::cout << "order " << order << " is not supported by DG4SBP3Shape1.cc" << std::endl;
      return NULL;
  }
}

} // end namespace apf
