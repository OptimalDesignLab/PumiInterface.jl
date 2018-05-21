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
    
        int countNodes() const {return 20;}
    
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
	     return 20;
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
            {xi = Vector3(0.30973030960236814, 0.30973030960236814, 0.0708090711928957); break; }
          case 1:
            {xi = Vector3(0.30973030960236814, 0.0708090711928957, 0.30973030960236814); break; }
          case 2:
            {xi = Vector3(0.30973030960236814, 0.30973030960236814, 0.30973030960236814); break; }
          case 3:
            {xi = Vector3(0.0708090711928957, 0.30973030960236814, 0.30973030960236814); break; }
          case 4:
            {xi = Vector3(0.06246817700935553, 0.06246817700935553, 0.8125954689719335); break; }
          case 5:
            {xi = Vector3(0.06246817700935553, 0.8125954689719335, 0.06246817700935553); break; }
          case 6:
            {xi = Vector3(0.06246817700935553, 0.06246817700935553, 0.06246817700935553); break; }
          case 7:
            {xi = Vector3(0.8125954689719335, 0.06246817700935553, 0.06246817700935553); break; }
          case 8:
            {xi = Vector3(0.06196811528023993, 0.5934702098815567, 0.2825935595579634); break; }
          case 9:
            {xi = Vector3(0.06196811528023993, 0.06196811528023993, 0.2825935595579634); break; }
          case 10:
            {xi = Vector3(0.5934702098815567, 0.06196811528023993, 0.2825935595579634); break; }
          case 11:
            {xi = Vector3(0.5934702098815567, 0.2825935595579634, 0.06196811528023993); break; }
          case 12:
            {xi = Vector3(0.06196811528023993, 0.2825935595579634, 0.06196811528023993); break; }
          case 13:
            {xi = Vector3(0.06196811528023993, 0.2825935595579634, 0.5934702098815567); break; }
          case 14:
            {xi = Vector3(0.06196811528023993, 0.5934702098815567, 0.06196811528023993); break; }
          case 15:
            {xi = Vector3(0.5934702098815567, 0.06196811528023993, 0.06196811528023993); break; }
          case 16:
            {xi = Vector3(0.06196811528023993, 0.06196811528023993, 0.5934702098815567); break; }
          case 17:
            {xi = Vector3(0.2825935595579634, 0.06196811528023993, 0.5934702098815567); break; }
          case 18:
            {xi = Vector3(0.2825935595579634, 0.06196811528023993, 0.06196811528023993); break; }
          case 19:
            {xi = Vector3(0.2825935595579634, 0.5934702098815567, 0.06196811528023993); break; }
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
    
        int countNodes() const {return 38;}
    
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
	     return 38;
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
            {xi = Vector3(0.060559038064470994, 0.060559038064470994, 0.818322885806587); break; }
          case 1:
            {xi = Vector3(0.060559038064470994, 0.818322885806587, 0.06055903806447105); break; }
          case 2:
            {xi = Vector3(0.060559038064470994, 0.060559038064470994, 0.060559038064470994); break; }
          case 3:
            {xi = Vector3(0.818322885806587, 0.060559038064470994, 0.060559038064470994); break; }
          case 4:
            {xi = Vector3(0.17995491327350105, 0.17995491327350105, 0.4601352601794968); break; }
          case 5:
            {xi = Vector3(0.17995491327350105, 0.4601352601794968, 0.17995491327350105); break; }
          case 6:
            {xi = Vector3(0.17995491327350105, 0.17995491327350105, 0.17995491327350105); break; }
          case 7:
            {xi = Vector3(0.4601352601794968, 0.17995491327350105, 0.17995491327350105); break; }
          case 8:
            {xi = Vector3(0.3585270272483152, 0.1414729727516848, 0.1414729727516848); break; }
          case 9:
            {xi = Vector3(0.1414729727516848, 0.3585270272483152, 0.1414729727516848); break; }
          case 10:
            {xi = Vector3(0.1414729727516848, 0.1414729727516848, 0.3585270272483152); break; }
          case 11:
            {xi = Vector3(0.3585270272483152, 0.3585270272483152, 0.1414729727516848); break; }
          case 12:
            {xi = Vector3(0.3585270272483152, 0.1414729727516848, 0.3585270272483152); break; }
          case 13:
            {xi = Vector3(0.1414729727516848, 0.3585270272483152, 0.3585270272483152); break; }
          case 14:
            {xi = Vector3(0.04406618399879214, 0.6122547631423183, 0.2996128688600974); break; }
          case 15:
            {xi = Vector3(0.04406618399879214, 0.04406618399879214, 0.2996128688600974); break; }
          case 16:
            {xi = Vector3(0.6122547631423183, 0.04406618399879214, 0.2996128688600974); break; }
          case 17:
            {xi = Vector3(0.6122547631423183, 0.2996128688600974, 0.04406618399879214); break; }
          case 18:
            {xi = Vector3(0.04406618399879214, 0.2996128688600974, 0.04406618399879214); break; }
          case 19:
            {xi = Vector3(0.04406618399879214, 0.2996128688600974, 0.6122547631423183); break; }
          case 20:
            {xi = Vector3(0.04406618399879214, 0.6122547631423183, 0.04406618399879214); break; }
          case 21:
            {xi = Vector3(0.6122547631423183, 0.04406618399879214, 0.04406618399879214); break; }
          case 22:
            {xi = Vector3(0.04406618399879214, 0.04406618399879214, 0.6122547631423183); break; }
          case 23:
            {xi = Vector3(0.2996128688600974, 0.04406618399879214, 0.6122547631423183); break; }
          case 24:
            {xi = Vector3(0.2996128688600974, 0.04406618399879214, 0.04406618399879214); break; }
          case 25:
            {xi = Vector3(0.2996128688600974, 0.6122547631423183, 0.04406618399879214); break; }
          case 26:
            {xi = Vector3(0.23441922019715833, 0.026246508568518623, 0.5049150510371647); break; }
          case 27:
            {xi = Vector3(0.23441922019715833, 0.23441922019715833, 0.5049150510371647); break; }
          case 28:
            {xi = Vector3(0.026246508568518623, 0.23441922019715833, 0.5049150510371647); break; }
          case 29:
            {xi = Vector3(0.026246508568518623, 0.5049150510371647, 0.23441922019715833); break; }
          case 30:
            {xi = Vector3(0.23441922019715833, 0.5049150510371647, 0.23441922019715833); break; }
          case 31:
            {xi = Vector3(0.23441922019715833, 0.5049150510371647, 0.026246508568518623); break; }
          case 32:
            {xi = Vector3(0.23441922019715833, 0.026246508568518623, 0.23441922019715833); break; }
          case 33:
            {xi = Vector3(0.026246508568518623, 0.23441922019715833, 0.23441922019715833); break; }
          case 34:
            {xi = Vector3(0.23441922019715833, 0.23441922019715833, 0.026246508568518623); break; }
          case 35:
            {xi = Vector3(0.5049150510371647, 0.23441922019715833, 0.026246508568518623); break; }
          case 36:
            {xi = Vector3(0.5049150510371647, 0.23441922019715833, 0.23441922019715833); break; }
          case 37:
            {xi = Vector3(0.5049150510371647, 0.026246508568518623, 0.23441922019715833); break; }
          default:
            {xi = Vector3(0, 0, 0); break; }
        }  // end case statement

    }  // end function getNodeXi
};  // class SBPQuartic





FieldShape* getDG4SBP3Shape(int order)
{
  static DG4SBP3Linear linear1;
  static DG4SBP3Quadratic quadratic1;
//  static DG4SBP3Cubic cubic1;
//  static DG4SBP3Quartic quartic1;
  switch (order) {
    case 1:
	  return &linear1;
	case 2:
	  return &quadratic1;
//	case 3:
//	  return &cubic1;
//	case 4:
//	  return &quartic1;
	default:
	  std::cout << "order " << order << " is not supported by DG4SBP3Shape1.cc" << std::endl;
	  return NULL;
    }
}

} // end namespace apf
