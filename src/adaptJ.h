// Wrappers for functions to access mesh adaptation functionality from Julia


#include <apf.h>
//#include <gmi_mesh.h>
//#include <gmi_null.h>
//#include <apfMDS.h>
//#include <apfMesh2.h>
//#include <PCU.h>
#include <apfNumbering.h>
//#include <apfShape.h>
#include <ma.h>
//#include <cmath>
#include <cassert>

#ifndef ADAPTJ_H
#define ADAPTJ_H

class IsotropicFunctionJ : public ma::IsotropicFunction
{

  public:

    IsotropicFunctionJ() {};  // default constructor
    IsotropicFunctionJ(apf::Mesh* _m, apf::Field* _f)
    {
      // this works because the FieldShapes are singletons
      assert( _m->getShape() == apf::getShape(_f) );
      m = _m;
      f = _f;
    }

    double getValue(apf::MeshEntity* e)
    {
      // retrieve the computed value from the field
      double val =  apf::getScalar(f, e, 0);
      return val;
    }

  private:

    apf::Mesh* m = NULL;
    // Field (same FieldShape as the coordinate) that stores the desired
    // edge length
    // Adjoint-based error estimators provide an elementwise error estimate,
    // so this field is usually the minimum size of the elements that are
    // upward adjacent of the entity
    apf::Field* f;

};
#endif
