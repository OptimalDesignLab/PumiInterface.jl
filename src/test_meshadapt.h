#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <ma.h>
#include <cmath>

#ifndef TEST_MESHADAPT_H
#define TEST_MESHADAPT_H

class IsotropicFunctionJ : public ma::IsotropicFunction
{

  public:

//    IsotropicFunctionJ() {};  // default constructor
    IsotropicFunctionJ(apf::Mesh* _m, apf::Numbering* _numberings[4],
                       apf::Field* _f)
    {
      m = _m;
      f = _f;

      for (int i=0; i < 4; ++i)
        numberings[i] = _numberings[i];
    }

    double getValue(apf::MeshEntity* vert)
    {

//      int type = m->getType(vert);
//      int dim = m->typeDimension[type];
//      std::cout << "getting size for entity of type " << type << " of dimension " << dim << std::endl;
//
      // retrieve the computed value from the field
      double val =  apf::getScalar(f, vert, 0);
/*
      // compute the exact value
      apf::Vector3 coords;
      m->getPoint(vert, 0, coords);
      // go from element size 2/10 at inner radius to 1/10 at outer radius
      const double h_inner = 2.0/10.0;
      const double h_outer = 1.0/30.0;
      const double r_slope = (h_outer - h_inner)/2.0;
      const double r_in = 1;
      double r = std::sqrt(coords.x()*coords.x() + coords.y()*coords.y());
      double val2 = 2.0/10.0 + r_slope*(r - r_in);
*/ 
//      std::cout << "val  = " << val << std::endl;
//      std::cout << "val2 = " << val2 << std::endl;
      return val;
//      return 2.0/10.0;
    }

  private:

    apf::Mesh* m = NULL;
    apf::Numbering* numberings[4] = {NULL, NULL, NULL, NULL};
    apf::Field* f = NULL;

};
#endif
