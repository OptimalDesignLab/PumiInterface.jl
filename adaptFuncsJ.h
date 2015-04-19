#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <ma.h>

// header guard
#ifndef ISOTROPICFUNC_h
#define ISOTROPICFUNC_h

// this header file defines classes derived from PUMI's ma::IsotropicFunction and ma::Anisotropic function to be used for mesh adapation from Julia
//
// dealing with julia callbacks will be interesting

class IsotropicFunctionJ : public ma::IsotropicFunction
{
  public: 

    IsotropicFunctionJ() { };  // default constructor
    IsotropicFunctionJ(apf::Mesh* m_): m(m_) { } // useful constructor

    apf::Mesh* m;

    double getValue(apf::MeshEntity* vert);



}; // end of class IsotropicfunctionJ



#endif
