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
    IsotropicFunctionJ(apf::Mesh2* m_, double(*sizefunc)(apf::MeshEntity*vert, apf::Mesh2* m_ptr, double* u), double* u_): m(m_), juliafunc(sizefunc), u(u_) { } // useful constructor

    // data members
    apf::Mesh2* m;

    double getValue(apf::MeshEntity* vert);

  private:
    double (*juliafunc)(apf::MeshEntity* vert, apf::Mesh2* m_ptr, double *u);
    double *u;  // pointer to array that holds u (avoids copying in constructor)




}; // end of class IsotropicfunctionJ



#endif
