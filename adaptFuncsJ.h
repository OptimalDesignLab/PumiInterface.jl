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

class AnisotropicFunctionJ : public ma::AnisotropicFunction
{
  public: 

    AnisotropicFunctionJ() { };  // default constructor
    AnisotropicFunctionJ(apf::Mesh2* m_, void(*sizefunc)(apf::MeshEntity*vert,double r[3][3], double h[3], apf::Mesh2* m_ptr, double* u), double* u_): m(m_), juliafunc(sizefunc), u(u_) { } // useful constructor

    // data members
    apf::Mesh2* m;

    void getValue(apf::MeshEntity* vert, ma::Matrix &r, ma::Vector &v);

  private:
    void (*juliafunc)(apf::MeshEntity* vert, double r[3][3], double h[3], apf::Mesh2* m_ptr, double *u);
    double *u;  // pointer to array that holds u (avoids copying in constructor)
    double r_array[3][3];  // array used to turn ma::Matrix into standard data type
    double h_array[3];  // array used to turn vector into standard data type





}; // end of class IsotropicfunctionJ




#endif
