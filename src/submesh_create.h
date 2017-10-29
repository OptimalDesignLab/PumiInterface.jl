// header file for submesh creation

#ifndef SUBMESH_CREATE_H
#define SUBMESH_CREATE_H
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <unordered_set>

// Pumi headers
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <ma.h>
//#include <stdlib.h>   // malloc, free, etc.
//#include <math.h>
//#include <string.h>

extern "C" {
apf::Mesh2* createSubMesh(apf::Mesh* m, int* el_list, int numel);
} // extern C

#endif
