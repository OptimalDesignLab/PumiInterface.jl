#ifndef TRIANGULATION_H
#define TRIANGULATION_H


#include <iostream>
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <ma.h>
#include <crv.h>  // curved mesh stuff
#include <stdlib.h>   // malloc, free, etc.
#include <math.h>
#include <string.h>



apf::Mesh2* createSubMesh(apf::Mesh* m, int numtriangles, int triangulation[][3], apf::Numbering* numberings[3]);

#endif
