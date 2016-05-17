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
#include <stdint.h>
#include "triangulationDG.h"

extern "C" {
apf::Mesh2* createSubMesh(apf::Mesh* m, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3]);


void transferField(apf::Mesh* m, apf::Mesh* m_new, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3], apf::Field* field_old, apf::Field* field_new);

}

#endif
