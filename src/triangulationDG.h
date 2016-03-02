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

extern "C" {

/* Create a new mesh where each element of the old mesh is divided into several 
 * triangles.  Only for visualization
 * Inputs:
 *   m: the old Mesh*
 *   numtriangles: the number of triangles each old mesh element is divided into
 *   triangulation: see explanation in transferField
 *   elementNodeOffsets: the offsets used to remap nodes on shared MeshEntityies
 *   typeOffsetsPerElement: the index where the nodes of each dimension MeshEntity
 *                          start in a list of all nodes on the element
 *   numberings[]: the array of Numbering* for verts, edges, and faces
 * Outputs:
 *   m_new: the new Mesh*
*/
apf::Mesh2* createSubMesh(apf::Mesh* m, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3]);


/* Copies a Field from the old mesh to the new mesh
 * Inputs:
 *   m: the old Mesh*
 *   m_new: the new Mesh*
 *   numtriangles: the number of triangles each old mesh element is divided into
 *   triangulation: the numtriangles x 3 array that specifies which vertices and
 *                  nodes of the old mesh are used as vertices of the new mesh
 *                  triangles.  The first 3 integers correspond to the vertices
 *                  of the old mesh element, then the element nodes are labelled.
 *                  1-based indexing
 *  elementNodeOffsets[]: the indices telling where the nodes on each type of 
 *                        MeshEntity start 
 *  numberings[]: the Numbering* of the verts, edges, and faces
 *  interp_op: the 3 x num_nodes_per_el array that interpoltes the solution values
 *             to the vertices of the old mesh element
 *  field_old: the Field* on the old mesh
 *  field_new: the field* on the new mesh (must already have been created)
*/
void transferField(apf::Mesh* m, apf::Mesh* m_new, const int numtriangles, const int triangulation[][3], uint8_t elementNodeOffsets[], int typeOffsetsPerElement[], apf::Numbering* numberings[3], double* interp_op, apf::Field* field_old, apf::Field* field_new);
}

#endif
