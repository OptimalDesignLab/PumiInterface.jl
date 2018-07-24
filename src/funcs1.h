#ifndef FUNCS1_H
#define FUNCS1_H

#include <stdint.h>
#include <cstdlib>
// header file for funcs1
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

//#include "adaptFuncsJ.h"
#include "adaptJ.h"
#include "apfSBPShape.h"
#include "apfSBPShape3.h"
#include "dgSBPShape1.h"
#include "dgSBPShape2.h"
#include "dgSBPShape4.h"
#include "dgSBPShape5.h"
#include "dgSBPShape6.h"
#include "dgSBP3Shape1.h"
#include "dgSBP3Shape2.h"
#include "dgSBP3Shape4.h"
#include "dgSBP3Shape5.h"
#include "triangulation.h"
#include "triangulationDG.h"
#include "submesh_create.h"

/*
 * Params[in]:
 * pointer to mesh m
 */
extern "C" {
// this function does not pass opaque pointers
//extern int initABC(char* dmg_name, char* smb_name, int downward_counts[4][4], int numberEntities[4], apf::Mesh2* m_ptr_array[1], apf::FieldShape* mshape_ptr_array[1]);


int initABC(char* dmg_name, char* smb_name, int number_entities[4], apf::Mesh2* m_ptr_array[1], apf::FieldShape* mshape_ptr_array[1], apf::Numbering* n_arr[], int order, int load_mesh, int shape_type );

apf::Mesh2* loadMesh(const char* dmg_name, const char* smb_name, int shape_type,
                     int order, int dim_ret[1]);
void initMesh(apf::Mesh* m, int number_entities[],
       apf::FieldShape* mshape_ptr_array[1], apf::Numbering* n_array[]);


// these functions are not user accessible
void cleanup(apf::Mesh* m_local);

void pushMeshRef(apf::Mesh* m);
void popMeshRef(apf::Mesh* m);

apf::FieldShape* getFieldShape(int shape_type, int order, int dim, bool& change_shape);

// these functions do pass pointers
extern apf::FieldShape* getMeshShapePtr(apf::Mesh* m);
extern apf::FieldShape* getConstantShapePtr(int dimension);

extern int count(apf::Mesh2* m_local, int dimension);

apf::MeshIterator* begin(apf::Mesh* m, int dim);
void end(apf::Mesh* m, apf::MeshIterator* it);
apf::MeshEntity* iterate(apf::Mesh* m, apf::MeshIterator* it);
void iteraten(apf::Mesh* m, apf::MeshIterator* it, int n);
apf::MeshEntity* deref(apf::Mesh2* m, apf::MeshIterator* it);

extern void writeVtkFiles(char* name, apf::Mesh2* m_local);

// geometric model functions

apf::ModelEntity* toModel(apf::Mesh* m_local, apf::MeshEntity* e);
int getModelType(apf::Mesh* m_local, apf::ModelEntity* e);
int getModelTag(apf::Mesh* m_local, apf::ModelEntity* e);

// these functions pass pointers
extern int getMeshDimension(apf::Mesh2* m_local);
extern int getType(apf::Mesh2* m_local, apf::MeshEntity* e);
extern int getDownward(apf::Mesh2* m_local, apf::MeshEntity* e, int dimension, apf::MeshEntity* downwards[12]);
extern int countAdjacent(apf::Mesh2* m_local, apf::MeshEntity* e, int dimension);
void getAdjacent(apf::MeshEntity* adjacencies_ret[]);


int countBridgeAdjacent(apf::Mesh* m_local, apf::MeshEntity* origin, int bridge_dimension, int target_dimension);
void getBridgeAdjacent(apf::MeshEntity* adjacencies_ret2[]);

void getAlignment(apf::Mesh* m_local, apf::MeshEntity* elem, apf::MeshEntity* boundary, int which[1], bool flip[1], int rotate[1]);

// FieldShape functions
bool hasNodesIn(apf::FieldShape* fshape_local, int dimension);
extern int countNodesOn(apf::FieldShape* mshape_ptr, int type);
extern apf::EntityShape* getEntityShape(apf::FieldShape* mshape_local, int type);
int getOrder(apf::FieldShape* fshape);

// MeshElement related functions
extern apf::MeshElement* createMeshElement( apf::Mesh2* m_local, apf::MeshEntity* e);
extern int countIntPoints(apf::MeshElement* e, int order);
extern void getIntPoint(apf::MeshElement* e, int order, int point, double coords[3]);
extern double getIntWeight(apf::MeshElement* e, int order, int point);
extern void getJacobian(apf::MeshElement* e, double coords[3], double jac[3][3]);

// Entityshape functions
extern int countNodes(apf::EntityShape* eshape_local);
extern void getValues(apf::Mesh* m, apf::EntityShape* eshape_local, double xi[3], double vals[]);
extern void getLocalGradients(apf::Mesh* m, apf::EntityShape* eshape_local, double xi[3], double vals[][3]);
void alignSharedNodes(apf::EntityShape* eshape_local, apf::Mesh* m_local, apf::MeshEntity* elem, apf::MeshEntity* shared, int order[]);


// these functinos pass pointers
extern void getVertCoords2(apf::Mesh* m, apf::MeshEntity* e, double coords[][3], int sx, int sy);
extern int getEdgeCoords2(apf::Mesh* m, apf::MeshEntity* e, double coords[2][3], int sx, int sy);
int getFaceCoords2(apf::Mesh* m, apf::MeshEntity* e, double coords[][3], int sx, int sy);
extern int getElCoords2(apf::Mesh* m, apf::MeshEntity* e, double coords[][3], int sx, int sy);

void getAllEntityCoords(apf::Mesh* m, apf::MeshEntity* e, double* coords);

// these function pass pointers
// create a generally defined numbering from julia
extern  apf::Numbering* createNumberingJ(apf::Mesh2* m_local, char* name, apf::FieldShape* field, int components);
extern void destroyNumbering(apf::Numbering* n);
apf::Numbering* findNumbering(apf::Mesh* m, const char* name);
extern void destroyNumberings(apf::Mesh* m, apf::Numbering* save_n[], int n_save);

apf::FieldShape* getNumberingShape(apf::Numbering* n);
extern int numberJ(apf::Numbering* n, apf::MeshEntity* e, int node, int component, int number);
extern  int getNumberJ(apf::Numbering* n, apf::MeshEntity* e, int node, int component);
extern bool isNumbered(apf::Numbering*n, apf::MeshEntity* e, int node, int component);

extern int getDofNumbers(apf::Numbering* n, apf::MeshEntity* entities[], uint8_t node_offsets[], uint8_t nodemap[], apf::MeshEntity* element, int dofnums[]);
void setNumberingOffset(apf::Numbering* num, int off);

void getElementNumbers(apf::Numbering* n, apf::MeshEntity*e, int num_nodes, int nums[]);
extern apf::Mesh* getMesh(apf::Numbering* n);
extern void printNumberingName(apf::Numbering* n);

extern apf::MeshTag* createDoubleTag(apf::Mesh2 * m_local, char* name, int size);

extern void setDoubleTag(apf::Mesh2* m_local, apf::MeshEntity* e, apf::MeshTag* tag,  double data[]);

extern void getDoubleTag(apf::Mesh2* m_local, apf::MeshEntity* e, apf::MeshTag* tag,  double* data);



// mesh adaptation functions

extern IsotropicFunctionJ* createIsoFunc(apf::Mesh* m, apf::Field* f);
extern void deleteIsoFunc(IsotropicFunctionJ* isofunc);
extern ma::SolutionTransfers* createSolutionTransfers();
extern void deleteSolutionTransfers(ma::SolutionTransfers* soltrans);
extern void addSolutionTransfer(ma::SolutionTransfers* soltrans, apf::Field* f);
extern ma::Input* configureMAInput(apf::Mesh2* m, IsotropicFunctionJ* isofunc, 
                            ma::SolutionTransfer* soltrans);

extern void runMA(ma::Input* in);
extern void getAvgElementSize(apf::Mesh* m, apf::Numbering* el_N, double* el_sizes);

// apf::Field functions (needed for automagical solution transfer)
apf::Field* createPackedField(apf::Mesh* m, char* fieldname, int numcomponents, apf::FieldShape* fshape);

void setComponents(apf::Field* f, apf::MeshEntity* e, int node,  double const components[]);

void getComponents(apf::Field* f, apf::MeshEntity*e, int node, double components[]);


void zeroField(apf::Field* f);
apf::Field* getCoordinateField(apf::Mesh* m_ptr);

apf::Field* findField(apf::Mesh* m, char* fieldname);
void destroyField(apf::Field* f);
void destroyFields(apf::Mesh* m, apf::Field* save_n[], int n_save);

apf::FieldShape* getSBPShapes(int type, int order);


std::size_t countPeers(apf::Mesh* m, int dim);
void getPeers(apf::Mesh*m, int part_nums[]);

int isShared(apf::Mesh* m, apf::MeshEntity* e);

std::size_t countRemotes(apf::Mesh* m, apf::MeshEntity* e);
void getRemotes(int part_nums[], apf::MeshEntity* entities[]);

// mesh warping functions
void setPoint(apf::Mesh2* m, apf::MeshEntity* e, int node, double* coords);
void acceptChanges(apf::Mesh2* m);
void Verify(apf::Mesh* m);
void getPoint(apf::Mesh* m, apf::MeshEntity* e,  int node, double* coords);

// matched entity functions
bool hasMatching(apf::Mesh* m);
apf::Sharing* getSharing(apf::Mesh* m);
bool isOwned(apf::Sharing* shr, apf::MeshEntity* e);
std::size_t countCopies(apf::Sharing* shr, apf::MeshEntity* e);
void getCopies(int part_nums[], apf::MeshEntity* entities[]);


std::size_t countMatches(apf::Mesh* m, apf::MeshEntity* e);
void getMatches(int part_nums[], apf::MeshEntity* entities[]);


void getTopologyMaps(int* tri_edge_verts_in, int* tet_edge_verts_in, int* tet_tri_verts_in);
}  // extern c


int getindex_c(const int i, const int j, const int si, const int sj);

#endif
