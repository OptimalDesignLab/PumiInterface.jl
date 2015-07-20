#ifndef FUNCS1_H
#define FUNCS1_H

// header file for funcs1


/*
 * Params[in]:
 * pointer to mesh m
 */
extern "C" {
// this function does not pass opaque pointers
//extern int initABC(char* dmg_name, char* smb_name, int downward_counts[4][4], int numberEntities[4], apf::Mesh2* m_ptr_array[1], apf::FieldShape* mshape_ptr_array[1]);


int initABC(char* dmg_name, char* smb_name, int number_entities[4], apf::Mesh2* m_ptr_array[1], apf::FieldShape* mshape_ptr_array[1], int order, int load_mesh, int shape_type );


int initABC2(char* dmg_name, char* smb_name, int number_entities[3], apf::Mesh2* m_ptr_array[1], apf::FieldShape* mshape_ptr_array[1], int order, int load_mesh, int shape_type );

// these functions are not user accessible
void cleanup(apf::Mesh* m_local);
void destroyNumberings(int dim); 

// these functions do pass pointers
extern apf::Mesh2* getMeshPtr();
extern apf::FieldShape* getMeshShapePtr();
extern apf::FieldShape* getConstantShapePtr(int dimension);
extern apf::Numbering* getVertNumbering();
extern apf::Numbering* getEdgeNumbering();
extern apf::Numbering* getFaceNumbering();
extern apf::Numbering* getElNumbering();

// these functions do not pass pointers
extern void resetVertIt();
extern void resetEdgeIt();
extern void resetFaceIt();
extern void incrementVertIt();
extern void incrementVertItn(int n);

extern void incrementEdgeIt();
extern void incrementEdgeItn(int n);

extern void incrementFaceIt();
extern void incrementFaceItn(int n);

extern void incrementElIt();
extern void incrementElItn(int n);

extern void resetElIt();
extern int count(apf::Mesh2* m_local, int dimension);
extern void writeVtkFiles(char* name, apf::Mesh2* m_local);

extern void setGlobalVertNumber(int val); // do not use
extern int getGlobalVertNumber(); // do not use
extern int getVertNumber();
extern int getEdgeNumber();
extern int getFaceNumber();
extern int getElNumber();

// these functions pass pointers to/from julia
extern apf::MeshEntity* getVert();
extern apf::MeshEntity* getEdge();
extern apf::MeshEntity* getFace();
extern apf::MeshEntity* getEl();

// these functions are deprecated, use getNumberJ insteady
extern int getVertNumber2(apf::MeshEntity* e);
extern int getEdgeNumber2(apf::MeshEntity* e);
extern int getFaceNumber2(apf::MeshEntity* e);
extern int getElNumber2(apf::MeshEntity* e);


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



void getAlignment(apf::Mesh* m_local, apf::MeshEntity* elem, apf::MeshEntity* boundary, int which[1], bool flip[1], int rotate[1]);

// FieldShape functions
bool hasNodesIn(apf::FieldShape* fshape_local, int dimension);
extern int countNodesOn(apf::FieldShape* mshape_ptr, int type);
extern apf::EntityShape* getEntityShape(apf::FieldShape* mshape_local, int type);

// MeshElement related functions
extern apf::MeshElement* createMeshElement( apf::Mesh2* m_local, apf::MeshEntity* e);
extern int countIntPoints(apf::MeshElement* e, int order);
extern void getIntPoint(apf::MeshElement* e, int order, int point, double coords[3]);
extern double getIntWeight(apf::MeshElement* e, int order, int point);
extern void getJacobian(apf::MeshElement* e, double coords[3], double jac[3][3]);

// Entityshape functions
extern int countNodes(apf::EntityShape* eshape_local);
extern void getValues(apf::EntityShape* eshape_local, double xi[3], double vals[]);
extern void getLocalGradients(apf::EntityShape* eshape_local, double xi[3], double vals[][3]);
void alignSharedNodes(apf::EntityShape* eshape_local, apf::Mesh* m_local, apf::MeshEntity* elem, apf::MeshEntity* shared, int order[]);

// these function do not pass pointers (they use the iterators internally)
extern void checkVars();
extern void checkNums();
extern void getVertCoords(double coords[][3], int sx, int sy);
extern int getEdgeCoords(double coords[2][3], int sx, int sy);
extern int getFaceCoords(double coords[][3], int sx, int sy);
extern int getElCoords(double coords[][3], int sx, int sy);
extern int getElCoords2(apf::MeshEntity* e, double coords[][3], int sx, int sy);

// these functinos pass pointers
extern void getVertCoords2(apf::MeshEntity* e, double coords[][3], int sx, int sy);
extern int getEdgeCoords2(apf::MeshEntity* e, double coords[2][3], int sx, int sy);
int getFaceCoords2(apf::MeshEntity* e, double coords[][3], int sx, int sy);



// these function pass pointers
// create a generally defined numbering from julia
extern  apf::Numbering* createNumberingJ(apf::Mesh2* m_local, char* name, apf::FieldShape* field, int components);
extern int numberJ(apf::Numbering* n, apf::MeshEntity* e, int node, int component, int number);
extern  int getNumberJ(apf::Numbering* n, apf::MeshEntity* e, int node, int component);

void getElementNumbers(apf::Numbering* n, apf::MeshEntity*e, int num_nodes, int nums[]);
extern apf::Mesh* getMesh(apf::Numbering* n);
extern void printNumberingName(apf::Numbering* n);

extern apf::MeshTag* createDoubleTag(apf::Mesh2 * m_local, char* name, int size);

extern void setDoubleTag(apf::Mesh2* m_local, apf::MeshEntity* e, apf::MeshTag* tag,  double data[]);

extern void getDoubleTag(apf::Mesh2* m_local, apf::MeshEntity* e, apf::MeshTag* tag,  double* data);



// mesh adaptation functions
extern void createIsoFunc(apf::Mesh2* m_local, double(*sizefunc)(apf::MeshEntity*vert, apf::Mesh2* m_local, double *u), double *u);
extern void runIsoAdapt(apf::Mesh2* m_local);
extern void createAnisoFunc(apf::Mesh2* m_local,  void (*sizefunc)(apf::MeshEntity* vert, double r[3][3], double h[3], apf::Mesh2* m_ptr, void *f_ptr, double *operator_ptr), apf::Field *f_ptr, double *operator_ptr);
void runAnisoAdapt(apf::Mesh2* m_local);


// apf::Field functions (needed for automagical solution transfer)
apf::Field* createPackedField(apf::Mesh* m, char* fieldname, int numcomponents);

void setComponents(apf::Field* f, apf::MeshEntity* e, int node,  double const components[]);

void getComponents(apf::Field* f, apf::MeshEntity*e, int node, double components[]);

}

#endif
