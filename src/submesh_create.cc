// creates a new mesh from an existing mesh and a list of elements to 
// retain on the new mesh


#include "submesh_create.h"

SubMeshData::SubMeshData(apf::Mesh* _m_old, apf::Numbering* _numberings[], int* _el_list, int numel)  // useful constructor
{
  m_old = _m_old;
  dim = m_old->getDimension();
  
  numberings.resize(dim);
  for (int i = 0; i <= dim; ++i)
    numberings[i] = _numberings[i];

  // create new mesh
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  bool ismatched = false; //TODO
  m_new = apf::makeEmptyMdsMesh(g, dim, ismatched);
  m_new->changeShape(m_old->getShape(), false);

  // convert the el_list to a 0-based std::vector
  el_list.resize(numel);
  for (int i = 0; i < numel; ++i)
    el_list[i] = _el_list[i] - 1;

  // turn el_list into a set of MeshEntity*
  auto el_N = numberings[dim];
//    el_set.reserve(numel);
  el_entities.resize(numel);
  int idx = 0; // current index in el_entities
  apf::MeshIterator* it = m_old->begin(dim);
  apf::MeshEntity* e;
  while ( (e = m_old->iterate(it)) )
  {
    int elnum = apf::getNumber(el_N, e, 0, 0);
    // if this element is in el_list
    if ( std::binary_search(el_list.begin(), el_list.end(), elnum) )
    {
      el_entities[idx] = e;
      idx += 1;
       
//        el_set.insert(e);
    }
  }

  assert(idx == numel);  // all entities found

  // fill verts, edges, faces, regions with NULL

  verts = std::vector<apf::MeshEntity*>(m_old->count(0), NULL);
  edges = std::vector<apf::MeshEntity*>(m_old->count(1), NULL);
  faces = std::vector<apf::MeshEntity*>(m_old->count(2), NULL);
  regions = std::vector<apf::MeshEntity*>(m_old->count(3), NULL);
  entities[0] = verts;
  entities[1] = edges;
  entities[2] = faces;
  entities[3] = regions;  

}  // constructor

//-----------------------------------------------------------------------------
// julia iterface


// el_list is the 1-based list of elements to preserve (must be sorted in
// ascending order)
// numberings[dim] is the array of zero-based apf::Numbering objects that n
// numbers the MeshEntities from dimension 0 to dim
SubMeshData* createSubMesh(apf::Mesh* m, apf::Numbering* numberings[],
                          int* el_list, int numel)
{
  SubMeshData* sdata = new SubMeshData(m, numberings, el_list, numel);
  createVertices(sdata);
  createEntities(sdata);

  //TODO: fix geometry classification

  return sdata;
}

void writeNewMesh(SubMeshData* sdata, const char* fname)
{
  sdata->writeNewMesh(fname);
}

apf::Mesh2* getNewMesh(SubMeshData* sdata)
{
  return sdata->getNewMesh();
}


//-----------------------------------------------------------------------------
// helper functions

// check that this element list can be used with this mesh
void checkInput(SubMeshData* sdata)
{
  // if an assert below triggers, this may leak memory, but thats ok

  // check there are fewer elements in the list than in the old mesh
  int max_el = 0;

  // verify el_list is sorted
  bool failflag = false;
  int prev_val = -1;
  for (std::vector<int>::iterator it = sdata->el_list.begin(); 
                                  it != sdata->el_list.end(); ++it)
  {
    failflag = failflag || (*it <= prev_val);
    prev_val = *it;

    if (*it > max_el)
      max_el = *it;
  }

  assert(max_el < (int)sdata->m_old->count(sdata->dim));
  assert(!failflag);
}  // checkInputs

// create vertocies on the new mesh, populates sdata->verts
void createVertices(SubMeshData* sdata)
{
  // Algorithm: loop over elements to be copied, see
  //            if vertex already created, if not create the entity

  apf::Downward down_verts;
  apf::Mesh* m_old = sdata->m_old;
  apf::MeshEntity* e;  // current vert
  apf::Numbering* vert_N = sdata->numberings[0];

  for (std::vector<apf::MeshEntity*>::iterator it = sdata->el_entities.begin();
       it != sdata->el_entities.end(); ++it)

  {
    e = *it;

    // get downard vertices
    int nverts = m_old->getDownward(e, 0, down_verts);
    
    for (int i = 0; i < nverts; ++i)
    {
      int vertnum = apf::getNumber(vert_N, down_verts[i], 0, 0);
      if ( sdata->verts[vertnum] == NULL) // vert not yet created
        sdata->verts[vertnum] = createVert(sdata, down_verts[i]);
    }
  }  // loop over elements

}  // createEntities

void createEntities(SubMeshData* sdata)
{
  // Algorithm: loop over elements, loop over dimensions 1 to dim,
  //            if entity not yet created, create it from lower adjacencies

  apf::Downward down_entities;
  apf::Mesh* m_old = sdata->m_old;
  auto dim = sdata->dim;

  for (std::vector<apf::MeshEntity*>::iterator it = sdata->el_entities.begin();
       it != sdata->el_entities.end(); it++)
    for (int d = 1; d <= dim; ++d)
    {
      int nentities = m_old->getDownward(*it, d, down_entities);
      for (int i = 0; i < nentities; ++i)
      {
        int entitynum = apf::getNumber(sdata->numberings[d], down_entities[i], 0, 0);
        if (sdata->entities[d][entitynum] == NULL)
          sdata->entities[d][entitynum] = createEntity(sdata, down_entities[i]);

      }  // loop i
    } // loop d

}  // createEntities


// create a vertex on the new mesh given a vertex on the old mesh
// uses same geometry classification as original mesh
apf::MeshEntity* createVert(SubMeshData* sdata, apf::MeshEntity* vert)
{

  // get the old geometry
  apf::ModelEntity* me = sdata->m_old->toModel(vert);
  int me_type = sdata->m_old->getModelType(me);
  int me_tag = sdata->m_old->getModelTag(me);

  // get/create the geometry on the new mesh
  // the new mesh was create with gmi_null geometry model, so it will create
  // geometric entities as needed

  apf::ModelEntity* me_new = sdata->m_new->findModelEntity(me_type, me_tag);
  apf::MeshEntity* vert_new = sdata->m_new->createVert(me_new);

  // set coordinates
  apf::Vector3 coords;
  sdata->m_old->getPoint(vert, 0, coords);
  sdata->m_new->setPoint(vert, 0, coords);

  return vert_new;
}  // createVert

// create all remaining entities from lower adjacencies
apf::MeshEntity* createEntity(SubMeshData* sdata, apf::MeshEntity* entity)
{
  // get the old geometry
  apf::ModelEntity* me = sdata->m_old->toModel(entity);
  int me_type = sdata->m_old->getModelType(me);
  int me_tag = sdata->m_old->getModelTag(me);

  // get/create the geometry on the new mesh
  apf::ModelEntity* me_new = sdata->m_new->findModelEntity(me_type, me_tag);

  // get the entity type and the lower adjacent entities that define it
  auto e_type = sdata->m_old->getType(entity);
  int dim = apf::Mesh::typeDimension[e_type];
  apf::Downward down_entities;
  int nentities = sdata->m_old->getDownward(entity, dim-1, down_entities);

  // get the corresponding entities on the new mesh
  for (int i = 0; i < nentities; ++i)
  {
    int e_num = apf::getNumber(sdata->numberings[dim-1], down_entities[i], 0, 0);
    down_entities[i] = sdata->entities[dim-1][e_num];
    assert(down_entities[i] != NULL); // entity already created
  }

  // create the entity
  apf::MeshEntity* e_new = sdata->m_new->createEntity(e_type, me_new, down_entities);

  // for curvilinear meshes, set the higher order node locations (if any)
  apf::FieldShape* fshape = sdata->m_old->getShape();

  if (fshape->hasNodesIn(dim))
  {
    int nnodes = fshape->countNodesOn(e_type);
    apf::Vector3 coords;
    for (int i = 0; i < nnodes; ++i)
    {
      sdata->m_old->getPoint(entity, i, coords);
      sdata->m_new->setPoint(e_new, i, coords);
    }
  }

  return e_new;
}  // function createEntity

