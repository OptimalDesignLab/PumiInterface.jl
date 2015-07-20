
// initilize global variables, used by all fucntions
// downward_counts = numDown
// number_entities = numEntity
// m_ptr = mesh pointer
// mshape_ptr = pointer to mesh shape (m->getShape())
int initABC(char* dmg_name, char* smb_name, int downward_counts[4][4], int number_entities[4], apf::Mesh2* m_ptr_array[1], apf::FieldShape* mshape_ptr_array[1] )
{
  std::cout << "Entered init\n" << std::endl;

  MPI_Init(0,NULL);  // initilize MPI 
  PCU_Comm_Init();   // initilize PUMI's communication
//  gmi_register_mesh();

  // load mesh
//  m = apf::loadMdsMesh("cube.dmg", "tet-mesh-1.smb");


  if (strcmp(dmg_name, ".null") == 0)
  {
    gmi_register_null();
    std::cout << "loading null geometric model" << std::endl;
    gmi_model* g = gmi_load(".null");
    std::cout << "finished loading geometric model" << std::endl;
    m = apf::loadMdsMesh(g, smb_name);
  } else {
    gmi_register_mesh();
    std::cout << "loading geometric model from file" << std::endl;
    m = apf::loadMdsMesh(dmg_name, smb_name);
  }


  apf::writeVtkFiles("output_check", m);
  std::cout << "finished loading mesh" << std::endl;
  m_ptr_array[0] = m;
  mshape_ptr_array[0] = m->getShape();
  std::cout << std::endl;
/* 
  // initilize iterators
  elIt = m->begin(3);
  its[2] = m->begin(2);
  its[1] = m->begin(1);
  its[0] = m->begin(0);
  
*/
  // initilize  number of each type of entity
  for (int i = 0; i < 4; ++i)
  {
    numEntity[i] = apf::countOwned(m, i);
    number_entities[i] = numEntity[i];
  }



  // initilize numberings
  numberings[0] = numberOwnedDimension(m, "vertNums", 0);
  numberings[1] = numberOwnedDimension(m, "edgeNums", 1);
  numberings[2] = numberOwnedDimension(m, "faceNums", 2);
  numberings[3] = numberOwnedDimension(m, "elNums", 3);

  // initalize tags
//  globalVertNums = m->createIntTag("globalNodeNumber", 1);


  // initilize iterators
  its[0] = m->begin(0);
  its[1] = m->begin(1);
  its[2] = m->begin(2);
  its[3] = m->begin(3);
  for (int i = 0; i < 4; ++i)
  {
    int type = m->getType(m->deref(its[i]));
    std::cout << "type of its[" << i << "] = " << type;
    std::cout << "index its[" << i << "] = " << type << std::endl;
  }

    
  // initilize number of downward adjacencies each type has, assuming all
  // elements have same number of downward adjacencies
  apf::MeshEntity* e_tmp;
  apf::Downward tmp;
  for (int i = 0; i < 4 ; ++i)
  {
    e_tmp = m->deref(its[i]);
    for (int j = 0; j < 4; ++j)
    {
      if ( j < i)
      {
        numDown[i][j] = m->getDownward(e_tmp, j, tmp);
        downward_counts[j][i] = numDown[i][j];
      }
      else
      {
        numDown[i][j] = 0;
        downward_counts[j][i] = 0;
      }
      std::cout << names[i] << " has " << numDown[i][j];
      std::cout << " downward adjacencies of type " << names[j] << std::endl;
    }
  }

  std::cout << std::endl;
/*
  numberings[0] = numberOwnedDimension(m, "elNums", 0);
  numberings[1] = numberOwnedDimension(m, "elNums", 1);
  numberings[2] = numberOwnedDimension(m, "elNums", 2);
  numberings[3] = numberOwnedDimension(m, "elNums", 3);
*/

  std::cout << "numV = " << numEntity[0] << " , numEdge = " << numEntity[1];
  std::cout << " , numFace = " << numEntity[2] << " , numEl = " << numEntity[3] << std::endl;
  std::cout << std::endl;

  apf::MeshEntity* e = m->deref(its[2]);
  apf::MeshElement* e_el = apf::createMeshElement(m, e);
  int numI = apf::countIntPoints(e_el, 3);
//  std::cout << numI << " integrations points required for 5th order accuracy" << std::endl;

  for ( int i = 0; i < numI; ++i)
  {
    apf::Vector3 coords;  // declare vector to hold coordinates
    apf::Matrix3x3 mat;
    apf::getIntPoint(e_el, 3, i, coords);
    apf::getJacobian(e_el, coords, mat);

    std::cout << "point " << i << " has coordinates " << coords << std::endl;
    std::cout << "  and jacobian = \n" << mat << std::endl;
  }
/*
  e = m->deref(its[1]);
  e_el = apf::createMeshElement(m,e);
  numI = apf::countIntPoints(e_el,5);
  std::cout << numI << " integration points required on edge for 5th order accuraccy" << std::endl;



  for ( int i = 0; i < numI; ++i)
  {
    apf::Vector3 coords;  // declare vector to hold coordinates
    apf::getIntPoint(e_el, 5, i, coords);
    std::cout << "point " << i << " has coordinates " << coords << std::endl;
  }

*/


  apf::writeVtkFiles("output_init", m);

  return 0;
}


