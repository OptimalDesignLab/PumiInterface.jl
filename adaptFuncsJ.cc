#include "adaptFuncsJ.h"
#include "math.h"
// define member functions for mesh adaptation functions
//


// return the size of element desired at a particular vertex
// this function can access any data members of the class
double IsotropicFunctionJ::getValue(apf::MeshEntity* vert)
{
/*
  apf::Vector3 coords;
  m->getPoint(vert, 0, coords);

  double x_coord = coords[0];
  double frac = (x_coord + 1.5)/2.0; // distance from left edge
  double h_value = 1.25*frac;

  std::cout << "coords = " << coords << " , h_value = " << h_value << std::endl;

  if (h_value <= 0.0)
    std::cout << "    Warning: h_value is negative " << std::endl;
*/
  double h_value = (*juliafunc)(vert, m, u);
  std::cout << "from c++, h_value = " << h_value << std::endl;
  return h_value;
}

// populates vector h with desired element size in each direction
// 
void AnisotropicFunctionJ::getValue(apf::MeshEntity* vert, ma::Matrix &r, ma::Vector &h)
{
  // copy ma::Matrix to C standard arrays
  r.toArray(r_array);
  std::cout << "from C++, r = " << r << std::endl;

  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      std::cout << "r_array " << i << " , " << j << " = " << r_array[i][j] << std::endl;
    }
  }
  
  // call julia function to populate h_array
  (*juliafunc)(vert, r_array, h_array, m, u);


  h.fromArray(h_array);
  std::cout << "after sizefunc, h = " << h << std::endl;
//  r.fromArray(r_array);
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      r[i][j] = r_array[i][j];
      std::cout << " r[i][j] = " << r[i][j] << std::endl;
    }
  }
 

  std::cout << "after size function, h = " << h << std::endl;
}
