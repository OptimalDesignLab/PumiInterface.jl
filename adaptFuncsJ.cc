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
