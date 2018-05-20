#ifndef APFDG1SBPSHAPE_H
#define APFDG1SBPSHAPE_H

#include"apf.h"
#include "apfNew.h"
#include <sstream>
#include <iostream>
#include "apfShape.h"

namespace apf {
  extern "C" {
    FieldShape* getDG1SBPShape(int order);
    FieldShape* getDG7SBPShape(int order);
  }
}

#endif
