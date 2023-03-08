#include "beam.h"
#include<vector>
#include<unordered_map>
#include<iostream>
#include<cassert>
using namespace beam_pp;
namespace beam_pp {
// default constructor  
beam::beam() { };
// default destructor place holder
beam::~beam() {
  if (el_) delete [] el_;
  if (q_) delete [] q_;
  if (M_) delete [] M_;
  if (K_) delete [] K_;
  if (C_) delete [] C_;
  if (F_) delete [] F_;
  if (G_) delete [] G_;
  if (elprop_) delete [] elprop_;
};
}
