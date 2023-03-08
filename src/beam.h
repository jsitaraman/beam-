#include<cstdio>
#include<cstdlib>
#include<vector>
#include<unordered_map>
namespace beam_pp {
// basic beam class
class beam
{
public:
  int nelem_{0};            // number of elements
  int dofs_per_elem_{0};    // number of dofs per element
  int ndofs_{0};            // total degrees of freedom
  int nprops_{0};           // number of properties given
  double L_{1};             // length of the beam
  double root_{0};          // root location of the beam
  double omega_{0};         // rotation rate
  // raw arrays
  // for GPU compatibility
  double * el_{nullptr};     // element length
  double * q_{nullptr};      // degrees of freedom {q,q_d}
  double * M_{nullptr};      // global mass matrix
  double * K_{nullptr};      // global stiffness matrix
  double * C_{nullptr};      // global damping matrix
  double * F_{nullptr};      // global forcing vector
  double * G_{nullptr};      // Iteration matrix for time integration
  double * elprop_{nullptr}; // element property data per quadrature pt

  beam();
  virtual ~beam();
};
}
// include child beam classes here
#include "EulerBernoulliBeam.h"
  


  
 

