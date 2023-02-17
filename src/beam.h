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
  // TODO change these to non-STL raw arrays
  // for GPU compatibility
  std::vector<double> el_;  // element length
  std::vector<double> q_;   // degrees of freedom {q,q_d}
  std::vector<double> M_;   // global mass matrix
  std::vector<double> K_;   // global stiffness matrix
  std::vector<double> C_;   // global damping matrix
  std::vector<double> F_;   // global forcing vector
  std::vector<double> G_;   // Iteration matrix for time integration
  std::vector<double> elprop_; // element property data per quadrature pt

  beam();
  virtual ~beam();
};
}
// include child beam classes here
#include "EulerBernoulliBeam.h"
  


  
 

