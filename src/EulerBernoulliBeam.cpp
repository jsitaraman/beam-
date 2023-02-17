#include "beam.h"
#include "basis_lib.h"
#include "quadrature.h"
#include<vector>
#include<unordered_map>
#include<iostream>
#include<cassert>
using namespace beam_pp;
namespace beam_pp {
// helper enum for picking properties from data array
// everything is prefixed by  "i" here for preventing
// conflict with other variables
enum pmap {iMass, iEIy, iGJ, iEIz, iEA, iKM1, iKM2, iEG, iCF, iXL, iXX};  
namespace element_matrices {
  // Mass matrix for an Euler-Bernoulli Beam
  void mass_matrix(const int elemID,      // element id
		   const int i,           // ith dof of element
		   const int j,           // jth dof of element
		   const int nelem,       // number of elements
		   const double xl,       // element length
		   const double *elprops, // element property data
		   const double omega,    // rotational speed
		   const int iloc,        // location in global matrix to fill in 
		   double *M)             // Mass matrix 
  {
    int ngauss=quadrature::nquadrature;
    // a small lambda to get properties
    auto getProp=[ngauss,nelem,elprops](int iprop,int idx)
    {
      return elprops[iprop*nelem*ngauss+idx];
    };
    // storage for basis functions
    double phiH[4],phiHd[4],phiHdd[4];
    double phiL[4],phiLd[4];
    // flapping terms
    if (i < 4 && j < 4) {
      double val =0 ;
      for(int n=0;n<ngauss;n++)
	{
	  int idx=elemID*ngauss+n;
	  double mass=getProp(iMass,idx);
	  double xl=getProp(iXL,idx);
	  double xp=quadrature::xloc[n];
	  basis_lib::Hermite(phiH,phiHd,phiHdd,xp,xl);
	  double wg = quadrature::weight[n]*0.5*xl;
	  val+=(mass*wg*phiH[i]*phiH[j]);
	}
      // replace by ops::atomic_Add(M[iloc],val)
      // for GPUs
      M[iloc]+=val;
    }
    // flap torsion terms
    else if (i < 4 && (j-4)*(j-6)<=0) {
      double val = 0;
      for(int n=0;n<ngauss;n++) {
	int idx=elemID*ngauss+n;
	double xl=getProp(iXL,idx);
	// TODO : template this to call only the needed
	// shape function
	double xp=quadrature::xloc[n];
	basis_lib::Hermite(phiH,phiHd,phiHdd,xp,xl);
	basis_lib::Lagrange(phiL,phiLd,xp,xl);
	//
	double mass=getProp(iMass,idx);
	double km1=getProp(iKM1,idx);
	double km2=getProp(iKM2,idx);
	double eg=getProp(iEG,idx);
	double wg=quadrature::weight[n]*0.5*xl;
	val+=wg*mass*eg*phiH[i]*phiL[j-4];
      }
      M[iloc]+=val;
    }
    // torsion terms
    else if (i >=4 && j >=4) {
      double val = 0;
      for(int n=0;n<ngauss;n++) {
	int idx=elemID*ngauss+n;
	double xl = getProp(iXL,idx);
	double mass= getProp(iMass,idx);
	double km1 = getProp(iKM1,idx);
	double km2 = getProp(iKM2,idx);
	double kmsq = km1+km2;
	double wg   = quadrature::weight[n]*0.5*xl;
	double xp   = quadrature::xloc[n];
	// TODO: template this to call only the
	// needed shape function	
	basis_lib::Lagrange(phiL,phiLd,xp,xl);
	val+=wg*mass*kmsq*phiL[i-4]*phiL[j-4];
      }
      M[iloc]+=val;
    }
  }
  // damping matrix for EB element
  void damping_matrix(const int elemID,      // element id
		      const int i,           // ith dof of element
		      const int j,           // jth dof of element
		      const int nelem,       // number of elements
		      const double xl,       // element length
		      const double *elprops, // element property data
		      const int iloc,        // location in global matrix to fill in
		      const double str_damp, // structural damping  
		      double *C)             // Mass matrix 
  {
    // just diagonal structural damping for now for stability
    if (i==j ) {
      C[iloc]+=str_damp;
    }
  }
  // stiffness matrix for EB element
  void stiffness_matrix(const int elemID,    // element id
			const int i,           // ith dof of element
			const int j,           // jth dof of element
			const int nelem,       // number of elements
			const double xl,       // element length
			const double *elprops, // element property data
         		const double omega,    // rotational speed 
			const int iloc,        // location in global matrix to fill in
			double *K)             // Stiffness matrix
  {
    int ngauss=quadrature::nquadrature;
    // a small lambda to get properties
    auto getProp=[ngauss,nelem,elprops](int iprop,int idx)
    {
      return elprops[iprop*nelem*ngauss+idx];
    };
    double phiH[4],phiHd[4],phiHdd[4];
    double phiL[4],phiLd[4];    
    // flapping
    if (i < 4 && j < 4) {
      double val =0 ;
      for(int n=0;n<ngauss;n++) {
	int idx=elemID*ngauss+n;
	double xl = getProp(iXL,idx);
	double cf = getProp(iCF,idx);
	double EIy = getProp(iEIy,idx);
	double xp = quadrature::xloc[n];
	basis_lib::Hermite(phiH,phiHd,phiHdd,xp,xl);
	// centrifugal stiffness term
	double wg=quadrature::weight[n]*0.5*xl;
	val+=wg*phiHd[i]*phiHd[j]*cf;
	// structural stiffness term
	val+=wg*EIy*phiHdd[i]*phiHdd[j];
      }
      K[iloc]+=val;
    }
    // flap torsion coupling terms
    if (i < 4 && (j-4)*(j-6)<=0) {
      double val=0;
      for(int n=0;n<ngauss;n++) {
	int idx=elemID*ngauss+n;
	double xl = getProp(iXL,idx);
	double xx = getProp(iXX,idx);
	double cf = getProp(iCF,idx);
	double mm = getProp(iMass,idx);
	double eg = getProp(iEG,idx);
	double xp = quadrature::xloc[n];
	basis_lib::Hermite(phiH,phiHd,phiHdd,xp,xl);
	basis_lib::Lagrange(phiL,phiLd,xp,xl);
	// centrifugal stiffness term
	double wg=quadrature::weight[n]*0.5*xl;
	val+=wg*phiHd[i]*phiHd[j]*cf;
	// structural stiffness term
	val+=wg*mm*xx*phiHd[i]*phiL[j-4]*omega*omega*eg;
      }
      K[iloc]+=val;
    }
    // torsion stiffness terms
    if (i >=4 && j >=4) {
	double val=0;
	for(int n=0;n<ngauss;n++) {
	  int idx=elemID*ngauss+n;
	  double  xl=getProp(iXL,idx);
	  double km1=getProp(iKM1,idx);
	  double km2=getProp(iKM2,idx);
	  double  wg=quadrature::weight[n]*0.5*xl;
	  double mass=getProp(iMass,idx);
	  double GJ  = getProp(iGJ,idx);
	  double kmsq=km1+km2;
	  val+=wg*mass*omega*omega*kmsq*phiL[i-4]*phiL[j-4];
	  val+=wg*GJ*phiLd[i-4]*phiLd[j-4];
	}
	K[iloc]+=val;
    }
  }  
}
using dict_properties = std::unordered_map<std::string,
					   std::vector<double>>;
// constructor for EulerBernoulliBeam  
EulerBernoulliBeam::EulerBernoulliBeam(dict_properties & input_data)
{  
  int n_required_data=6;
  int n_possible_data=4;
  std::string required_data[6]={"x","Mass","root","nelem","omega","L"};
  for(int idx=0;idx<n_required_data;idx++) {
    if (input_data.find(required_data[idx])==input_data.end())
      {
	std::cout << "Property " << required_data[idx] <<
	  " is not provided\n";
	exit(0);
      }
  }   
  root_=input_data["root"].data()[0];
  nelem_=input_data["nelem"].data()[0];
  omega_=input_data["omega"].data()[0];
  L_=input_data["L"].data()[0];
  nprops_=1;
  std::string possible_data[4]=
    {"EIy","GJ","EIz","EA"};
  std::string geometric_data[3]={ "yk1","yk2","eg"};
  int dof_count[4]={4,3,4,3};
  dofs_per_elem_=0;
  for(int idx=0;idx<4;idx++) {
    if (input_data.find(possible_data[idx])!=input_data.end())
      {
	dofs_per_elem_+=dof_count[idx];
	nprops_++;
      }
  }
  // make uniform size elements
  double deltal = (L_-root_)/nelem_;
  std::cout << "deltal" << " " << deltal << std::endl;
  int ngauss=quadrature::nquadrature;
  // uniform size elements for now
  el_.resize(nelem_,deltal);
  elprop_.resize(nelem_*ngauss*nprops_max_,0);
  double *x_by_l=input_data["x"].data();
  // non-dimensionalize the length
  int ndata=input_data["x"].size();
  for(int i=0;i<ndata;i++) x_by_l[i]/=L_;
  //
  std::unordered_map<std::string, int> prop_map;
  // create a helper map for filling in property data
  // based on input strings
  prop_map["Mass"]=iMass;
  prop_map["EIy"]=iEIy;
  prop_map["EIz"]=iEIz;
  prop_map["GJ"]=iGJ;
  prop_map["km1"]=iKM1;
  prop_map["km2"]=iKM2;
  prop_map["eg"]=iEG;
  // a small lambda to get properties
  auto setProp=[&](int iprop,int idx,double value)
   {
     elprop_[iprop*nelem_*ngauss+idx]=value;
   };
  //
  double totalmass=0;
  double xi=root_;
  for(int e=0;e<nelem_;e++) {
    for(int n=0;n<ngauss;n++)
      {
	double wg=quadrature::weight[n]; // gaussian weight
	double xp=quadrature::xloc[n];   // quadrature locatio
	double xx=xi+el_[e]*(1+xp)*0.5;  // x location
	double xx_by_L=xx/L_;            // non-dimensional location
	int idx = e*ngauss+n;            // index for storage
	int k;
	// search to find linear interpolation basis
	for(k=0;k<ndata;k++) {
	  if (xx_by_L < x_by_l[k]) break;
	}
	double scal=(xx_by_L - x_by_l[k-1])/(x_by_l[k]-x_by_l[k-1]);
	// build a guard here to make sure data is ok
	assert(k > 0);
	assert(xx_by_L < x_by_l[ndata-1]);
	std::cout << xx ;
	for(auto &prop : input_data) {
	  if (prop.first=="x" || prop.first=="root" ||
	      prop.first=="nelem"|| prop.first=="L"  || prop.first =="omega") continue;
	  std::cout << " " << prop.first;
	  int iprop=prop_map[prop.first];
	  double *data=prop.second.data();
	  // linear array to store the element properties
	  // at quadrature locations
	  double value = data[k-1]*(1-scal)+data[k]*scal;
	  std::cout << " " << value;
	  setProp(iprop,idx,value);
	  if (prop.first == "Mass") {
	    totalmass += 0.5*value*wg*el_[e];
	  }
	}
	std::cout << '\n';
	setProp(iXL,idx,el_[e]);
	setProp(iXX,idx,xx);
      }
    xi+=el_[e];
  }
  std::cout << "Beam Total Mass :"<< totalmass << std::endl;
  centrifugal_stiffening();
}
void EulerBernoulliBeam::centrifugal_stiffening()
{
  int ngauss=quadrature::nquadrature;
  for(int e=nelem_-1;e>=0;e--) {
    double cf=0;
    double xi=L_;
    // integrate the centrifugal force with quadrature
    // within each element
    for(int f=nelem_-1;f>=e+1;f--) {
      for(int n=0;n<ngauss;n++) {
	double wg=quadrature::weight[n]; // gaussian weight
	double xp=quadrature::xloc[n];   // quadrature locatio
	double xx=xi+el_[f]*(1+xp)*0.5;  // x location
	double mass=elprop_[iMass*nelem_*ngauss+f*ngauss+n];
	cf+=wg*0.5*el_[e]*xx*omega_*omega_;
      }
      xi-=el_[e];
    }
    // assume linear property distribution within the element
    // this integration needs to be changed if we spline the
    // properties
    double xend = xi;
    // mass property at last gauss point
    double x2=xi+(1+quadrature::xloc[ngauss-1])*0.5;
    double m2=elprop_[iMass*nelem_*ngauss+e*ngauss+ngauss-1];
    // mass property at first gauss point
    double x1=xi+(1+quadrature::xloc[0])*0.5;
    double m1=elprop_[iMass*nelem_*ngauss+e*ngauss];
    // slope of mass distribution
    double b=(m2-m1)/(x2-x1);
    // constant term
    double a=m1-b*x1;
    // m = a+bx
    // \int_x^xend (a+b*x)*x*dx
    // a*(xend^2-x^2)/2 + b*(xend^3-x^3)/3
    for(int n=ngauss-1;n>=0;n--) {
      double xp=quadrature::xloc[n];   
      double xx=xi+el_[e]*(1+xp)*0.5;  // x location
      elprop_[iCF+nelem_*ngauss+e*ngauss + n]=  cf +
	a*(xend*xend - xx*xx)*0.5 +
	b*(xend*xend*xend - xx*xx*xx)/3.0;
    }  
  }
}  
}
                                
				     
				     
				     

				     

				     

				    
