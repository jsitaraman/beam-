namespace beam_pp {
namespace basis_lib
{
  void Hermite(double phi[4],  double phid[4], double phidd[4],
		 const double xp,const double xl)
  {
    double x=(1+xp)*0.5; // change to [0,1] space
    phi[0]=2*x*x*x-3*x*x+1;
    phi[1]=(x*x*x-2*x*x+x)*xl;
    phi[2]=(-2*x*x*x+3*x*x);
    phi[3]=(x*x*x-x*x)*xl;

    phid[0]=(6*x*x-6*x)/xl;
    phid[1]=(3*x*x-4*x+1);
    phid[2]=(-6*x*x+6*x)/xl;
    phid[3]=(3*x*x-2*x);

    phidd[0]=(12*x-6)/xl/xl;
    phidd[1]=(6*x-4)/xl;
    phidd[2]=(-12*x+6)/xl/xl;
    phidd[3]=(6*x-2)/xl;      
  }

  void  Lagrange(double phi[3], double phid[3],
		 const double xp, const double xl)
  {
    double x=(1+xp)*0.5; // change to [0,1] space
    phi[0]=2*x*x-3*x+1;
    phi[1]=-4*x*x+4*x;
    phi[2]=2*x*x-x;

    phid[0]=(4*x-3)/xl;
    phid[1]=(-8*x+4)/xl;
    phid[2]=(4*x-1)/xl;
  }
}
}
  
