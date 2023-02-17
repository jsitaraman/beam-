%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEM for beams
% rotation terms not added (part of homework)
%
% Jay Sitaraman
% 3/8/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% nrel 5MW turbine data
%
EI0=[240E6];     % stiffness of beam
m0=[200];        % mass per unit length of beam
Radius=1;    % radius of beam
%
nelements=20; % number of elements
ll=ones(nelements,1)*Radius/nelements;
EI=ones(nelements,1)*EI0;
m=ones(nelements,1)*m0;
%
ngauss=6;
xgauss=[0.03376524300000
	0.16939530700000
        0.38069040700000
        0.61930959300000
        0.83060469300000
        0.96623475700000];

weights=[0.08566224600000
         0.18038078650000
         0.23395696750000
         0.23395696750000
         0.18038078650000
         0.08566224600000];
%
ndof=4;                 % number of degrees of freedom in each element
nvar=2*nelements+2;    % total number of degrees of freedom
%
M=zeros(nvar,nvar); % global mass matrix
K=zeros(nvar,nvar); % global stiffness matrix
%
for el=1:nelements  
  Kelem=zeros(4,4);
  Melem=zeros(4,4);
  for i=1:ndof
   for j=1:ndof
     for n=1:6
      x=xgauss(n);
      Hi=shapeFunc(x,i,ll(el));
      Hj=shapeFunc(x,j,ll(el));
      dHi=shapeFuncDerivative(x,i,ll(el));
      dHj=shapeFuncDerivative(x,j,ll(el));
      ddHi=shapeFuncDoubleDerivative(x,i,ll(el));
      ddHj=shapeFuncDoubleDerivative(x,j,ll(el));
      Kelem(i,j)=Kelem(i,j)+weights(n)*EI(el)*ddHi*ddHj*ll(el);
      Melem(i,j)=Melem(i,j)+weights(n)*m(el)*Hi*Hj*ll(el);
     end
    end
   end
   %
   % location of elements in 
   % the global stiffness matrix
   %
   iloc=2*el-1;
   jloc=2*el-1;  
   %
   K(iloc:iloc+3,jloc:jloc+3)=K(iloc:iloc+3,jloc:jloc+3)+Kelem;
   M(iloc:iloc+3,jloc:jloc+3)=M(iloc:iloc+3,jloc:jloc+3)+Melem;
   %imagesc(K);
   %input('go on?')
end
%
% apply boundary conditions
%
K=K(3:nvar,3:nvar);
M=M(3:nvar,3:nvar);
%
% find natural frequencies and modes
%
[v,wn]=eig(K,M);
wn=sqrt(wn);
