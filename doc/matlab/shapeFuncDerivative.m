function dH=shapeFuncDerivative(x,i,l)
% x  !< location along the element
% i  !< row number
% l !< length of the element
% H  !< value of shape function at the given location
  
l1=1./l;
%
if (i==1) 
dH=(6*x^2-6*x)*l1;
elseif (i==2)
   dH=(3*x^2-4*x+1);
elseif (i==3) 
   dH=(-6*x^2+6*x)*l1;
elseif (i==4) 
   dH=(3*x^2-2*x);
end
