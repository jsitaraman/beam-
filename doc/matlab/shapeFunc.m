function H=shapeFunc(x,i,l)
% x  !< location along the element
% i  !< row number
% l !< length of the element
% H  !< value of shape function at the given location
if (i==1) 
 H=2*x^3-3*x^2+1;
elseif (i==2) 
 H=(x^3-2*x^2+x)*l;
elseif (i==3) 
 H=(-2*x^3+3*x^2);
elseif (i==4)
 H=(x^3-x^2)*l;
end

