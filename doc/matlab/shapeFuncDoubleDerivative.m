function ddH=shapeFuncDoubleDerivative(x,i,l)
% x  !< location along the element
% i  !< row number
% l !< length of the element
% ddH  !< double derivative
%      !< value of shape function at the given location

l1=1./l;
l2=l1*l1;

if (i == 1) 
   ddH=(12*x-6)*l2;
elseif (i == 2) 
   ddH=(6*x-4)*l1;
elseif (i == 3) 
   ddH=(-12*x+6)*l2;
elseif (i == 4) 
   ddH=(6*x-2)*l1;
end
