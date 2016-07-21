function [Ke, Ce, Se,Pe]=elementMatrices(V,lambda,alpha,Tinf,b,h,c,ro)
if nargin < 6
  error('input_example :  a is a required input')
elseif nargin < 8
  c=0; ro=0; 
end
x1=V(1,1); x2=V(2,1); x3=V(3,1);
y1=V(1,2); y2=V(2,2); y3=V(3,2);

A=abs((x1*y2+x2*y3+x3*y1-x3*y2-x2*y1-x1*y3)/2);

D=[lambda, 0; 0, lambda];
B=matrixB(x1,y1,x2,y2,x3,y3);

Ke=A*h*B'*D*B+2*alpha*A/12*[4,0,0; 0,4,0; 0,0,4];
Ce=c*ro*h*A/12*[2 1 1; 1 2 1; 1 1 2];
Se=2*alpha*A*Tinf/3*[1;1;1];
Pe=h*A*b/3*[1;1;1];
end

function B=matrixB(x1,y1,x2,y2,x3,y3)
A=(x1*y2+x2*y3+x3*y1-x3*y2-x2*y1-x1*y3)/2;

b1=(y2-y3)/(2*A);
b2=(y3-y1)/(2*A);
b3=(y1-y2)/(2*A);
c1=(x3-x2)/(2*A);
c2=(x1-x3)/(2*A);
c3=(x2-x1)/(2*A);

B=[b1,b2,b3; c1,c2,c3];
end