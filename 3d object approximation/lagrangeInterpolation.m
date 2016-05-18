function ZZ=lagrangeInterpolation(Z, y, x, yy, xx)
nxx=length(xx); nyy=length(yy);
ZZ=zeros(1,nxx*nyy); k=0;
for i=1:length(x) 
    LX=lagrangePolynomial(x,i,xx);
    for j=1:length(y) 
        LY=lagrangePolynomial(y,j,yy); 
        k=k+1; 
        ZZ=ZZ+reshape(LY'*LX,1,nxx*nyy)*Z(k);
    end
end
ZZ=reshape(ZZ,nyy,nxx);
end

function L=lagrangePolynomial(X,j,x)
n=length(X);
L=1;
for k=1:n, if k ~= j, L=L.*(x-X(k))/(X(j)-X(k)); end, end
return
end