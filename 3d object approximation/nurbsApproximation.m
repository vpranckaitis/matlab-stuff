function [XN, YN, ZN]=nurbsApproximation(X, Y, Z, W, divy, divx)

X=[X X(:,2:4)]; Y=[Y Y(:,2:4)]; Z=[Z Z(:,2:4)]; W=[W W(:,2:4)];

nx=size(X,2);ny=size(X,1);

% vaizdavimo(smulkus) etaloninis tinklelis
xxb=[0:1/divx:1];
yyb=[0:1/divy:1];
nxx=length(xxb);nyy=length(yyb);

%Bernsteino daugianariu masyvai
BX=Bernsteino_UBS_daugianaris(xxb);
BY=Bernsteino_UBS_daugianaris(yyb); 

XN=zeros((ny-3)*divy+1, (nx-3)*divx+1);
YN=zeros((ny-3)*divy+1, (nx-3)*divx+1);
ZN=zeros((ny-3)*divy+1, (nx-3)*divx+1);

for iii=0:nx-4, for jjj=0:ny-4

    Xr=reshape(X(jjj+1:jjj+4,iii+1:iii+4),16,1); Yr=reshape(Y(jjj+1:jjj+4,iii+1:iii+4),16,1); 
    Zr=reshape(Z(jjj+1:jjj+4,iii+1:iii+4),16,1); Wr=reshape(W(jjj+1:jjj+4,iii+1:iii+4),16,1);
    XX=zeros(1,nxx*nyy);   YY=zeros(1,nxx*nyy);  ZZ=zeros(1,nxx*nyy); WW=zeros(1,nxx*nyy);
    k=0;
    for i=1:4, for j=1:4
        k=k+1; 
        XX=XX+reshape(BY(j,:)'*BX(i,:),1,nxx*nyy)*Xr(k)*Wr(k);
        YY=YY+reshape(BY(j,:)'*BX(i,:),1,nxx*nyy)*Yr(k)*Wr(k);
        ZZ=ZZ+reshape(BY(j,:)'*BX(i,:),1,nxx*nyy)*Zr(k)*Wr(k);
        WW=WW+reshape(BY(j,:)'*BX(i,:),1,nxx*nyy)*Wr(k);
    end, end
    xFrom=iii*divx+1;
    xTo=xFrom+divx;
    yFrom=jjj*divy+1;
    yTo=yFrom+divy;
    XN(yFrom:yTo,xFrom:xTo)=reshape(XX./WW,nyy,nxx);
    YN(yFrom:yTo,xFrom:xTo)=reshape(YY./WW,nyy,nxx);
    ZN(yFrom:yTo,xFrom:xTo)=reshape(ZZ./WW,nyy,nxx);

end, end
end

%*****************************************************************
function B=Bernsteino_UBS_daugianaris(x)
    % n - valdanciuju tasku skaicius
    % x - vaizdavimo taskai parametro reiksmiu intervale [0,1]
      T=[x'.^3, x'.^2, x',x'.^0];
        CC=[ -1,  3, -3, 1;
              3, -6,  3, 0;
             -3,  0,  3, 0;
              1,  4,  1, 0]/6;
    B=(T*CC)';
end