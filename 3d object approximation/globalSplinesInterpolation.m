function ZZ=globalSplinesInterpolation(Z, y, x, divy, divx)
nxx=divx*(length(x)-1);nyy=divy*(length(y)-1); 

%Dvigubu Lagranzo daugianariu apskaiciavimas ir sumos kaupimas:
Zr=reshape(Z,length(x)*length(y),1);
ZZ=zeros(1,nxx*nyy); k=0;
for i=1:length(x), [GX,~]=gsBaseFunction(x,i,divx,1);
    for j=1:length(y), [GY,~]=gsBaseFunction(y,j,divy,0); 
        k=k+1; ZZ=ZZ+reshape(GY'*GX,1,nxx*nyy)*Zr(k);
    end
end
ZZ=reshape(ZZ,nyy,nxx);
end
%*****************************************************************

function [GS,sss]=gsBaseFunction(X,i,nnn,iopt)
nP=length(X);
GS=[];sss=[];
    Y=zeros(1,nP); Y(i)=1; DDF=splineCoefficients(X,Y,iopt);
    for iii=1:nP-1  %------  ciklas per intervalus tarp gretimu tasku
        [GS1,sss1]=spline(X(iii:iii+1),Y(iii:iii+1),DDF(iii:iii+1),nnn); % reiksmes ir vaizdavinmo abscises
%         plot(sss1,GS1)
%         pause
        GS=[GS,GS1];sss=[sss,sss1];
    end %-----------------ciklas per intervalus pabaiga
return
end

function DDF=splineCoefficients(X,Y,iopt)
% apskaiciuojamos antros isvestines splaino mazguose
% iopt=1 - periodinis splainas

n=length(X);
A=zeros(n);b=zeros(n,1);
d=X(2:n)-X(1:(n-1));
 for i=1:n-2
     A(i,i:i+2)=[d(i)/6, (d(i)+d(i+1))/3,d(i+1)/6];
     b(i)=(Y(i+2)-Y(i+1))/d(i+1)-(Y(i+1)-Y(i))/d(i);
 end
 
if iopt == 0,  A(n-1,1)=1;A(n,n)=1;
else, A(n-1,[1,2,n-1,n])=[d(1)/3, d(1)/6, d(n-1)/6,d(n-1)/3]; 
      A(n,[1,n])=[1,-1];  
      b(n-1)=(Y(2)-Y(1))/d(1)-(Y(n)-Y(n-1))/d(n-1);
end

DDF=A\b;
 
return
end

function [S,sss]=spline(X,Y,DDF,nnn)
% splaino intervale tarp dvieju tasku apskaiciavimas
% nnn - vaizdavio tzku skaicius
% S - splaino reiksmes
% sss - vaizdavimo abscises
d=X(2)-X(1);
sss=X(1):(X(2)-X(1))/(nnn-1):X(2);
S=DDF(1)/2*(sss-X(1)).^2+(DDF(2)-DDF(1))/(6*d)*(sss-X(1)).^3+(sss-X(1))*((Y(2)-Y(1))/d-DDF(1)*d/3-DDF(2)*d/6) +Y(1);

return
end