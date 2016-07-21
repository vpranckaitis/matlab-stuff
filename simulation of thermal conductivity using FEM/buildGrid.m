function [V,G,lambda,alpha,b,IS,Tc,Q]=buildGrid(lambda1, lambda2, alpha1, alpha2, b1, b2, b3,t1,t2,q,d)
x=1:d:4; y=1:d:7;
[Y,X]=meshgrid(y, x);
V=[X(:), Y(:)];


Tc=zeros(size(V,1),1); Q=zeros(size(Tc)); IS=false(size(Tc));
[~,indx]=ismember(V, [4 1],'rows'); Tc(indx==1)=t1; IS(indx==1)=true;
[~,indx]=ismember(V, [1 7],'rows'); Tc(indx==1)=t2; IS(indx==1)=true;
[~,indx]=ismember(V, [4 4],'rows'); Q(indx==1)=q;

G = zeros(sum(size(X) - [1 1]),3);
lambda = zeros(size(G,1), 1);
alpha = zeros(size(lambda));
b = zeros(size(lambda));
k=1;
for i=1:(length(y)-1)
  for j=1:(length(x)-1)
    nx=length(x);
    if y(i) < 4
      t1=(i-1)*nx+j; t2=t1+1;
      t4=i*nx+j; t3=t4+1;
      G(k,:)=[t1 t2 t4];  G(k+1,:)=[t2 t3 t4];
      lambda(k)=lambda1; lambda(k+1)=lambda1;
      alpha(k)=alpha1; alpha(k+1)=alpha1;
      if y(i) + x(j) < 5
          b(k)=b1;
      else
          b(k)=b2;
      end;
      if y(i+1) + x(j+1) <= 5
          b(k+1)=b1;
      else
          b(k+1)=b2;
      end;
    else
      t1=(i-1)*nx+j; t2=t1+1;
      t4=i*nx+j; t3=t4+1;
      G(k,:)=[t1 t2 t3];  G(k+1,:)=[t3 t4 t1];
      if y(i) - x(j+1) < 3
          lambda(k)=lambda2;
          alpha(k)=alpha2;
          b(k)=b1;
      else
          lambda(k)=lambda1;
          alpha(k)=alpha1;
          b(k)=b3;
      end;
      if y(i+1) - x(j) <= 3
          lambda(k+1)=lambda2;
          alpha(k+1)=alpha2;
          b(k+1)=b1;
      else
          lambda(k+1)=lambda1;
          alpha(k+1)=alpha1;
          b(k+1)=b3;
      end;
    end
    k=k+2;
  end
end
end