function u1(d)
if (nargin < 1)
  d=0.5; 
end

[lambda1, lambda2, alpha1, alpha2, b1, b2, b3, t1, t2, q, h, Tinf] = parameters;
c=42;
ro=100;

[V,G,lambda,alpha,b,IS,Tc,Q]=buildGrid(lambda1, lambda2, alpha1, alpha2, b1, b2, b3,t1,t2,q,d);
Tc(Tc==0)=273*ones(size(Tc(Tc==0)));
n=size(V,1);
m=size(G,1);

figure; hold on; axis equal; axis([0 1 0 1]*8);

K=zeros(n,n);
C=zeros(n,n);
S=zeros(n,1);
P=zeros(n,1);

GG = [G G(:,1)];
for i=1:m
  VV=V(G(i, :), :);
  [Ke,Ce,Se,Pe]=elementMatrices(VV,lambda(i),alpha(i),Tinf,b(i),h,c,ro);
  Gi=G(i, :);
  K(Gi,Gi)=K(Gi,Gi)+Ke;
  C(Gi,Gi)=C(Gi,Gi)+Ce;
  S(Gi)=S(Gi)+Se;
  P(Gi)=P(Gi)+Pe;
end

dt=1;
b0=dt^2/4; b1=dt/2; b2=1;
X=Tc; DX=zeros(size(S)); DDX=zeros(size(S));
for t=dt:dt:1000
    
    q0=X(~IS)+dt*DX(~IS)+dt^2/2*DDX(~IS); q1=DX(~IS)+dt*DDX(~IS); q2=DDX(~IS);
    
    F=(S+P)+Q;
    
    M=zeros(n,n);
    Kis=K(~IS,~IS); Cis=C(~IS,~IS); Mis=M(~IS,~IS); Fis=F(~IS);
    
    deltaDDX=zeros(n,1);
    deltaDDX(~IS)=(b2*Mis+b1*Cis+b0*Kis)\(Fis-K(~IS,IS)*Tc(IS)-(Mis*q2+Cis*q1+Kis*q0));
    X(~IS)=q0+b0*deltaDDX(~IS); DX(~IS)=q1+b1*deltaDDX(~IS); DDX(~IS)=q2+b2*deltaDDX(~IS);

    cla; hold on;
    for i=1:m
      VV=V(GG(i, :), :);
      XX=X(GG(i, :));
      fill(VV(:,1), VV(:, 2), XX);
    end
    colorbar; caxis([0, 1000]); drawnow;
    pause(0.01);
    if (sum(sum(abs(DX))) < 1), break; end
end
end
