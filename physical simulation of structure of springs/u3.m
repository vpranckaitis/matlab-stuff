function pos=u3
m1=1; m2=2;
k1=100; k2=120;
% m1=1; m2=2; k1=12000; k2=10;

dU = [-2 -1];

% konstrukcijos parametrai
mass=[m2 m2 m1 m2 m1];        % mases
rad=[0.1 0.1 0.1 0.1 0.1] * 2; % spinduliai
ind=[1 2; 1 3; 1 4; 1 5; 2 5; 3 4; 4 5];   % elementu mazgu globalieji numeriai   
cor=[3 4; 5 3; 1 2; 3 1; 5 1]; % mazgu koordinates
g=9.8;                      % laisvo kritimo pagreitis
F=zeros(length(mass)*2,1); % pridetos jegos
F(2:2:end)=F(2:2:end)-g*mass';
IS= [1 1 0 1 1 0 0 0 0 0];        % sablonas duotiems greiciams
stiff=[k1 k1 k2 k1 k2 k1 k2];        % standumai


llsk=2; nmz=length(mass); NN=nmz*llsk;
U=zeros(NN,1); U([1 2])=U([1 2])+dU';

h=figure; axis equal;axis ([-1 6 -1 5]);grid on; set(gcf,'Name','Niutono-Rafsono metodas');

for it=1:300
    psi=netiktis(U,cor,ind,stiff,F);
    KT=jakobian(U,cor,ind,stiff);
    
    deltaU=zeros(2*nmz,1);
    deltaU(~IS)=KT(~IS,~IS) \ psi(~IS);
    
    U=U+deltaU;

    figure(h); cla; hold on; vaizdavimas(U,ind,cor,rad,F,IS); hold off; drawnow; pause(1);
    
    eps=norm(deltaU)/(norm(U) + norm(deltaU));
    
    if eps<1e-8, break,end
end
pos=cor+reshape(U, size(cor'))';

M=diag(mass(floor((2:1:nmz*2+1)/2)));
[Y,om]=eig(KT(~IS,~IS),M(~IS,~IS));
om=diag(om).^(0.5);
om, relativeOm=max(om)/min(om)
end

function psi=netiktis(U,cor,ind,k,F)
nel=size(ind,1); nmz=size(cor,1);
psi=zeros(2*nmz,1);

for i=1:nel
    ri=ind(i,1); si=ind(i,2);
    r=[2*ri-1, 2*ri]; s=[2*si-1, 2*si];
   
    cr=cor(ri,:); cs=cor(si,:);
    Lvec=cs+U(s)'-cr-U(r)'; L=norm(Lvec); n=Lvec/L; L0=norm(cs-cr);
    T=k(i)*(L-L0);
    psi([r,s])=psi([r,s])+[n';-n']*T;
end
psi=psi+F;
end

function KT=jakobian(U,cor,ind,k)
nel=size(ind,1); nmz=size(cor,1);
KT=zeros(2*nmz,2*nmz);
for i=1:nel
    ri=ind(i,1); si=ind(i,2);
    r=[2*ri-1, 2*ri]; s=[2*si-1, 2*si];
   
    cr=cor(ri,:); cs=cor(si,:);
    Lvec=cs+U(s)'-cr-U(r)'; L=norm(Lvec); n=Lvec/L; L0=norm(cs-cr);
    
    z1=1-L0/L+Lvec(1)^2*L0/L^3;
    z23=Lvec(1)*Lvec(2)*L0/L^3;
    z4=1-L0/L+Lvec(2)^2*L0/L^3;
    
    KTe=k(i)*[z1 z23; z23 z4];
    
    KT([r,s],[r,s])=KT([r,s],[r,s])+[ KTe -KTe;
                                     -KTe  KTe];
end
end

function vaizdavimas(U,ind,cor,rad,F,IS)
nmz=length(rad);llsk=2;nel=size(ind,1);

xlim=get(gca,'XLim'); ylim=get(gca,'YLim');
xn=xlim(2)-xlim(1);yn=ylim(2)-ylim(1);
range=min(xn,yn);
maxForce=max(abs(F));
mast= range/maxForce*0.1;
constrLength=range/17;      

for i=1:nmz
    r=[(i-1)*llsk+1,i*llsk]; c=cor(i,:)'+U(r);
    rd=rad(i);
    rectangle('Position',[c'-rd,2*rd,2*rd],'Curvature',[1,1],'FaceColor',[0.4 0.6 1]);
        
    f=F(r)*mast; 
    x1=c(1);x2=c(1)+f(1);y1=c(2);y2=c(2)+f(2); 
    line([x1,x2],[y1,y2],'Color','red','LineWidth',1);
    varr=[x1-x2;y1-y2]; varr =varr/norm(varr)*range/40; 
    alf=pi/6; transf = [cos(alf) sin(alf);-sin(alf) cos(alf)];
    varr1=transf*varr; line([x2, x2+varr1(1)],[y2, y2+varr1(2)],'Color','red','LineWidth',1);
    varr1=transf'*varr;line([x2, x2+varr1(1)],[y2, y2+varr1(2)],'Color','red','LineWidth',1);
                       
    constr=IS(r);
    if constr(1) ~= 0, line(([c(1), c(1)]),([c(2)-constrLength/2, c(2)+constrLength/2]),'Color',[ 0.2 0.2 0.2],'LineWidth',3);end  
    if constr(2) ~= 0, line(([c(1)-constrLength/2, c(1)+constrLength/2]),([c(2), c(2)]),'Color',[ 0.2 0.2 0.2],'LineWidth',3);end  
end

for i=1:nel
    ri=ind(i,1);si=ind(i,2); 
    r=[(ri-1)*llsk+1,ri*llsk];s=[(si-1)*llsk+1,si*llsk];
    cr=cor(ri,:)'+U(r); cs=cor(si,:)'+U(s);
    plot([cr(1),cs(1)] , [cr(2),cs(2)],'b-');
end

return
end
