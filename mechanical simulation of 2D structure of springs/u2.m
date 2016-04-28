function pos=u2(dt)
m1=1; m2=2;
k1=100; k2=120;
c1=20; c2=10;
% m1=1; m2=2; k1=12000; k2=10;

dU = [-2 -1];
tu0 = 0.2; dtu = 0.9;

% konstrukcijos parametrai
mass=[m2 m2 m1 m2 m1];        % mases
rad=[0.1 0.1 0.1 0.1 0.1] * 2; % spinduliai
ind=[1 2; 1 3; 1 4; 1 5; 2 5; 3 4; 4 5];   % elementu mazgu globalieji numeriai   
cor=[3 4; 5 3; 1 2; 3 1; 5 1]; % mazgu koordinates
g=9.8;                      % laisvo kritimo pagreitis
F=zeros(length(mass)*2,1); % pridetos jegos
F(2:2:end)=F(2:2:end)-g*mass';
IS=[1 1 0 1 1 0 0 0 0 0];

stiff=[k1 k1 k2 k1 k2 k1 k2];        % standumai
damp_rel=[c1 c1 c2 c1 c2 c1 c2];

llsk=2; nmz=length(mass); NN=nmz*llsk;

U=zeros(NN,1); DU=zeros(NN,1); DDU=zeros(NN,1);

M=diag(mass(floor((2:1:nmz*2+1)/2))); Mis=M(~IS,~IS);

if ~exist('dt','var'), dt=0.01; end; TT=10;
Trez=[0]; Urez=U; DUrez=DU; 

h=figure; set(gcf,'Name',sprintf('Niumarko metodas, dt=%0.3f', dt));
h1=subplot(2,2,[1 3]); axis equal; axis([-1 6 -1 5]);grid on; title('Sistemos iðsidëstymas');
h2=subplot(2,2,2); hold on; axis([0 TT -2 1]); title('Mazgø poslinkiai'); 
h3=subplot(2,2,4); hold on; axis([0 TT -5 2]); title('Mazgø greièiai');

b0=dt^2/4; b1=dt/2; b2=1;

filename = 'animation2.gif';

for t=dt:dt:TT
    if (t > tu0 + dtu), 
        U([1 2]) = dU';
        DU([1 2]) = [0 0];
        DDU([1 2]) = [0 0];
    elseif (t > tu0), 
        omega=(pi/2)/dtu;
        U([1 2]) = dU' * sin(omega*(t - tu0));
        DU([1 2]) = dU' * omega * cos(omega*(t - tu0));
        DDU([1 2]) = -dU' * omega^2 * cos(omega*(t - tu0));
    end
    
    q0=U(~IS)+dt*DU(~IS)+dt^2/2*DDU(~IS); q1=DU(~IS)+dt*DDU(~IS); q2=DDU(~IS);
    
    deltaDDU=zeros(2*nmz,1);
    for it=1:300
        U(~IS)=q0+b0*deltaDDU(~IS); DU(~IS)=q1+b1*deltaDDU(~IS); DDU(~IS)=q2+b2*deltaDDU(~IS);
        
        KT=jakobian(U,cor,ind,stiff);
        TU=stiffnessVector(U,cor,ind,stiff);
        CU=dampingVector(U,cor,ind,damp_rel,DU);
        C=dampingMatrix(U,cor,ind,damp_rel);
        
        KTis=KT(~IS,~IS); TUis=TU(~IS); CUis=CU(~IS); Cis=C(~IS,~IS);
        % Cis=Mis*(-4); CUis=Cis*DU(~IS);
        deltadeltaDDU=zeros(2*nmz,1);
        deltadeltaDDU(~IS)=(b2*Mis+b1*Cis+b0*KTis)\(CUis+TUis+F(~IS)-(Mis*DDU(~IS)));
        deltaDDU=deltaDDU+deltadeltaDDU;
        if (norm(deltadeltaDDU)<1e-6), break; end
    end
    
    Trez(end+1)=t; Urez(:,end+1)=U; DUrez(:,end+1)=DU;

    figure(h);subplot(h2); cla; plot(Trez,Urez); subplot(h3); cla; plot(Trez,DUrez);
    subplot(h1); cla; hold on; vaizdavimas(U,ind,cor,rad,F,IS); hold off; drawnow; pause(dt);
    if (sum(sum(abs(DU))) < 1e-6), break; end
end
pos=cor+reshape(U, size(cor'))';

subplot(2,2,2); plot(Trez,Urez); subplot(2,2,4); plot(Trez,DUrez);
end

function C=dampingMatrix(U,cor,ind,c)
nel=size(ind,1); nmz=size(cor,1);
C=zeros(2*nmz,2*nmz);
for i=1:nel
    ri=ind(i,1); si=ind(i,2);
    r=[2*ri-1, 2*ri]; s=[2*si-1, 2*si];
    cr=cor(ri,:); cs=cor(si,:);
    
    Lvec=cs+U(s)'-cr-U(r)'; L=norm(Lvec); n=Lvec/L;
    
    a=atan2(n(2), n(1));
    ci=c(i);
    Kel=[ci 0 -ci 0;
         0 0  0 0;
        -ci 0  ci 0;
         0 0  0 0];
    T=zeros(4,4);
    T(1:2,1:2)=[cos(a), sin(a);
                -sin(a), cos(a)];
    T(3:4,3:4)=[cos(a), sin(a);
                -sin(a), cos(a)];
    Ce=T'*Kel*T;
    
    C([r,s],[r,s])=C([r,s],[r,s])+Ce;
end
end

function TU=stiffnessVector(U,cor,ind,k)
nel=size(ind,1); nmz=size(cor,1);
TU=zeros(2*nmz,1);
for i=1:nel
    ri=ind(i,1); si=ind(i,2);
    r=[2*ri-1, 2*ri]; s=[2*si-1, 2*si];
   
    cr=cor(ri,:); cs=cor(si,:);
    Lvec=cs+U(s)'-cr-U(r)'; L=norm(Lvec); n=Lvec/L; L0=norm(cs-cr);
    TUe=k(i)*(L - L0);

    TU([r,s])=TU([r,s])+[n';-n']*TUe;
end
end

function CU=dampingVector(U,cor,ind,c,DU)
nel=size(ind,1); nmz=size(cor,1);
CU=zeros(2*nmz,1);
for i=1:nel
    ri=ind(i,1); si=ind(i,2);
    r=[2*ri-1, 2*ri]; s=[2*si-1, 2*si];
   
    cr=cor(ri,:); cs=cor(si,:);
    dur=DU(r); dus=DU(s);
    Lvec=cs+U(s)'-cr-U(r)'; L=norm(Lvec); n=Lvec/L;
    
    CUe=c(i)*dot(n,dus-dur);

    CU([r,s])=CU([r,s])+[n';-n']*CUe;
end
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
    
    KTeblokas=k(i)*[z1 z23; z23 z4];
    
    KT([r,s],[r,s])=KT([r,s],[r,s])+[ KTeblokas -KTeblokas;
                                     -KTeblokas  KTeblokas];
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
