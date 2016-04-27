function pos=u1(dt)
m1=1; m2=2;
k1=100; k2=120;
c1=20; c2=10;
% m1=1; m2=2; k1=12000; k2=10;

dU = [-2 -1];
tu0 = 0.2; dtu = 0.9;

mass=[m2 m2 m1 m2 m1];        % mases
rad=[0.1 0.1 0.1 0.1 0.1] * 2; % spinduliai
ind=[1 2; 1 3; 1 4; 1 5; 2 5; 3 4; 4 5];   % elementu mazgu globalieji numeriai   
cor=[3 4; 5 3; 1 2; 3 1; 5 1]; % mazgu koordinates
g=9.8;                      % laisvo kritimo pagreitis
F=-[(0*mass); mass]*g; % pridetos jegos
F=F(:);
IS= [1 1 0 1 1 0 0 0 0 0];        % sablonas duotiems greiciams

stiff=[k1 k1 k2 k1 k2 k1 k2];        % standumai
damp_rel=[c1 c1 c2 c1 c2 c1 c2];   % reliatyviojo mazgu judesio slopintuvai

llsk=2; nmz=length(mass); NN=nmz*llsk; 
U=zeros(NN,1); DU=zeros(NN,1);   

if ~exist('dt','var'), dt=0.01; end; TT=10;   
Trez=[0]; Urez=U; DUrez=DU; 

h=figure; set(gcf,'Name',sprintf('Centriniø skirtumø metodas, dt=%0.3f', dt));
h1=subplot(2,2,[1 3]); axis equal; axis([-1 6 -1 5]);grid on; title('Sistemos iðsidëstymas');
h2=subplot(2,2,2); hold on; axis([0 TT -2 1]); title('Mazgø poslinkiai'); 
h3=subplot(2,2,4); hold on; axis([0 TT -5 2]); title('Mazgø greièiai');

for t=dt:dt:TT
    DDU=pagreiciai(U,DU,t,mass,stiff,damp_rel,F,ind,cor);
    DU=DU+dt*DDU'; DU(find(IS))=0; U=U+dt*DU; 
    if (t > tu0 + dtu), 
        U([1 2])=dU'; DU([1 2])=[0 0]; DDU([1 2])=[0 0];
    elseif (t > tu0), 
        omega=(pi/2)/dtu;
        U([1 2])=dU' * sin(omega*(t - tu0));
        DU([1 2])=dU' * omega * cos(omega*(t - tu0));
        DDU([1 2])=-dU' * omega^2 * cos(omega*(t - tu0));
    end
    
    Trez(end+1)=t; Urez(:,end+1)=U; DUrez(:,end+1)=DU;
    
    figure(h);subplot(h2); cla; plot(Trez,Urez); subplot(h3); cla; plot(Trez,DUrez);
    subplot(h1); cla; hold on; vaizdavimas(U,ind,cor,rad,F,IS); hold off; drawnow; pause(dt);
    if (sum(sum(abs(DU))) < 1e-6), break; end
end
pos=cor+reshape(U, size(cor'))';
return
end

function DDU=pagreiciai(U,DU,t,mass,stiff,damp_rel,F,ind,cor)
llsk=2; nmz=length(mass); nel=size(ind,1);
 
T=F;
for i=1:nel
    ri=ind(i,1); si=ind(i,2);
    r=[(ri-1)*llsk+1,ri*llsk];s=[(si-1)*llsk+1,si*llsk];
    cr=cor(ri,:)'+U(r); cs=cor(si,:)'+U(s);
    dur=DU(r);dus=DU(s);
    
    l0=norm(cor(si,:)-cor(ri,:)); lrs=norm(cs-cr); n=(cs-cr)/lrs;
    
    Trs=stiff(i)*(lrs-l0)+damp_rel(i)*dot( n, dus-dur);
     
    T(r)=T(r)+Trs*n; T(s)=T(s)-Trs*n;
end

for i=1:nmz 
    r=[(i-1)*llsk+1,i*llsk];
    DDU(r)=T(r)/mass(i);
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