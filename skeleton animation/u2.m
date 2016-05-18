function u2
FV = stlread('T2_8.stl');

T=eye(4); T=rotationMatrix([0, 0, -90]*pi/180);
T(1:3,4)=[8 0 -20];
FV.vertices = [FV.vertices ones(size(FV.vertices, 1), 1)]*T';
FV.vertices = FV.vertices(:, 1:3);


master_l=[24 8];
master_r=[0 -10 0; 0 -25 0] * pi/180;
master_rconst=[0 0 0; 0 0 0] * pi/180;
master_rsin=[0 0 0; 0 0 15] * pi/180;
master_rphase=[0 0 0; 0 0 0] * pi/180;
master_rspeed=[1 1 1; 1 1 2];

for i=1:length(master_l)
  master_xyz{i}=[0 master_l(i); 0 0; 0 0; 1 1];
  master_v{i}=[];
end

slave_l{1}=[7 14 8];
slave_l{2}=[7 14 8];
slave_l{3}=[5 16 10];
slave_l{4}=[5 16 10];
slave_l{5}=[20];

slave_r{1}=[-110 0 90; 0 0 90; 25 -60 -25] * pi/180;
slave_r{2}=[110 0 -90; 0 0 -90; -25 -60 25] * pi/180;
slave_r{3}=[-100 0 90; 0 0 85; 0 10 10] * pi/180;
slave_r{4}=[100 0 -90; 0 0 -85; 0 10 -10] * pi/180;
slave_r{5}=[0 225 0] * pi/180;

slave_rconst{1}=[0 0 0; 0 0 0; 0 0 0] * pi/180;
slave_rconst{2}=[0 0 0; 0 0 0; 0 0 0] * pi/180;
slave_rconst{3}=[0 0 0; 0 0 0; 0 15 0] * pi/180;
slave_rconst{4}=[0 0 0; 0 0 0; 0 15 0] * pi/180;
slave_rconst{5}=[0 0 0] * pi/180;

slave_rsin{1}=[0 0 0; 0 15 0; 0 15 0] * pi/180;
slave_rsin{2}=[0 0 0; 0 -15 0; 0 -15 0] * pi/180;
slave_rsin{3}=[0 0 0; 0 -15 0; 0 30 0] * pi/180;
slave_rsin{4}=[0 0 0; 0 15 0; 0 -30 0] * pi/180;
slave_rsin{5}=[0 5 5] * pi/180;

slave_rphase{1}=[0 0 0; 0 0 0; 0 45 0] * pi/180;
slave_rphase{2}=[0 0 0; 0 0 0; 0 45 0] * pi/180;
slave_rphase{3}=[0 0 0; 0 0 0; 0 45 0] * pi/180;
slave_rphase{4}=[0 0 0; 0 0 0; 0 45 0] * pi/180;
slave_rphase{5}=[0 0 45] * pi/180;

slave_rspeed{1}=[1 1 1; 1 1 1; 1 1 1];
slave_rspeed{2}=[1 1 1; 1 1 1; 1 1 1];
slave_rspeed{3}=[1 1 1; 1 1 1; 1 1 1];
slave_rspeed{4}=[1 1 1; 1 1 1; 1 1 1];
slave_rspeed{5}=[1 2 1];

target0=[-25; 0; 20];
targetD=[0; 20; 10];
targetP=[0; 0; pi/4];

for i=1:5 
  for j=1:length(slave_l{i})
    slave_xyz{i}{j}=[0 slave_l{i}(j); 0 0; 0 0; 1 1];
    slave_v{i}{j}=[];
  end
end

slave_j=[1 1 2 2 1];

cla;
patch(FV, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.3); 
axis equal; grid on; hold on; xlabel('x');ylabel('y'); axis([-25 40 -20 20 -20 26]); %view([-1 1 1]);

rr=master_r;
[TT, SS, master_xyz1] = placeBones(rr, master_l, master_xyz);
for i=1:length(master_l)
  for j=find(slave_j==i)
    rr=slave_r{j};
    [~,~,slave_xyz1{j}]=placeBones(rr, slave_l{j}, slave_xyz{j}, TT{i}, SS{i}); 
  end
end

for i=1:size(FV.vertices, 1)
  mn=1e10; mnj=0; mnk=0;
  v=FV.vertices(i,:); v=v';
  for j=1:length(slave_xyz1)
    for k=1:length(slave_xyz1{j})
      p1=slave_xyz1{j}{k}(:,1);
      p2=slave_xyz1{j}{k}(:,2);
      d=distance(p1,p2, v);
      if d < mn
         mn = d; mnj=j; mnk=k; 
      end
    end
  end
  for k=1:length(master_xyz1)
    p1=master_xyz1{k}(:,1);
    p2=master_xyz1{k}(:,2);
    d=distance(p1,p2, v);
    if d < mn
       mn = d; mnj=0; mnk=k; 
    end
  end
  
  if mnj == 0
    master_v{mnk}=[master_v{mnk} i];  
  else
    slave_v{mnj}{mnk}=[slave_v{mnj}{mnk} i];
  end
end

rr=master_r;
for k=1:length(master_v)
  xyz{k}=FV.vertices(master_v{k},:);
  xyz{k}=[xyz{k}'; ones(1, size(xyz{k},1))];
end
[TT, SS, master_xyz1] = unplaceBones(rr, master_l, xyz);

for i=1:length(master_l)
  for j=find(slave_j==i)
    rr=slave_r{j};
    for k=1:length(slave_v{j})
      xyz{k}=FV.vertices(slave_v{j}{k},:);
      xyz{k}=[xyz{k}'; ones(1, size(xyz{k},1))];
    end
    [~,~,slave_xyz1{j}]=unplaceBones(rr, slave_l{j}, xyz, TT{i}, SS{i});

  end
end

v=zeros(size(FV.vertices));

for t=0:0.1:10000
  target=target0 + targetD .* sin(t*0.1 + targetP);
    
  ik=inverseKinematic(slave_l{5}, target);
  
  rr=master_r + master_rconst + master_rsin .* sin(master_rphase + master_rspeed*t);
  [TT, SS, xyz1] = placeBones(rr, master_l, master_xyz1);
  for k=1:length(master_v)
    xyz1{k}=xyz1{k}(1:3,:);
    v(master_v{k},:)=xyz1{k}';
  end
  
  for i=1:length(master_l)
    for j=find(slave_j==i)
      rr=slave_r{j} + slave_rconst{j} + slave_rsin{j} .* sin(slave_rphase{j} + slave_rspeed{j}*t);
      if j==5, rr=ik; end
      [~,~,xyz1]=placeBones(rr, slave_l{j}, slave_xyz1{j}, TT{i}, SS{i}); 
      for k=1:length(slave_v{j})
        xyz1{k}=xyz1{k}(1:3,:);
        v(slave_v{j}{k},:)=xyz1{k}';
      end
    end
  end
  
  
  patch('Faces',FV.faces,'Vertices',v, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.75); 
  plot3(target(1), target(2), target(3), 'r.','markersize',20);
  drawnow; pause(0.01);
  cla;
end

end

function d=distance(a, b, v)
  a=a(1:3); b=b(1:3); 
  l2=sum((a - b).^2);
  t = max(0, min(1, dot(v - a, b - a) / l2));
  projection = a + t * (b - a);
  d=sqrt(sum((v - projection).^2));
end



function [ TT, SS, xyz1 ]=placeBones(r, l, xyz, Tcum, S)
  if nargin < 5
    S=zeros(3,1);
  end
  if nargin < 4
    Tcum=eye(4);
  end
  
  for i=1:length(l)
    TT{i} = Tcum;
    SS{i} = S;
    T=eye(4);
    T=rotationMatrix(r(i,:));
    T(1:3,4)=S;
    Tcum=Tcum*T;
    A=Tcum*xyz{i};
    xyz1{i}=A;
    S=[l(i); 0; 0];
    
    A=Tcum*[0 l(i);0 0; 0 0; 1 1];
    c=[1 0 0];
    plot3(A(1,:), A(2,:), A(3,:),'linewidth',2,'color',c);
    A=Tcum*[0 0; 0 0; 0 2; 1 1];
    plot3(A(1,:), A(2,:), A(3,:),'linewidth',2,'color',c);
  end
end

function [ TT, SS, xyz1 ]=unplaceBones(r, l, xyz, Tcum, S)
  if nargin < 5
    S=zeros(3,1);
  end
  if nargin < 4
    Tcum=eye(4);
  end
  
  for i=1:length(l)
    TT{i} = Tcum;
    SS{i} = S;
    T=eye(4);
    T=rotationMatrix(r(i,:));
    T(1:3,4)=S;
    Tcum=Tcum*T;
    TcumI=inv(Tcum);
    A=TcumI*xyz{i};
    xyz1{i}=A;
    S=[l(i); 0; 0];
  end
end


function R=rotationMatrix(angles)
    
    c=cos(angles(1)); s=sin(angles(1));
    Rx = [1 0 0 0;
          0 c -s 0;
          0 s c 0;
          0 0 0 1];
    c=cos(angles(2)); s=sin(angles(2));
    Ry = [c 0 s 0;
          0 1 0 0;
          -s 0 c 0;
          0 0 0 1];
    c=cos(angles(3)); s=sin(angles(3));
    Rz = [c -s 0 0;
          s c 0 0;
          0 0 1 0;
          0 0 0 1];
    
    R=Rz*Ry*Rx;
return
end

function R=inverseKinematic(l, target)
  srigs=cell(length(l)+1,1);
  srigs{1}=[0;0;0];
  Wa=eye(length(l)); Wb=eye(length(l)); Wg=eye(length(l));
  k=30;
  
  for i=1:length(l), srigs{i+1}=[l(i);0;0]; end;
  angles=zeros(length(l), 3);
  
  step = 0.00001;
  
  for iii=1:50
    R=cell(length(l),1);
    for i=1:length(l)
      R{i}=rotationMatrix(angles(i,:)); 
      R{i}=R{i}(1:3,1:3);
    end
    P=[srigs{end};1]; for i=length(l):-1:1, P=[R{i}, srigs{i}; 0 0 0 1]*P; end
    PP=(P(1:3)-target);
      
    [DPDa,DPDb,DPDg]=DPDangles(R,angles,srigs);  
    DpsiDa=angles(:,1)'*Wa+k*PP'*DPDa;  
    DpsiDb=angles(:,2)'*Wb+k*PP'*DPDb;
    DpsiDg=angles(:,3)'*Wg+k*PP'*DPDg;
    
    angles = angles - step*[DpsiDa DpsiDb DpsiDg];
    angles(:,2)=min(-110*pi/180, max(angles(:,2), -190*pi/180));
    
  end
  R=angles;
end

function [DPDa,DPDb,DPDg]=DPDangles(R,Kampai,srigs)

nrigs=length(srigs)-1; % kaulu skaicius

DPDg=zeros(4,nrigs); DPDb=zeros(4,nrigs); DPDa=zeros(4,nrigs);

for iii=1:nrigs
  kampai=Kampai(iii,:);sa=sin(kampai);ca=cos(kampai);
  DPDg(:,iii)=[srigs{end};1];
  for i=nrigs:-1:1
    if  i ~= iii,   DPDg(:,iii)=[R{i},srigs{i};0 0 0 1]*DPDg(:,iii);
    else,     DR=[-sa(3)  -ca(3)   0   ;  ca(3)  -sa(3)    0   ;   0     0     0  ]*...
                 [ ca(2)   0     sa(2) ;  0        1       0   ; -sa(2)  0    ca(2)]*...
                 [    1    0       0   ;  0      ca(1)  -sa(1) ;  0    sa(1)  ca(1)];
              DPDg(:,iii)=[DR,zeros(3,1);0 0 0 0]*DPDg(:,iii);
    end
  end
      
  DPDb(:,iii)=[srigs{end};1];
  for i=nrigs:-1:1
    if  i~=iii,   DPDb(:,iii)=[R{i},srigs{i};0 0 0 1]*DPDb(:,iii);
    else,     DR=[  ca(3)  -sa(3)   0   ;  sa(3)   ca(3)     0   ;    0     0     1    ]*...
                 [ -sa(2)   0     ca(2) ;   0        0       0   ; -ca(2)  0    -sa(2) ]*...
                 [    1    0       0    ;   0      ca(1)  -sa(1) ;    0    sa(1)  ca(1)];
              DPDb(:,iii)=[DR,zeros(3,1);0 0 0 0]*DPDb(:,iii);
    end
  end
   
  DPDa(:,iii)=[srigs{end};1];
  for i=nrigs:-1:1
    if  i~=iii,   DPDa(:,iii)=[R{i},srigs{i};0 0 0 1]*DPDa(:,iii);
    else,     DR=[ca(3)  -sa(3)   0    ;  sa(3)  ca(3)     0    ;   0     0     1    ]*...
                 [ ca(2)   0     sa(2) ;  0        1       0    ; -sa(2)  0    ca(2) ]*...
                 [    0    0       0   ;  0      -sa(1)  -ca(1) ;  0    ca(1)  -sa(1)];
              DPDa(:,iii)=[DR,zeros(3,1);0 0 0 0]*DPDa(:,iii);
    end
  end
        
end
DPDa(4,:)=[];  DPDb(4,:)=[]; DPDg(4,:)=[]; % panaikinama "homogenine" vektoriu komponente

return
end

