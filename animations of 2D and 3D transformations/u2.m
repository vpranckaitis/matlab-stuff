function u2
clf;
%
omg1 = 10;
omg2 = 30;
S = [0.5 2];
Z = [0 5];
tMax = 2;
d = 15;
B = [ 0  0;
      10  0;
      10  0];
%  

points = [-1  1  0 -1 1  0;
           1  1 -1  1 1 -1;
          -1 -1 -1  1 1  1;
           1  1  1  1 1  1];
          
f = [1 2 3 1;
     1 2 5 4;
     1 3 6 4;
     2 3 6 5;
     4 5 6 4]; % faces
      
subplot(1, 2, 1); grid on, hold on, axis equal,axis(7*[-1 1 -1 1 -1 1]),view([1 -1 1]),xlabel('x'),ylabel('y'),zlabel('z'); 
plt_coord();
subplot(1, 2, 2); grid on, hold on, axis equal,axis(7*[-1 1 -1 1]);

v = [2; 0; Z(2); 1];
dt = 0.01;

w = [];

for t=0:dt:100
  angl = omg1 * t;
  
  vT = normalizeCoords(rotationZ(angl) * v);
  
  c = vT(1:3) / vT(3) * min(Z(1) + t / tMax * (Z(2) - Z(1)), Z(2));
    
  ra = cross(vT(1:3)', [0 0 1]); % pasukimo aðis
  r = acos(dot(vT(1:3)', [0 0 1]) / norm(vT(1:3)')); % pasukimo kampas
  
  rz = omg2 * t; % savasis prizmës sukimasis
  
  s1 = min(S(1) + t / tMax * (S(2) - S(1)), S(2));
  s2 = max(S(2) - t / tMax * (S(2) - S(1)), S(1));
  
  pointsT = normalizeCoords(translation(c(1), c(2), c(3)) * rotation(r, ra) * rotationZ(rz) * scale(s1, s1, s2) * points);
  
  x = reshape(pointsT(1, f), size(f))';
  y = reshape(pointsT(2, f), size(f))';
  z = reshape(pointsT(3, f), size(f))';
  
  b = B(:,2)-B(:,1);
  ra = cross(b', [0 0 -1]); % pasukimo aðis
  r = acos(dot(b, [0 0 -1]) / norm(b)); %pasukimo kampas
  
  Tp = perspective(d) * rotation(-r, ra) * translation(-B(1,1), -B(2,1), -B(3,1));
  vTp = normalizeCoords(Tp * vT);
  
  pointsTp = normalizeCoords(Tp * pointsT);
  
  xp = reshape(pointsTp(1, f), size(f))';
  yp = reshape(pointsTp(2, f), size(f))';
  
  delete(w);
  subplot(1, 2, 1);
  w = [quiver3(0, 0, 0, vT(1), vT(2), vT(3), 'b', 'linewidth', 3);
       quiver3(B(1,1), B(2,1), B(3,1), b(1), b(2), b(3), 'r', 'linewidth', 3);
       fill3(x, y, z, 1:5)];
  
  subplot(1, 2, 2);
  w = [w;
       quiver(0, 0, vTp(1), vTp(2), 'b', 'linewidth', 3);
       fill(xp, yp, 'r')];
  drawnow;
  pause(dt);
end
end

function M=normalizeCoords(L) 
M = L ./ L(4*ones(1,4), :);
end

function T=scale(sx, sy, sz) 
T = [sx  0  0 0;
      0 sy  0 0;
      0  0 sz 0;
      0  0  0 1];
end

function T=translation(xc, yc, zc) 
T = [1 0 0 xc;
     0 1 0 yc;
     0 0 1 zc;
     0 0 0 1];
end

function T=rotationZ(angle)
c = cos(angle); s = sin(angle);
T = [ c -s  0  0;
      s  c  0  0;
      0  0  1  0;
      0  0  0  1];
end


function T=rotationX(angle)
c = cos(angle); s = sin(angle);
T = [ 1  0  0  0;
      0  c -s  0;
      0  s  c  0;
      0  0  0  1];
end

function T=rotationY(angle)
c = cos(angle); s = sin(angle);
T = [ c  0  s  0;
      0  1  0  0;
     -s  0  c  0;
      0  0  0  1];
end

function T=rotation(angle, axis)
if angle==0 | norm(axis)==0, 
  T = eye(4);
else
  phi = atan2(axis(1), axis(2));
  psi = acos(norm(axis(1:2)) / norm(axis));
  T = rotationZ(-phi) * rotationX(-psi) * rotationY(-angle) * rotationX(psi) * rotationZ(phi);
end
end

function T=perspective(d) 
g = 1 / d;
T = [ 1 0 0 0;
      0 1 0 0;
      0 0 0 0;
      0 0 g 0];
end

function plt_coord()
% braizo koordinaciu plokstumas aktyviam paveikslui
xx=axis; 
fill3([xx(1),xx(1),xx(2),xx(2)],[xx(3),xx(4),xx(4),xx(3)],[1 1 1 1]*0,'c','FaceAlpha',0.2,'EdgeColor','c');
fill3([xx(1),xx(1),xx(2),xx(2)],[1 1 1 1]*0,[xx(5),xx(6),xx(6),xx(5)],'c','FaceAlpha',0.2,'EdgeColor','c');
fill3([1 1 1 1]*0,[xx(3),xx(3),xx(4),xx(4)],[xx(5),xx(6),xx(6),xx(5)],'c','FaceAlpha',0.2,'EdgeColor','c');
line([xx(1),xx(2)],[1 1]*0,[1 1]*0,'Color','c','LineStyle','--','Linewidth',1.5);
line([1 1]*0,[xx(3),xx(4)],[1 1]*0,'Color','c','LineStyle','--','Linewidth',1.5);
line([1 1]*0,[1 1]*0,[xx(5),xx(6)],'Color','c','LineStyle','--','Linewidth',1.5);
end