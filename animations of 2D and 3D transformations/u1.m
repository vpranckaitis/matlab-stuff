%
omega = 2*pi;
rMax = 20;
rMin = 3;
%

circle = 0:0.1:2*pi;
figure(1); grid on, hold on, axis equal,axis([-(rMax+1) (rMax+1) -(rMax+1) (rMax+1)]);

wing = [0 0 0.3;
        0 1 0.3;
        1 1 1];     
Tr = [ 0 -1  0
       1  0  0
       0  0  1];
       
wings = [wing Tr*wing Tr^2*wing Tr^3*wing];

plot(rMax*sin([circle 0]), rMax*cos([circle 0]), '--');
w = fill(0, 0, 'r');

dt = 1 / 60; rotation = 0;

for t=0:dt:100 
  r = (rMax + rMin) / 2 - cos(t) * (rMax - rMin) / 2;
  
  omg = (r - rMax) / r * omega;
  
  rotation = rotation + omg * dt;
  
  T1 = [r 0 0; 
        0 r 0; 
        0 0 1];
        
  T2 = [cos(rotation) -sin(rotation) 0;
        sin(rotation) cos(rotation)  0;
        0             0              1];
        
  xc = (rMax - r)*cos(omega * t);
  yc = (rMax - r)*sin(omega * t);
  T3 = [1 0 xc;
        0 1 yc;
        0 0 1];
  
  wingT =  T3 * T2 * T1 * wings;
  wingT = wingT ./ wingT([3 3 3], :);
  delete(w);
  w = fill(reshape(wingT(1, :), 4, 3)', reshape(wingT(2, :), 4, 3)', [1 5 7 9]);
  
  drawnow; pause(dt);
end