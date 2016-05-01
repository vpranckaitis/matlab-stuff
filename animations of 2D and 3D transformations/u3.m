function u3
figure(1), grid on, hold on, axis equal, axis(7*[-1 1 -1 1 -0.5 1.5]), view([1 -1 1]), xlabel('x'), ylabel('y'), zlabel('z'), caxis([0 2]);
R = 5;
O = 1;

w = [];

dt = 0.01;

for t=0:dt:100
  s = sin(t);
  c = cos(t);
    
  if cos(t + dt) - cos(t) > 0
      colormap('spring');
  else
      colormap('autumn');
  end
  th = ((0:20:360)*pi/180)' + c * O;
  rr = R + (1 - c) * (R / 2);
  r = [0:0.25:rr rr];
  V = th(:,ones(1,size(r,2)));
  U = r(ones(1,size(th,1)), :);
  X = U .* cos(V);
  Y = U.*sin(V);
  Z = sin(U) * s;
  C = 2 - (Z - min(min(Z)));
  
  delete(w);
  w = surf(X, Y, Z, C);
  drawnow;
  pause(dt);
end