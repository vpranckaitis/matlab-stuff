function comparison
pos=u3
dts=[0.01:0.01:0.1 0.1:0.1:1];
for dt=dts
  pos1=u1(dt); pos2=u2(dt);
  e1=error(pos-pos1); e2=error(pos-pos2);
  fprintf('dt=%f  e1=%e  e2=%e\n', dt, e1, e2);
end
end

function e=error(pos)
e=0;
for i=1:size(pos,1)
    e=e+norm(pos(i,:));
end
end