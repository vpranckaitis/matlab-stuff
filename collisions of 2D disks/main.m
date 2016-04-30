function main
figure; axis equal; axis([0 20 -2 18]); hold on;

k=1000;
c0=2;
cf=1;
f=0.5;

elems=[
    Boundary([0 18], [0 8])
    Boundary([0 8], [4 4])
    Boundary([4 4], [6 6])
    Boundary([6 6], [10 2])
    Boundary([10 2], [12 4])
    Boundary([12 4], [16 0])
    Boundary([16 0], [20 0])
    Boundary([20 0], [20 18])
    Disk(0.5, 0.1, [0.5 20], [0 0])
];

h=[];

dt=0.005;

for l=1:20000
    if mod(l, 500)==0
        rr=rand();
        r=0.5 + rr*1;
        m=0.1 + rr*1;
        d=Disk(r, m, [r 20], [0 0]);
        elems=[elems; d];
    end
    
    for i=1:length(elems)
        for j=(i+1):length(elems)
            [u N]=checkCollision(elems(i), elems(j));
            Tau=[N(2) -N(1)];
            if u > 0
                er=elems(i); es=elems(j);
                Vr=er.V - Tau*er.Vr;
                Vs=es.V + Tau*es.Vr;
                
                c=damping(er, es, k, c0);
                
                Trs=k*u + dot(Vr-Vs, N)*c; Trs=max([Trs 0]);
                er.applyForce(-Trs*N); es.applyForce(Trs*N);
                
                Frs=cf*dot(Tau, Vr-Vs);
                if abs(Frs) > f*(Trs), Frs=f*(Trs)*sign(Frs); end
                er.applyForce(-Frs*Tau); es.applyForce(Frs*Tau);
                
                if strcmp(class(er), 'Disk')
                    Sr=N*er.r; Mr=cross([Sr 0], [(-Frs*Tau) 0]);
                    er.applyForceMoment(Mr(3));
                end
                
                if strcmp(class(es), 'Disk')
                    Ss=-N*es.r; Ms=cross([Ss 0], [(Frs*Tau) 0]);
                    es.applyForceMoment(Ms(3));
                end
            end
        end
    end
    
    for i=1:length(elems)
        elems(i).processForces(dt);
    end
    
    if mod(l, 10)==0 
        delete(h);
        for i=1:length(elems)
            h(i) = elems(i).plot();
        end
        pause(0.01);
    end
end
end

function [d N]=checkCollision(er, es)
if strcmp(class(er), 'Disk') && strcmp(class(es), 'Disk')
    c=sqrt(sum((er.P - es.P).^2));
    d=er.r + es.r - c;
    d=max([0 d]);
    N=es.P - er.P;
    N=N/norm(N);
elseif strcmp(class(er), 'Boundary') && strcmp(class(es), 'Disk')
    V1=er.B - er.A;
    V2=es.P - er.A;
    proj=dot(V1, V2)/norm(V1);
    C=V1/norm(V1)*proj + er.A;
    if proj < 0, C=er.A;
    elseif proj > norm(V1), C=er.B;
    end
    d=norm(es.P - C);
    d=es.r - d;
    d=max([0 d]);
    N=es.P-C;
    N=N/norm(N);
elseif strcmp(class(er), 'Disk') && strcmp(class(es), 'Boundary')
    [d N]=checkCollision(es, er);
    N=N * (-1);
else
    d=0; N=[0 0];
end
end

function c=damping(er, es, k, c0)
mm=1; ms=0;
if ~er.isStatic()
   mm = mm*er.m;
   ms = ms+er.m;
end
if ~es.isStatic()
   mm = mm*es.m;
   ms = ms+es.m;
end
ms=max([1 ms]);
c=c0*sqrt(mm*k/ms);
end
