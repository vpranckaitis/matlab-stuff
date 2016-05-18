function u1
close all; clc;

FV = stlread('T1_8.stl');

midX=mean(FV.vertices(FV.vertices(:,3)==0, 1));
midY=mean(FV.vertices(FV.vertices(:,3)==0, 2));

FV.vertices(:,1) = FV.vertices(:,1) - midX;
FV.vertices(:,2) = FV.vertices(:,2) - midY;

EPS = 1e-6;

nLayers=11;
nAngles=10;
divLayers=4;
divAngles=1;
nLayersInt=(nLayers-1) * divLayers + 1;
nAnglesInt=(nAngles-1) * divAngles + 1;

h=[min(FV.vertices(:,3))+EPS 10+EPS 25+EPS max(FV.vertices(:,3))-EPS];
h=[h(1):(h(2)-h(1))/(floor(nLayers/2) - 1):h(2) h(3):(h(4)-h(3))/(ceil(nLayers/2) - 1):h(4)];
a=0:2*pi/(nAngles - 1):2*pi;

hh=zeros(1,nLayersInt);
for i=1:length(h)-1
    from=(i-1)*divLayers+1;
    to=from+divLayers;
    hh(from:to)=h(i):(h(i+1)-h(i))/divLayers:h(i+1);
end
aa=0:2*pi/(nAnglesInt - 1):2*pi;

[X, Y, Z]=buildMesh(FV, h, a);
[XX, YY, ZZ]=buildMesh(FV, hh, aa);

[XL, YL, ZL]=interpolateLagrange(X, Y, Z, h, a, hh, aa);
[XGS, YGS, ZGS]=interpolateGlobalSplines(X, Y, Z, h, a, divLayers, divAngles);

W=ones(size(X)+[2 0]);  
W([1 end],:)=ones(size(W([1 end],:)))*10;
W([6 end],:)=ones(size(W([6 end],:)))*10;
[XN, YN, ZN]=nurbsApproximation(duplicateHorizontalEdge(X), duplicateHorizontalEdge(Y), duplicateHorizontalEdge(Z), W, divLayers, divAngles);

handles = [];

figure; set(gcf,'numbertitle','off','name','Originalus objektas'); 
handles(end+1)=subplot(1, 1, 1); hold on;
patch(FV, 'FaceColor', [0.5 0.5 0.5]); 
axis equal; grid on; view([1 1 1]);

handles(end+1)=mySurf(X, Y, Z, 'Sparse mesh');
handles(end+1)=mySurf(XX, YY, ZZ, 'Dense mesh');
handles(end+1)=mySurf(XL, YL, ZL, 'Lagrange interpolation');
handles(end+1)=mySurf(XGS, YGS, ZGS, 'Global splines interpolation');
handles(end+1)=mySurf(XN, YN, ZN, 'NURBS aproximation');

drawnow;
Link=linkprop(handles, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);

fprintf('Errors\n');
fprintf('  Lagrange: %.5f\n', error(FV, XL, YL, ZL));
fprintf('  Global splines: %.5f\n', error(FV, XGS, YGS, ZGS));
fprintf('  NURBS: %.5f\n', error(FV, XN, YN, ZN));

end

function [X, Y, Z]=buildMesh(FV, h, a)

nh = length(h); na=length(a);

X=zeros(nh, na);
Y=zeros(nh, na);
Z=zeros(nh, na);

for i=1:nh
    hi=h(i);
    
    F = [];
    for k=1:size(FV.faces,1) 
        V = FV.vertices(FV.faces(k, :), :);
        minZ = min(V(:, 3));
        maxZ = max(V(:, 3));
        if minZ <= hi && maxZ >= hi 
            F = [F k];
        end
    end
    
    for j=1:na 
        Z(i, j) = hi;
        aj = a(j);
        
        P1=[0 0 hi]; P2=[100*cos(aj) 100*sin(aj) hi]; 
        
        dist = 0;
        for k=F 
            V = FV.vertices(FV.faces(k, :), :);
            
            [inSegment, x, y, z] = lineIntersectsTriangle(P1, P2, V(1,:), V(2,:), V(3,:));
            if inSegment == true && norm([x y]) > dist
                X(i, j) = x;
                Y(i, j) = y;
                Z(i, j) = z;
                dist = norm([x y]); 
            end
        end
    end
    
end
end

function h=mySurf(X, Y, Z, title)
figure; set(gcf,'numbertitle','off','name',title); h=subplot(1, 1, 1); hold on;
surf(X, Y, Z, 'FaceColor', [rand rand rand]*0.5+0.25); 
axis equal; grid on; hold on; view([1 1 1]);
end

function [XL, YL, ZL]=interpolateLagrange(X, Y, Z, h, a, hh, aa) 
  XL=lagrangeInterpolation(X, h, a, hh, aa);
  YL=lagrangeInterpolation(Y, h, a, hh, aa);
  ZL=lagrangeInterpolation(Z, h, a, hh, aa);
end

function [XGS, YGS, ZGS]=interpolateGlobalSplines(X, Y, Z, h, a, divLayers, divAngles) 
  tl=divLayers+1; ta=divAngles+1;

  XGS=globalSplinesInterpolation(X, h, a, tl, ta);
  YGS=globalSplinesInterpolation(Y, h, a, tl, ta);
  ZGS=globalSplinesInterpolation(Z, h, a, tl, ta);
  
  selectorH=[find(mod([1:size(XGS,1)],divLayers+1)) size(XGS,1)];
  selectorA=[find(mod([1:size(XGS,2)],divAngles+1)) size(XGS,2)];
  XGS=XGS(selectorH, selectorA);
  YGS=YGS(selectorH, selectorA);
  ZGS=ZGS(selectorH, selectorA);
end

function e=error(FV, X, Y, Z)
  [XX, YY, ZZ]=buildMesh(FV, Z(:,1), atan2(Y(1,:), X(1,:)));
  e=squareError(X, XX) + squareError(Y, YY) + squareError(Z, ZZ);
end

function e=squareError(A, B)
  e=sum(sum((A-B).^2));
end

function AA=duplicateHorizontalEdge(A) 
  AA=A([1 1:end end],:);
end