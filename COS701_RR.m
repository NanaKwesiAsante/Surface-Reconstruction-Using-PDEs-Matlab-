% Vivian Montiforte
% Kwesi Acheampong

% COS 701
% Fall 2019
% Research Project

clear variables;

%%% Loading/Storing Data
ptCloud = pcread('vase.ply');
xyz = ptCloud.Location;
normals = pcnormals(ptCloud);
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
N = length(x);

%%% Cleaning Data (x and z)
zz = 0;
xx = 0;
ind = 0;
for i = 1:N
    zi = z(i);
    xi = x(i);
    if zi < 4 && zi > -4 && xi < 4 && xi > -4
        zz = [zz; zi];
        xx = [xx; xi];
        ind = [ind; i];
    end
    
end
ind = ind(2:end,:);
zz = zz(2:end,:);
xx = xx(2:end,:);

%%% Matching Indices with Cleaned Data (x and z)
yy = zeros(length(ind),1);
norms1 = zeros(length(ind),3);
for j = 1:length(ind)
    in = ind(j);
    yy(j) = y(in);
    norms1(j) = normals(in);
end
    
%%% Cleaning Data (y)
indy = 0;
yyy = 0;
for i = 1:length(yy)
    zi = zz(i);
    if zi >= -4 && zi <= 0
        yi = yy(i);
        if yi > -1.6
            yyy = [yyy; yi];
            indy = [indy; i];
        end
    end
    if zi > 0 && zi <= 4
        yi = yy(i);
        if yi > -1.6
            yyy = [yyy; yi];
            indy = [indy; i];
        end
    end
end
yyy = yyy(2:end,:);
indy = indy(2:end,:);

%%% Matching Indices with Cleaned Data (y)
xxx = zeros(length(indy),1);
zzz = zeros(length(indy),1);
norms = zeros(length(indy),3);
for j = 1:length(indy)
    in = indy(j);
    xxx(j) = xx(in);
    zzz(j) = zz(in);
    norms(j) = norms1(in);
end

%%% Choosing every h points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 20;     % test different values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xxx = xxx(1:h:end,:);
yyy = yyy(1:h:end,:);
zzz = zzz(1:h:end,:);
norms = norms(1:h:end,:);

%%% Checking normal vector direction
sensorCenter = [0, 5, -1];
for j = 1:length(xxx)
    p1 = sensorCenter - [xxx(j),yyy(j),zzz(j)];
    p2 = [norms(j,1),norms(j,2),norms(j,3)];
    ang = atan2(norm(cross(p1,p2)),p1*p2');
    if ang > pi/2 || ang < -pi/2
        norms(j,1) = -norms(j,1);
        norms(j,2) = -norms(j,2);
        norms(j,3) = -norms(j,3);
    end
end
 
%%% Creating intpts using normal vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp = .25;     % test different values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bdpts = [xxx, yyy, zzz];
intpts = bdpts(1:20:end,:) + sp*norms(1:20:end,:);

%%% Plotting Boundary and Interior Points
figure('Name','Boundary Points and Interior Points',...
    'ToolBar','none','NumberTitle','off')
plot3(intpts(:,1),intpts(:,2),intpts(:,3),'r.',...
    'MarkerSize',10);
hold on;
plot3(bdpts(:,1),bdpts(:,2),bdpts(:,3),'b.',...
    'MarkerSize',3);
axis equal
view(1,70)
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
xlabel('x','FontSize',18,'FontWeight','bold');
ylabel('y','FontSize',18,'FontWeight','bold');
zlabel('z','FontSize',18,'FontWeight','bold');

%%% xyz intpts and xyz bdpts; centers
ni = length(intpts)
nb = length(bdpts)
n = ni + nb;
xbd = bdpts(:,1); ybd = bdpts(:,2); zbd = bdpts(:,3);
xint = intpts(:,1); yint = intpts(:,2); zint = intpts(:,3);
ctrs(1:ni,:) = intpts; ctrs(ni+1:n,:) = bdpts;

%%% Euclidean Distance between points and centers
DMint = DistanceMatrix([xint yint zint],ctrs);
DMbd = DistanceMatrix([xbd, ybd, zbd],ctrs);
b = ones(n,1);  % Right hand side

%%% Second-Order Matrix
%MULTIQUADRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lam = 100;  % test different values
%c = 1000;   % test different values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rbf = @(r,c) sqrt(1 +r.^2*c^2);
% A(1:ni,:) = (3*c^2+2*c^4*(DMint).^2)./(rbf(DMint,c).^3)-lam*rbf(DMint,c);
% %A(1:ni,:) =(c^4*(DMint.^2)+2*c^2)./(rbf(DMint,c).^3)-lam*rbf(DMint,c);
% A(ni+1:n,:) = rbf(DMbd,c);
% coef = A\b;

% %%% Fourth-Order Matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lam = 500;  % test different values
% c = 200;    % test different values
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rbf = @(r,c) sqrt(1 +r.^2*c^2);
% rbf_l2 = @(r,c) -15./(1+r.^2*c^2).^(7/2);
% A(1:ni,:) = rbf_l2(DMint,c) - lam*rbf(DMint,c);
% A(ni+1:n,:) = rbf(DMbd,c);
% coef = A\b;

% %%% Inverse MultiQuadric RBF
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %lam = 1000;  % test different values
   %c = 0.3;    % test different values
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rbf = @(r,c) 1./sqrt(r.^2+c^2);
%  %A(1:ni,:) = rbf(DMint,c).*(1 - c./(c+sqrt(DMint.^2+c^2)))-lam*rbf(DMint,c);
% A(1:ni,:) =((DMint.^2)-2*c^2).*(rbf(DMint,c).^3)-lam*rbf(DMint,c);
% A(ni+1:n,:) = rbf(DMbd,c);
% coef = A\b;

%Gaussian Rbf
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  lam = 100000;  % test different values
%  c = 5;    % test different values
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rbf = @(r,c) exp(-(c.^2*r.^2));
% A(1:ni,:) = (4*c^2*rbf(DMint,c).*((c^2*DMint.^2)-1))-lam*rbf(DMint,c);
% A(ni+1:n,:) = rbf(DMbd,c);
% coef = A\b;



%%% Creating Evaluation Points
xmin = min(xbd); xmax = max(xbd);
ymin = min(ybd); ymax = max(ybd);
zmin = min(zbd); zmax = max(zbd);
neval = 40;
xe = linspace(xmin-0.01,xmax+0.01,neval);
ye = linspace(ymin-0.01,ymax+0.01,neval);
ze = linspace(zmin-0.01,zmax+0.01,neval);
[XE, YE, ZE] = meshgrid(xe,ye,ze);
xev = XE(:); yev = YE(:); zev = ZE(:);
epts = [ xev yev zev];

%%% Reconstructing Surface
nn = length(epts)/20;
for i = 1:20
    DMeval = DistanceMatrix(epts((i-1)*nn+1:i*nn,:,:),ctrs); % create r for current interval
    u_hat((i-1)*nn+1:i*nn,1) = rbf(DMeval,c)*coef;  % store approx for current interval
end
u_hat = reshape(u_hat,neval,neval,neval);   

%%% Plotting Surface Reconstruction
black = [ 0 0 0 ];
figure('Name','Surface Reconstruction',...
    'ToolBar','none','NumberTitle','off','Color',black)
isosurface(XE,YE,ZE,u_hat,1)
axis equal
light('Position',[0 0 1],'Style','infinite');
light('Position',[0 0 -1],'Style','infinite');
colormap autumn
material shiny
axis off
view(1,70)
            
        