% Vivian Montiforte
% Kwesi Acheampong

% COS 701
% Fall 2019
% Research Project

clear variables;

load fig14b2.txt
bdpts = fig14b2;
intpts = [ 36 14 -10 ];
ni = 1

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
% ni = length(intpts);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lam = 100;  % test different values
%c = 1000;   % test different values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rbf = @(r,c) sqrt(1 +r.^2*c^2);
% A(1:ni,:) = (3*c^2+2*c^4*(DMint).^2)./(rbf(DMint,c).^3)-lam*rbf(DMint,c);
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
   %lam = 100;  % test different values
   %c = 3;    % test different values
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
  lam = 10000000;  % test different values
  c = 0.4;    % test different values
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rbf = @(r,c) exp(-(c.^2*r.^2));
A(1:ni,:) = (4*c^2*rbf(DMint,c).*((c^2*DMint.^2)-1))-lam*rbf(DMint,c);
A(ni+1:n,:) = rbf(DMbd,c);
coef = A\b;





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