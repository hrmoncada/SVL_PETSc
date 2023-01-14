%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all   % close all open such as : figures, fuctions, etc
clc         % clear the command prompt
clear all   % clear all variables

%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Loading Binary Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
NPHI = dlmread('OUTPUT_PHI.mat');
NS   = dlmread('OUTPUT_S.mat');
NUC  = dlmread('OUTPUT_UC.mat');

m = length(NPHI);

PHI = nonzeros(NPHI(3:m));
S   = nonzeros(NS(3:m));
UC  = nonzeros(NUC(3:m));

m = length(PHI);
M = cbrt(length(PHI));
%M =  round(nthroot(length(PHI), 3))

PHI = reshape(PHI,[M,M,M]);
S   = reshape(S,[M,M,M]);
UC  = reshape(UC,[M,M,M]);

[m,n,p] = size(PHI)

% Grid size
Lx = 1; % x-axis unit cell grid size
Ly = 1; % y-axis unit cell grid size
Lz = 1; % y-axis unit cell grid size

% Step size
dx = Lx/(M-1); % x-axis step size
dy = Ly/(M-1); % x-axis step size
dz = Lz/(M-1); % x-axis step size

x = 0 : dx : 1;
y = 0 : dy : 1;
z = 0 : dz : 1;

[X, Y, Z] = meshgrid(x, y, z);

%figure (1)
%  %subplot(1,3,1);
%  [T, p] = isosurface(X, Y, Z, PHI, .5);
%  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
%  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
%  %set(pa,'EdgeColor', 'none')
%  %set(pa,'FaceColor','interp');
%  axis equal; 
%  view(-30, 30)
%  xlabel('x');
%  ylabel('y');
%  zlabel('z');
%  title('PHI')
%  grid on
%  % SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%%saveas (2,"test1.eps")
%%print (2,"test2.eps")
%%print('figure.eps','-deps')
%print -deps -color OUTPUT_PHI.eps
%
%figure (2)
%  %subplot(1,3,2);
%  [T, p] = isosurface(X, Y, Z, S, .5);
%  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
%  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
%  %set(pa,'EdgeColor', 'none')
%  %set(pa,'FaceColor','interp');
%  axis equal; 
%  view(-30, 30)
%  xlabel('x');
%  ylabel('y');
%  zlabel('z');
%  title('S')
%  grid on
%% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%%saveas (2,"test1.eps")
%%print (2,"test2.eps")
%%print('figure.eps','-deps')
%print -deps -color OUTPUT_S.eps
%
%
%figure (3)
%  %subplot(2,3,1);
%  [T, p] = isosurface(X, Y, Z, UC, .5);
%  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
%  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
%  %set(pa,'EdgeColor', 'none')
%  %set(pa,'FaceColor','interp');
%  axis equal; 
%  view(-30, 30)
%  xlabel('x');
%  ylabel('y');
%  zlabel('z');
%  title('UC (-30, 30)')
%  grid on
%% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%%saveas (2,"test1.eps")
%%print (2,"test2.eps")
%%print('figure.eps','-deps')
%print -deps -color OUTPUT_UC.eps

figure (4)
  subplot(2,3,1);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC (-30, 30)')
  grid on
  
  subplot(2,3,2);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(30, -30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC (30, -30)')
  grid on
  
  subplot(2,3,3);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, -30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC (-30, -30)')
  grid on
  
  subplot(2,3,4);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-0, 0)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC (-0, 0)')
  grid on
  
  subplot(2,3,5);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-90, 90)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC(-90, 90)')
  grid on
  
  subplot(2,3,6);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-180, 90)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC(-180, 90)')
  grid on
  
  
% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_OUT.eps
pause()
close all

%%box("off")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELIMINATION_GRATING_AMPLITUD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ('OUTPUT_KAMN')
%
%M1 = length(OUTPUT_KAMN(:,1));      %  matrix  (NM x NN)
%
%for j = 1 : M1    U_KAMN(:,j) =  OUTPUT_KAMN(:,j); end
%
%% Step size
%dx = Lx/(M1-1); % x-axis step size
%dy = Ly/(M1-1); % x-axis step size
%x = 0 : dx : 1;
%y = 0 : dx : 1;
%
%h =  figure (2)
%  imagesc(x, y, U_KAMN)     %  matrix  (NM x NN)
%  axis equal tight
%  xlabel('x')
%  ylabel('y')
%  title('IMPROMENT 1: ELIMINATION GRATING ACCORDING TO THEIR AMPLITUD')
%  grid on
%
%W = 2.5; H = 2.5;
%set(h,'PaperUnits','inches')
%set(h,'PaperOrientation','portrait');
%set(h,'PaperSize',[H,W])
%set(h,'PaperPosition',[0,0,W,H])
%print(h,'-deps','-color','OUTPUT_ELIMINATION_GRATING_AMPLITUD.eps');
%pause() 
%print -deps -color OUTPUT_ELIMINATION_GRATING_AMPLITUD.eps
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IDENTIFIED_COLLINEAR_PLANAR_GRATING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ('OUTPUT_KPQ')
%
%M2 = length(OUTPUT_KPQ(:,1));    %  matrix  (NM x NN)
%
%for j = 1 : M2    U_KPQ(:,j) =  OUTPUT_KPQ(:,j); end
%
%% Step size
%dx = Lx/(M2-1); % x-axis step size
%dy = Ly/(M2-1); % x-axis step size
%x = 0 : dx : 1;
%y = 0 : dx : 1;
%
%h = figure (3)
%  imagesc(x, y, U_KPQ)      %  matrix  (NM x NN)
%  axis equal tight
%  xlabel('x')
%  ylabel('y')
%  title('IMPROMENT 2: IDENTIFIED COLLINEAR PLANAR GRATING')
%  %grid on
%
%W = 2.5; H = 2.5;
%set(h,'PaperUnits','inches')
%set(h,'PaperOrientation','portrait');
%set(h,'PaperSize',[H,W])
%set(h,'PaperPosition',[0,0,W,H])
%print(h,'-deps','-color','OUTPUT_COLLINEAR_PLANAR_GRATING.eps');
%print -deps -color OUTPUT_COLLINEAR_PLANAR_GRATING.eps
%pause() 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IDENTIFIED_COLLINEAR_PLANAR_GRATING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ('OUTPUT_KC')
%
%M3 = length(OUTPUT_KC(:,1));    %  matrix  (NM x NN)
%
%for j = 1 : M3    U_KC(:,j) =  OUTPUT_KC(:,j); end
%
%% Step size
%dx = Lx/(M3-1); % x-axis step size
%dy = Ly/(M3-1); % x-axis step size
%x = 0 : dx : 1;
%y = 0 : dx : 1;
%
%h = figure (4)
%  imagesc(x, y, U_KC)      %  matrix  (NM x NN)
%  xlabel('x')
%  ylabel('y')
%  title('IMPROMENTS 1 & 2')
%  %grid on
%
%W = 2.5; H = 2.5;
%set(h,'PaperUnits','inches')
%set(h,'PaperOrientation','portrait');
%set(h,'PaperSize',[H,W])
%set(h,'PaperPosition',[0,0,W,H])
%print(h,'-deps','-color','OUTPUT_IMPROVEMENTS.eps');
%print -deps -color OUTPUT__IMPROVEMENTS.eps
%pause() 


