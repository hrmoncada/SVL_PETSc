%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all   % close all open such as : figures, fuctions, etc
clc         % clear the command prompt
clear all   % clear all variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LOADING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NPHI0 = dlmread('OUTPUT_PHI_0.mat');
NS0   = dlmread('OUTPUT_S_0.mat');
NUC0  = dlmread('OUTPUT_UC_0.mat');

NPHI = dlmread('OUTPUT_PHI.mat');
NS   = dlmread('OUTPUT_S.mat');
NUC  = dlmread('OUTPUT_UC.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  RESHAPE DATA INTO 3D ARRAYS (ZERO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = length(NPHI0);

PHI0 = nonzeros(NPHI0(3:m));
S0   = nonzeros(NS0(3:m));
UC0  = nonzeros(NUC0(3:m));

m = length(PHI0);
M = cbrt(length(PHI0)); % cubic root
%M0 =  round(nthroot(length(PHI0), 3))

PHI0 = reshape(PHI0,[M,M,M]);
S0   = reshape(S0,[M,M,M]);
UC0  = reshape(UC0,[M,M,M]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  RESHAPE DATA INTO 3D ARRAYS (LAST)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = length(NPHI);

PHI = nonzeros(NPHI(3:m));
S   = nonzeros(NS(3:m));
UC  = nonzeros(NUC(3:m));

m = length(PHI);
M = cbrt(length(PHI)); % cubic root
%M =  round(nthroot(length(PHI), 3))

PHI = reshape(PHI,[M,M,M]);
S   = reshape(S,[M,M,M]);
UC  = reshape(UC,[M,M,M]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SUBSTRACT 3D ARRAYS  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PHI_RES = PHI - PHI0;
S_RES   = S - S0;
UC_RES  = UC - UC0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PLOT DATA 3D ARRAYS SIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n,p] = size(PHI0);
[m,n,p] = size(PHI);
[m,n,p] = size(PHI_RES);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PLOT DATA 3D ARRAYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure (1)
  subplot(2,3,1);
  [T, p] = isosurface(X, Y, Z, PHI0, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('PHI_0')
  grid on

  subplot(2,3,2);
  [T, p] = isosurface(X, Y, Z, S0, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');  
  ylabel('y');
  zlabel('z');
  title('S_0')
  grid on

  subplot(2,3,3);
  [T, p] = isosurface(X, Y, Z, UC0, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC_0 (-30, 30)')
  grid on
 
  subplot(2,3,4);
  [T, p] = isosurface(X, Y, Z, PHI, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('PHI')
  grid on

  subplot(2,3,5);
  [T, p] = isosurface(X, Y, Z, S, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('S')
  grid on

  subplot(2,3,6);
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
  
% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_UC_0.eps

figure (2)
  subplot(2,3,1);
  [T, p] = isosurface(X, Y, Z, UC0, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC_0 (-30, 30)')
  grid on
  
  subplot(2,3,2);
  [T, p] = isosurface(X, Y, Z, UC0, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(30, -30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC_0 (30, -30)')
  grid on

  subplot(2,3,3);
  [T, p] = isosurface(X, Y, Z, UC0, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, -30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC_0 (-30, -30)')
  grid on
  
  subplot(2,3,4);
  [T, p] = isosurface(X, Y, Z, UC0, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-0, 0)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC_0 (-0, 0)')
  grid on
  
  subplot(2,3,5);
  [T, p] = isosurface(X, Y, Z, UC0, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',figure (5)T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-90, 90)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC_0 (-90, 90)')
  grid on
  
  subplot(2,3,6);
  [T, p] = isosurface(X, Y, Z, UC0, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-180, 90)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC_0 (-180, 90)')
  grid on
  
  
% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_OUT_0.eps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Loading Binary Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure (3)
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
  title('UC (-90, 90)')
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
  title('UC (-180, 90)')
  grid on
  
% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_OUT.eps

figure (4)
  subplot(1,3,1);
  [T, p] = isosurface(X, Y, Z, PHI_RES, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('PHI - PHI_0 (-30, 30)')
  grid on
  
  subplot(1,3,2);
  [T, p] = isosurface(X, Y, Z, S_RES, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('S - S_0 (-30, 30)')
  grid on
  
  subplot(1,3,3);
  [T, p] = isosurface(X, Y, Z, UC_RES, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('UC - UC_0 (-30, 30)')
  grid on
  
% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_OUT.eps
pause()
close all

