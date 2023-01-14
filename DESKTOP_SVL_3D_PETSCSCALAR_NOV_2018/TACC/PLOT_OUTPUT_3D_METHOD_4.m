%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all   % close all open such as : figures, fuctions, etc
clc         % clear the command prompt
clear all   % clear all variables

%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Loading Binary Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the numeric values in the file. Specify a space delimiter, a row offset of 1, and a column offset of 0.
% filename = 'dlmlist.txt';
% M = dlmread(filename,' ',1,0)
%NPHI = dlmread('OUTPUT_PHI.mat',' ',4,0);
%NS   = dlmread('OUTPUT_S.mat',' ',4,0);
%NUC  = dlmread('OUTPUT_UC.mat',' ',4,0);

NPHI = dlmread('OUTPUT_PHI.mat');
NS   = dlmread('OUTPUT_S.mat');
NUC  = dlmread('OUTPUT_UC.mat');

%NPHI(1:20,1)
m = length(NPHI)
j = 1;
for i = 1:m
  if (NPHI(i,1) == 0) 
     printf("%d %f %f %f\n", i,NPHI(i,1),NS(i,1),NUC(i,1))
  else   
     PHI(j,1) = NPHI(i,1);
     S(j,1) = NS(i,1);
     UC(j,1) = NUC(i,1);
     j++;
  end   
end
%NPHI1 
m = length(PHI)  
M = cbrt(length(PHI))
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

figure (1)
  %subplot(1,3,1);
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
  % SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_PHI.eps

figure (2)
  %subplot(1,3,2);
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
% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_S.eps


figure (3)
  %subplot(2,3,1);
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
clear all
