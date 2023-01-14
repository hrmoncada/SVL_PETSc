%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all   % close all open such as : figures, fuctions, etc
clc         % clear the command prompt
clear all   % clear all variables

%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Loading Binary Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
PHI = dlmread('OUTPUT_PHI.dat');
S   = dlmread('OUTPUT_S.dat');
UC  = dlmread('OUTPUT_UC.dat');

m = length(PHI);
M = nthroot(m, 3);
%M =  round(nthroot(length(PHI), 3));

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
  subplot(1,3,1);
  [T, p] = isosurface(X, Y, Z, PHI, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  view(3)
  xlabel('x','FontSize',18);
  ylabel('y','FontSize',18');
  zlabel('z','FontSize',18,'Rotation',0);
  set(gca,'FontSize',18);
  axis vis3d tight
  camlight left
  colormap('parula');
  lighting gouraud
  grid on
  box on
  ax = gca;
  ax.BoxStyle = 'full';
  title('PHI','FontSize',18)

  subplot(1,3,2);
  [T, p] = isosurface(X, Y, Z, S, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x','FontSize',18);
  ylabel('y','FontSize',18');
  zlabel('z','FontSize',18,'Rotation',0);
  set(gca,'FontSize',18);
  axis vis3d tight
  camlight left
  colormap('parula');
  lighting gouraud
  grid on
  box on
  ax = gca;
  ax.BoxStyle = 'full';
  title('S','FontSize',18)
  grid on

  subplot(1,3,3);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x','FontSize',18);
  ylabel('y','FontSize',18');
  zlabel('z','FontSize',18,'Rotation',0);
  set(gca,'FontSize',18);
  axis vis3d tight
  camlight left
  colormap('parula');
  lighting gouraud
  grid on
  box on
  ax = gca;
  ax.BoxStyle = 'full';
  title('UC (-30, 30)','FontSize',18)
  grid on

figure (2)
  subplot(2,3,1);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  %view(3)
  view(-30, 30)
  xlabel('x','FontSize',18);
  ylabel('y','FontSize',18');
  zlabel('z','FontSize',18,'Rotation',0);
  set(gca,'FontSize',18);
  axis vis3d tight
  camlight left
  colormap('parula');
  lighting gouraud
  grid on
  box on
  ax = gca;
  ax.BoxStyle = 'full';
  title('UC (-30, 30)','FontSize',18)

  subplot(2,3,2);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  %view(3)
  view(30, -30)
  xlabel('x','FontSize',18);
  ylabel('y','FontSize',18');
  zlabel('z','FontSize',18,'Rotation',0);
  set(gca,'FontSize',18);
  axis vis3d tight
  camlight left
  colormap('parula');
  lighting gouraud
  grid on
  box on
  ax = gca;
  ax.BoxStyle = 'full';
  title('UC (30, -30)','FontSize',18)
  
  subplot(2,3,3);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, -30)
  %view(3)
  xlabel('x','FontSize',18);
  ylabel('y','FontSize',18');
  zlabel('z','FontSize',18,'Rotation',0);
  set(gca,'FontSize',18);
  axis vis3d tight
  camlight left
  colormap('parula');
  lighting gouraud
  grid on
  box on
  ax = gca;
  ax.BoxStyle = 'full';
  title('UC (-30, -30)','FontSize',18)
  
  subplot(2,3,4);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  %view(3)
  view(-0, 0)
  xlabel('x','FontSize',18);
  ylabel('y','FontSize',18','Rotation',0);
  zlabel('z','FontSize',18,'Rotation',0);
  set(gca,'FontSize',18);
  axis vis3d tight
  camlight left
  colormap('parula');
  lighting gouraud
  grid on
  box on
  ax = gca;
  ax.BoxStyle = 'full';
  title('UC (-0, 0)','FontSize',18)
  
  subplot(2,3,5);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  %view(3)
  view(-90, 90)
  xlabel('x','FontSize',18);
  ylabel('y','FontSize',18');
  zlabel('z','FontSize',18,'Rotation',0);
  set(gca,'FontSize',18);
  axis vis3d tight
  camlight left
  colormap('parula');
  lighting gouraud
  grid on
  box on
  ax = gca;
  ax.BoxStyle = 'full';
  title('UC(-90, 90)','FontSize',18)
  
  subplot(2,3,6);
  [T, p] = isosurface(X, Y, Z, UC, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-180, 90)
  %view(3)
  xlabel('x','FontSize',18);
  ylabel('y','FontSize',18');
  zlabel('z','FontSize',18,'Rotation',0);
  set(gca,'FontSize',18);
  axis vis3d tight
  camlight left
  colormap('parula');
  lighting gouraud
  grid on
  box on
  ax = gca;
  ax.BoxStyle = 'full';
  title('UC(-180, 90)','FontSize',18)

   pause
   close all

