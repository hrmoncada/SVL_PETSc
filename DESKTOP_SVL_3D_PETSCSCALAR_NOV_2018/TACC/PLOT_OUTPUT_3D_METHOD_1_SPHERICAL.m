%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all   % close all open such as : figures, fuctions, etc
clc         % clear the command prompt
clear all   % clear all variables
clf         % clear functions

% Grid size
Lx = 1; % x-axis unit cell grid size
Ly = 1; % y-axis unit cell grid size
Lz = 1; % y-axis unit cell grid size

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Loading Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
load ('OUTPUT_UNIT_CELL_REAL')
load ('OUTPUT_UNIT_CELL_IMAG')
load ('OUTPUT_FFTW_REAL')
load ('OUTPUT_FFTW_IMAG')
load ('OUTPUT_INV_FFTW_REAL')
load ('OUTPUT_INV_FFTW_IMAG')
load ('OUTPUT_SWAP_REAL')
load ('OUTPUT_SWAP_IMAG')
load ('OUTPUT_TRUNC_FFTW_REAL')
load ('OUTPUT_TRUNC_FFTW_IMAG')
load ('OUTPUT_KX')
load ('OUTPUT_KY')
load ('OUTPUT_KZ')
load ('OUTPUT_X')
load ('OUTPUT_Y')
load ('OUTPUT_Z')
load ('OUTPUT_THETA')
load ('OUTPUT_RSQ')
load ('OUTPUT_VARPHI')
load ('OUTPUT_PER')
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Length Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
M1_rows = length(OUTPUT_UNIT_CELL_REAL(:,1));       % Unit Cell matrix    	(Nx x Ny x Nz)
M1_cols = length(OUTPUT_UNIT_CELL_REAL(1,:));       % Unit Cell matrix    	(Nx x Ny x Nz)

M2_rows = length(OUTPUT_FFTW_REAL(:,1));       	    % FFTW REAL matrix    	(Nx x Ny x Nz)
M2_cols = length(OUTPUT_FFTW_REAL(1,:));            % FFTW REAL matrix    	(Nx x Ny x Nz)

M3_rows = length(OUTPUT_INV_FFTW_REAL(:,1));        % FFTW REAL matrix    	(Nx x Ny x Nz)
M3_cols = length(OUTPUT_INV_FFTW_REAL(1,:));        % FFTW REAL matrix    	(Nx x Ny x Nz)

M4_rows = length(OUTPUT_SWAP_REAL(:,1));            % FFTW REAL matrix    	(Nx x Ny x Nz)
M4_cols = length(OUTPUT_SWAP_REAL(1,:));            % FFTW REAL matrix    	(Nx x Ny x Nz)

M5_rows = length(OUTPUT_TRUNC_FFTW_REAL(:,1));      % TRUC FFTW REAL matrix     (NM x NN x NP)
M5_cols = length(OUTPUT_TRUNC_FFTW_REAL(1,:));      % TRUC FFTW REAL matrix     (NM x NN x NP)

M6_rows = length(OUTPUT_KX(:,1)); 		  %  KX matrix      	(NM x NN x NP)
M6_cols = length(OUTPUT_KX(1,:)); 		  %  KX matrix      	(NM x NN x NP)

M7_rows = length(OUTPUT_KY(:,1)); 		  %  KY matrix      	(NM x NN x NP)
M7_cols = length(OUTPUT_KY(1,:)); 		  %  KY matrix      	(NM x NN x NP)

M8_rows = length(OUTPUT_KZ(:,1)); 		  %  KZ matrix      	(NM x NN x NP)
M8_cols = length(OUTPUT_KZ(1,:)); 		  %  KZ matrix      	(NM x NN x NP)

M9_rows = length(OUTPUT_X(:,1)); 		  %  X matrix                 (New_Nx x New_Ny x New_Nz)
M9_cols = length(OUTPUT_X(1,:)); 		  %  X matrix                 (New_Nx x New_Ny x New_Nz)

M10_rows = length(OUTPUT_Y(:,1)); 		  %  Y matrix                 (New_Nx x New_Ny x New_Nz)
M10_cols = length(OUTPUT_Y(1,:)); 		  %  Y matrix                 (New_Nx x New_Ny x New_Nz)

M11_rows = length(OUTPUT_Z(:,1)); 		  %  Z matrix                 (New_Nx x New_Ny x New_Nz)
M11_cols = length(OUTPUT_Z(1,:)); 		  %  Z matrix                 (New_Nx x New_Ny x New_Nz)

M12_rows = length(OUTPUT_RSQ(:,1));   	          %  RSQ matrix               (New_Nx x New_Ny x New_Nz)
M12_cols = length(OUTPUT_RSQ(1,:));   	          %  RSQ matrix    	      (New_Nx x New_Ny x New_Nz)

M13_rows = length(OUTPUT_THETA(:,1)); 	          %  THETA matrix             (New_Nx x New_Ny x New_Nz)
M13_cols = length(OUTPUT_THETA(1,:)); 	          %  THETA matrix    	      (New_Nx x New_Ny x New_Nz)

M14_rows = length(OUTPUT_VARPHI(:,1)); 		  %  X matrix                 (New_Nx x New_Ny x New_Nz)
M14_cols = length(OUTPUT_VARPHI(1,:)); 		  %  X matrix                 (New_Nx x New_Ny x New_Nz)

M15_rows = length(OUTPUT_PER(:,1));   	          %  PER matrix    	      (New_Nx x New_Ny x New_Nz)
M15_cols = length(OUTPUT_PER(1,:));   	          %  PER matrix    	      (New_Nx x New_Ny x New_Nz)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Array Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1
for j = 1 : M1_cols   U1(:,:,j) = OUTPUT_UNIT_CELL_REAL((j-1)*M1_cols+1 : j*M1_cols,:); end
for j = 1 : M1_cols   U2(:,:,j) = OUTPUT_UNIT_CELL_IMAG((j-1)*M1_cols+1 : j*M1_cols,:); end
% figure 2
for j = 1 : M2_cols   U3(:,:,j) = OUTPUT_FFTW_REAL((j-1)*M2_cols+1 : j*M2_cols,:); end
for j = 1 : M2_cols   U4(:,:,j) = OUTPUT_FFTW_IMAG((j-1)*M2_cols+1 : j*M2_cols,:); end
for j = 1 : M3_cols   U5(:,:,j) = OUTPUT_INV_FFTW_REAL((j-1)*M3_cols+1 : j*M3_cols,:); end
for j = 1 : M3_cols   U6(:,:,j) = OUTPUT_INV_FFTW_IMAG((j-1)*M3_cols+1 : j*M3_cols,:); end
for j = 1 : M4_cols   U7(:,:,j) = OUTPUT_SWAP_REAL((j-1)*M4_cols+1 : j*M4_cols,:); end
for j = 1 : M4_cols   U8(:,:,j) = OUTPUT_SWAP_IMAG((j-1)*M4_cols+1 : j*M4_cols,:); end
% figure 3
for j = 1 : M5_cols   U9(:,:,j)  = OUTPUT_TRUNC_FFTW_REAL((j-1)*M5_cols+1 : j*M5_cols,:); end
for j = 1 : M5_cols   U10(:,:,j) = OUTPUT_TRUNC_FFTW_IMAG((j-1)*M5_cols+1 : j*M5_cols,:); end
for j = 1 : M6_cols   U11(:,:,j) = OUTPUT_KX((j-1)*M6_cols+1 : j*M6_cols,:);  end
for j = 1 : M7_cols   U12(:,:,j) = OUTPUT_KY((j-1)*M7_cols+1 : j*M7_cols,:);  end
for j = 1 : M8_cols   U13(:,:,j) = OUTPUT_KZ((j-1)*M8_cols+1 : j*M8_cols,:);  end
% figure 4
for j = 1 : M9_cols   U14(:,:,j) = OUTPUT_X((j-1)*M9_cols+1 : j*M9_cols,:);    end
for j = 1 : M10_cols  U15(:,:,j) = OUTPUT_Y((j-1)*M10_cols+1 : j*M10_cols,:);  end
for j = 1 : M11_cols  U16(:,:,j) = OUTPUT_Z((j-1)*M11_cols+1 : j*M11_cols,:);  end
% figure 5
for j = 1 : M12_cols  U17(:,:,j) =  OUTPUT_RSQ((j-1)*M12_cols+1 : j*M12_cols,:);    end
for j = 1 : M13_cols  U18(:,:,j) =  OUTPUT_THETA((j-1)*M13_cols+1 : j*M13_cols,:);  end
for j = 1 : M14_cols  U19(:,:,j) =  OUTPUT_VARPHI((j-1)*M14_cols+1 : j*M14_cols,:);  end
for j = 1 : M15_cols  U20(:,:,j) =  OUTPUT_PER((j-1)*M15_cols+1 : j*M15_cols,:);    end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subplot Array Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
% Step size
dx = Lx/(M1_cols-1); % x-axis step size
dy = Ly/(M1_cols-1); % y-axis step size
dz = Lz/(M1_cols-1); % z-axis step size

x = 0 : dx : 1;
y = 0 : dy : 1;
z = 0 : dz : 1;

[X, Y, Z] = meshgrid(x, y, z);

  subplot(1,2,1);
  [T, p] = isosurface(X, Y, Z, U1, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal tight; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('3D UNIT CELL REAL')
  grid on

  subplot(1,2,2);
  [T, p] = isosurface(X, Y, Z, U2, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  %axis equal tight; 
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('3D UNIT CELL IMAG')
  grid on

% SAVE PLOTS: saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_3D_1.eps

figure(2)
% Step size
dx = Lx/(M3_cols-1); % x-axis step size
dy = Ly/(M3_cols-1); % y-axis step size
dz = Lz/(M3_cols-1); % z-axis step size

x = 0 : dx : 1;
y = 0 : dy : 1;
z = 0 : dz : 1;

[X, Y, Z] = meshgrid(x, y, z);
  subplot(3,2,1);
  [T, p] = isosurface(X, Y, Z, U3, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal tight; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('REAL FFTW')
  grid on

  subplot(3,2,2);
  [T, p] = isosurface(X, Y, Z, U4, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal tight; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('IMAG FFTW')
  grid on

  subplot(3,2,3);
  [T, p] = isosurface(X, Y, Z, U5, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal tight; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('REAL INV FFTW')
  grid on

  subplot(3,2,4);
  [T, p] = isosurface(X, Y, Z, U6, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal tight; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('IMAG INV FFTW')
  grid on

  subplot(3,2,5);
  [T, p] = isosurface(X, Y, Z, U7, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal tight; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('REAL SWAP')
  grid on

  subplot(3,2,6);
  [T, p] = isosurface(X, Y, Z, U8, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal tight; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('IMAG SWAP')
  grid on

% SAVE PLOTS: saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_3D_2.eps

figure(3) 
% Step size
dx = Lx/(M5_cols-1); % x-axis step size
dy = Ly/(M5_cols-1); % y-axis step size
dz = Lz/(M5_cols-1); % z-axis step size

x = 0 : dx : 1;
y = 0 : dy : 1;
z = 0 : dz : 1;

[X, Y, Z] = meshgrid(x, y, z);

  subplot(3,2,1);
  [T, p] = isosurface(X, Y, Z, U9, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal tight; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('TRUNC REAL FFTW')
  grid on

  subplot(3,2,2);
  [T, p] = isosurface(X, Y, Z, U10, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('TRUNC IMAG FFTW')
  grid on


  subplot(3,2,3);
  [T, p] = isosurface(X, Y, Z, U11, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal tight; 
  axis([0 1 0 1 0 1])
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('KX')
  grid on

  subplot(3,2,4);
  [T, p] = isosurface(X, Y, Z, U12, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal tight; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('KY')
  grid on

  subplot(3,2,5);
  [T, p] = isosurface(X, Y, Z, U13, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  daspect([1,1,1])
  axis equal tight; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('KZ')
  grid on

% SAVE PLOTS: saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_3D_3.eps


figure(4)
% Step size
% Step size
dx = Lx/(M9_cols-1); % x-axis step size
dy = Ly/(M9_cols-1); % y-axis step size
dz = Lz/(M9_cols-1); % z-axis step size

x = 0 : dx : 1;
y = 0 : dy : 1;
z = 0 : dz : 1;

[X, Y, Z] = meshgrid(x, y, z);

  subplot(1,3,1);
  [T, p] = isosurface(X, Y, Z, U14, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('X')
  grid on


  subplot(1,3,2);
  [T, p] = isosurface(X, Y, Z, U15, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('Y')
  grid on

  subplot(1,3,3);
  [T, p] = isosurface(X, Y, Z, U16, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('Z')
  grid on

% SAVE PLOTS: saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_3D_4.eps

figure(5)
% Step size
% Step size
dx = Lx/(M12_cols-1); % x-axis step size
dy = Ly/(M12_cols-1); % y-axis step size
dz = Lz/(M12_cols-1); % z-axis step size

x = 0 : dx : 1;
y = 0 : dy : 1;
z = 0 : dz : 1;

[X, Y, Z] = meshgrid(x, y, z);

  subplot(2,2,1);
  [T, p] = isosurface(X, Y, Z, U17, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('RSQ')
  grid on
  
  subplot(2,2,2);
  [T, p] = isosurface(X, Y, Z, U18, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('THETA')
  grid on

  subplot(2,2,3);
  [T, p] = isosurface(X, Y, Z, U19, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('VARPHI')
  grid on
  
  subplot(2,2,4);
  [T, p] = isosurface(X, Y, Z, U20, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  axis equal; 
  view(-30, 30)
  xlabel('x');
  ylabel('y');
  zlabel('z');
  title('PER')
  grid on

% SAVE PLOTS: saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_3D_5.eps


%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Loading Binary Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
NPHI = dlmread('OUTPUT_PHI.mat');
NS   = dlmread('OUTPUT_S.mat');
NUC  = dlmread('OUTPUT_UC.mat');

m = length(NPHI);

PHI = nonzeros(NPHI(4:m)); % Return a vector of the nonzero values of the sparse matrix s.
S   = nonzeros(NS(4:m));
UC  = nonzeros(NUC(4:m));

m = length(PHI);
M = cbrt(length(PHI));
%M =  round(nthroot(length(PHI), 3))

PHI = reshape(PHI,[M,M,M]);
S   = reshape(S,[M,M,M]);
UC  = reshape(UC,[M,M,M]);

[m,n,p] = size(PHI);


figure(6)
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

  subplot(1,3,1);
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

  subplot(1,3,2);
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
  
  subplot(1,3,3);
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
  title('UC')
  grid on
  

% SAVE PLOTS: saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_3D_6.eps

pause() 
closer all
