%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all   % close all open such as : figures, fuctions, etc
clc         % clear the command prompt
clear all   % clear all variables
clf         % clear functions

% Grid size
Sx = 1; % x-axis unit cell grid size
Sy = 1; % y-axis unit cell grid size
Sz = 1; % y-axis unit cell grid size

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Loading Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
load ('OUTPUT_UNIT_CELL_REAL')
load ('OUTPUT_KX')
load ('OUTPUT_KY')
load ('OUTPUT_KZ')
load ('OUTPUT_X')
load ('OUTPUT_Y')
load ('OUTPUT_Z')
load ('OUTPUT_THETA')
load ('OUTPUT_RSQ')
load ('OUTPUT_PER')
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Length Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
M1_rows = length(OUTPUT_UNIT_CELL_REAL(:,1));       % Unit Cell matrix    	(Nx x Ny x Nz)
M1_cols = length(OUTPUT_UNIT_CELL_REAL(1,:));       % Unit Cell matrix    	(Nx x Ny x Nz)

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

M15_rows = length(OUTPUT_PER(:,1));   	          %  PER matrix    	      (New_Nx x New_Ny x New_Nz)
M15_cols = length(OUTPUT_PER(1,:));   	          %  PER matrix    	      (New_Nx x New_Ny x New_Nz)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Array Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1
for j = 1 : M1_cols   U1(:,:,j) = OUTPUT_UNIT_CELL_REAL((j-1)*M1_cols+1 : j*M1_cols,:); end
% figure 2
for j = 1 : M6_cols   U11(:,:,j) = OUTPUT_KX((j-1)*M6_cols+1 : j*M6_cols,:);  end
for j = 1 : M7_cols   U12(:,:,j) = OUTPUT_KY((j-1)*M7_cols+1 : j*M7_cols,:);  end
for j = 1 : M8_cols   U13(:,:,j) = OUTPUT_KZ((j-1)*M8_cols+1 : j*M8_cols,:);  end
% % figure 3
for j = 1 : M9_cols   U14(:,:,j) = OUTPUT_X((j-1)*M9_cols+1 : j*M9_cols,:);    end
for j = 1 : M10_cols  U15(:,:,j) = OUTPUT_Y((j-1)*M10_cols+1 : j*M10_cols,:);  end
for j = 1 : M11_cols  U16(:,:,j) = OUTPUT_Z((j-1)*M11_cols+1 : j*M11_cols,:);  end
% figure 4
for j = 1 : M12_cols  U17(:,:,j) =  OUTPUT_RSQ((j-1)*M12_cols+1 : j*M12_cols,:);    end
for j = 1 : M13_cols  U18(:,:,j) =  OUTPUT_THETA((j-1)*M13_cols+1 : j*M13_cols,:);  end
for j = 1 : M15_cols  U20(:,:,j) =  OUTPUT_PER((j-1)*M15_cols+1 : j*M15_cols,:);    end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%  UNIT CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)

% Step size
dx = Sx/(M1_cols-1); % x-axis step size
dy = Sy/(M1_cols-1); % y-axis step size
dz = Sz/(M1_cols-1); % z-axis step size

x = 0 : dx : 1;
y = 0 : dy : 1;
z = 0 : dz : 1;

[X, Y, Z] = meshgrid(x, y, z);

  %subplot(1,1,1);
  [T, p] = isosurface(X, Y, Z, U1, .5);
  pa = patch('Faces',T,'Vertices',p,'FaceVertexCData',p,'FaceColor','interp', 'EdgeColor', 'none');
  %pa = patch('Faces',T,'Vertices',p); # draws the surface with standard settings
  %set(pa,'EdgeColor', 'none')
  %set(pa,'FaceColor','interp');
  %set(gca,'fontsize',14)
  %axis equal tight; 
  %view(-30, 30)
  view(3)
  xlabel('x','FontSize',18);
  ylabel('y','FontSize',18');
  zlabel('z','FontSize',18,'Rotation',0);
  set(gca,'FontSize',18);
  title('3D UNIT CELL','FontSize',18)
  axis vis3d tight
  camlight left
  colormap('parula');
  lighting gouraud
  grid on
  box on
  ax = gca;
  ax.BoxStyle = 'full';
  
%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GRATING
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step size
dx = Sx/(M6_cols-1); % x-axis step size
dy = Sy/(M6_cols-1); % y-axis step size
dz = Sz/(M6_cols-1); % z-axis step size

x = [-37.70   -31.42   -25.13   -18.85   -12.57   -6.28    0.00    6.28   12.57   18.85   25.13   31.42   37.70 ];
y = [-37.70   -31.42   -25.13   -18.85   -12.57   -6.28    0.00    6.28   12.57   18.85   25.13   31.42   37.70 ];
z = [-37.70   -31.42   -25.13   -18.85   -12.57   -6.28    0.00    6.28   12.57   18.85   25.13   31.42   37.70 ];

[X, Y, Z] = meshgrid(x, y, z);

figure(2) 
 subplot(1,3,1);
 slice(X,Y,Z,U13,[ ],[ ],37.70);  
 hold on
 slice(X,Y,Z,U13,[],-37.70,[ ]); 
 slice(X,Y,Z,U13,-37.70,[],[ ]); 
 set(gca,'fontsize',14)
 axis equal tight; 
 axis([-37.70 37.70 -37.70 37.70 -37.70 37.70])
 view(-30, 30)
 xlabel('x','FontSize',18);
 ylabel('y','FontSize',18');
 zlabel('z','FontSize',18,'Rotation',0);
 set(gca,'FontSize',18);
 title('Kx','FontSize',18)
 grid on
 
 subplot(1,3,2);
 slice(X,Y,Z,U12,[ ],[ ],37.70);  
 hold on
 slice(X,Y,Z,U12,[],-37.70,[ ]); 
 slice(X,Y,Z,U12,-37.70,[],[ ]); 
 set(gca,'fontsize',14)
 axis equal tight; 
 axis([-37.70 37.70 -37.70 37.70 -37.70 37.70])
 view(-30, 30)
 xlabel('x','FontSize',18);
 ylabel('y','FontSize',18');
 zlabel('z','FontSize',18,'Rotation',0);
 set(gca,'FontSize',18);
 title('Ky','FontSize',18)
 grid on
 
 subplot(1,3,3);
 slice(X,Y,Z,U11,[ ],[ ],37.70);  
 hold on
 slice(X,Y,Z,U11,[],-37.70,[ ]); 
 slice(X,Y,Z,U11,-37.70,[],[ ]); 
 set(gca,'fontsize',14)
 axis equal tight; 
 axis([-37.70 37.70 -37.70 37.70 -37.70 37.70])
 view(-30, 30)
 xlabel('x','FontSize',18);
 ylabel('y','FontSize',18');
 zlabel('z','FontSize',18,'Rotation',0);
 set(gca,'FontSize',18);
 title('Kz','FontSize',18)
 grid on
 %colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%
%  XYZ Array
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step size
[m n p] = size(U16);

dx = Sx/(M9_cols-1); % x-axis step size
dy = Sy/(M9_cols-1); % y-axis step size
dz = Sz/(M9_cols-1); % z-axis step size

x = 0 : U16(m,2,1) : U16(1,m,1);
y = 0 : U16(m,2,1) : U16(1,m,1);
z = 0 : U16(m,2,1) : U16(1,m,1);

[X, Y, Z] = meshgrid(x, y, z);

figure(3)
 subplot(1,3,1);
 slice(X,Y,Z,U16,[ ],[ ],U16(1,m,1));  
 hold on
 slice(X,Y,Z,U16,[],0,[ ]); 
 slice(X,Y,Z,U16,0,[],[ ]); 
 axis equal tight; 
 axis([0 U16(1,m,1) 0 U16(1,m,1) 0 U16(1,m,1)])
 view(-30, 30)
 xlabel('x','FontSize',18);
 ylabel('y','FontSize',18');
 zlabel('z','FontSize',18,'Rotation',0);
 set(gca,'FontSize',18);
 title('X','FontSize',18)
 grid on
 
 subplot(1,3,2);
 slice(X,Y,Z,U15,[ ],[ ],U16(1,m,1));  
 hold on
 slice(X,Y,Z,U15,[],0,[ ]); 
 slice(X,Y,Z,U15,0,[],[ ]); 
 axis equal tight; 
 axis([0 U16(1,m,1) 0 U16(1,m,1) 0 U16(1,m,1)])
 view(-30, 30)
 xlabel('x','FontSize',18);
 ylabel('y','FontSize',18');
 zlabel('z','FontSize',18,'Rotation',0);
 set(gca,'FontSize',18);
 title('Y','FontSize',18)
 grid on
 
 subplot(1,3,3);
 slice(X,Y,Z,U14,[ ],[ ],U16(1,m,1));  
 hold on
 slice(X,Y,Z,U14,[],0,[ ]); 
 slice(X,Y,Z,U14,0,[],[ ]); 
 axis equal tight; 
 axis([0 U16(1,m,1) 0 U16(1,m,1) 0 U16(1,m,1)])
 view(-30, 30)
 xlabel('x','FontSize',18);
 ylabel('y','FontSize',18');
 zlabel('z','FontSize',18,'Rotation',0);
 set(gca,'FontSize',18);
 title('Z','FontSize',18)
 grid on
 %colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%  RSQ, THETA, PER
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step size
dx = Sx/(M12_cols-1); % x-axis step size
dy = Sy/(M12_cols-1); % y-axis step size
dz = Sz/(M12_cols-1); % z-axis step size

x = 0 : dx : 1;
y = 0 : dy : 1;
z = 0 : dz : 1;

[X, Y, Z] = meshgrid(x, y, z);

 figure(4)
 subplot(1,3,1);
 slice(X,Y,Z,U17,[ ],[ ],1);  
 hold on
 slice(X,Y,Z,U17,[],0,[ ]); 
 slice(X,Y,Z,U17,0,[],[ ]); 
 axis equal; 
 view(-30, 30)
 xlabel('x','FontSize',18);
 ylabel('y','FontSize',18');
 zlabel('z','FontSize',18,'Rotation',0);
 set(gca,'FontSize',18);
 title('RSQ','FontSize',18)
 grid on

 subplot(1,3,2);
 slice(X,Y,Z,U18,[ ],[ ],1);  
 hold on
 slice(X,Y,Z,U18,[],0,[ ]); 
 slice(X,Y,Z,U18,0,[],[ ]); 
 axis equal; 
 view(-30, 30)
 xlabel('x');
 xlabel('x','FontSize',18);
 ylabel('y','FontSize',18');
 zlabel('z','FontSize',18,'Rotation',0);
 set(gca,'FontSize',18);;
 title('THETA')
 grid on

 subplot(1,3,3);
 slice(X,Y,Z,U20,[ ],[ ],1);  
 hold on
 slice(X,Y,Z,U20,[],0,[ ]); 
 slice(X,Y,Z,U20,0,[],[ ]); 
 axis equal; 
 view(-30, 30)
 xlabel('x','FontSize',18);
 ylabel('y','FontSize',18');
 zlabel('z','FontSize',18,'Rotation',0);
 set(gca,'FontSize',18);
 title('PER','FontSize',18)
 grid on
 %colorbar
 
 pause
 close all
 