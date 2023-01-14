import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from skimage import measure

def pause():
    program_Pause = raw_input("Press the <ENTER> key to continue...\n")
    plt.cla()   # Clear axis
    plt.clf()   # Clear figure
    plt.close() # Close a figure window
    exit()      # exit program

###################################################################
#                         PATHS
###################################################################
# Path Data Files
Filename_1  = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_UNIT_CELL_REAL"
Filename_2  = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_UNIT_CELL_IMAG"
Filename_3  = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_FFTW_REAL"
Filename_4  = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_FFTW_IMAG"
Filename_5  = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_INV_FFTW_REAL"
Filename_6  = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_INV_FFTW_IMAG"
Filename_7  = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_SWAP_REAL"
Filename_8  = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_SWAP_IMAG"

Filename_9  = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_TRUNC_FFTW_REAL"
Filename_10 = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_TRUNC_FFTW_IMAG"
Filename_11 = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_KX"
Filename_12 = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_KY"
Filename_13 = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_KZ"

Filename_14 = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_RSQ"
Filename_15 = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_THETA"
Filename_16 = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_VARPHI"
Filename_17 = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_PER"

Filename_18 = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_PHI.mat"
Filename_19 = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_S.mat"
Filename_20 = "/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_UC.mat"

###################################################################
#                         LOADS
###################################################################
# Load Data Files
U1 = np.loadtxt(Filename_1)
U2 = np.loadtxt(Filename_2)
U3 = np.loadtxt(Filename_3)
U4 = np.loadtxt(Filename_4)
U5 = np.loadtxt(Filename_5)
U6 = np.loadtxt(Filename_6)
U7 = np.loadtxt(Filename_7)
U8 = np.loadtxt(Filename_8)

U9  = np.loadtxt(Filename_9)
U10 = np.loadtxt(Filename_10)
U11 = np.loadtxt(Filename_11)
U12 = np.loadtxt(Filename_12)
U13 = np.loadtxt(Filename_13)

U14 = np.loadtxt(Filename_14)
U15 = np.loadtxt(Filename_15)
U16 = np.loadtxt(Filename_16)
U17 = np.loadtxt(Filename_17)

# Load files
U18 = open(Filename_18, "r")
U19 = open(Filename_19, "r")
U20 = open(Filename_20, "r")

###################################################################
#                         FIGURES
###################################################################
# figure 1
print '\nfigure 1 : number of elemets = ', np.size(U1)
print "2D Array size, (x, y) = ", U1.shape
rows = U1.shape[0]
columns = U1.shape[1]
# Reshape 2D data file to 3D
U1 = U1.reshape((columns, columns, columns))
U2 = U2.reshape((columns, columns, columns))
U3 = U3.reshape((columns, columns, columns))
U4 = U4.reshape((columns, columns, columns))
U5 = U5.reshape((columns, columns, columns))
U6 = U6.reshape((columns, columns, columns))
U7 = U7.reshape((columns, columns, columns))
U8 = U8.reshape((columns, columns, columns))
print "3D Array size, (x, y, z) = ", U8.shape

# figure 2
print '\nfigure 2 : number of elemets size = ', np.size(U9)
print "2D Array size, U9(x,y) = ", U9.shape
rows = U9.shape[0]
columns = U9.shape[1]
# Reshape 2D data file to 3D
U9  = U9.reshape((columns, columns, columns))
U10 = U10.reshape((columns, columns, columns))
U11 = U11.reshape((columns, columns, columns))
U12 = U12.reshape((columns, columns, columns))
U13 = U13.reshape((columns, columns, columns))
print "3D Array size, (x, y, z) = ", U13.shape

# figure 3
print '\nfigure 3 : number of elemets size = ', np.size(U14)
print "2D Array size, U14(x,y) = ", U14.shape
rows = U14.shape[0]
columns = U14.shape[1]
# Reshape 2D data file to 3D
U14 = U14.reshape((columns, columns, columns))
U15 = U15.reshape((columns, columns, columns))
U16 = U16.reshape((columns, columns, columns))
U17 = U17.reshape((columns, columns, columns))
print "3D Array size, (x, y, z) = ", U16.shape

"""
# figure 4
print '\nfigure 4 : number of elemets size = ', np.size(U18)
print "2D Array size, U17(x,y) = ", U18.shape
rows = U18.shape[0]
columns = U18.shape[1]
U18 = U18.reshape((columns, columns, columns))
U19 = U19.reshape((columns, columns, columns))
U20 = U20.reshape((columns, columns, columns))
print "3D Array size, (x, y, z) = ", U20.shape
"""
# 3D PLOT
# figure 1
fig = plt.figure(1)
ax = fig.add_subplot(221, projection='3d')
verts, faces = measure.marching_cubes(U1, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('3D UNIT CELL REAL',fontsize = 10)
ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

ax = fig.add_subplot(222, projection='3d')
verts, faces = measure.marching_cubes(U3, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('REAL FFTW',fontsize = 10)
ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

ax = fig.add_subplot(223, projection='3d')
verts, faces = measure.marching_cubes(U5, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('REAL INV FFTW',fontsize = 10)
ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

ax = fig.add_subplot(224, projection='3d')
verts, faces = measure.marching_cubes(U7, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('REAL SWAP',fontsize = 10)
ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

# Save Figure
plt.tight_layout(0.5) 
plt.savefig("figure1.eps", bbox_inches='tight')

# figure 2
fig = plt.figure(2)
ax = fig.add_subplot(231, projection='3d')
verts, faces = measure.marching_cubes(U9, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('OUTPUT_TRUNC_FFTW_REAL',fontsize = 10)
ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

ax = fig.add_subplot(232, projection='3d')
verts, faces = measure.marching_cubes(U11, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('KX',fontsize = 10)
#ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

ax = fig.add_subplot(233, projection='3d')
verts, faces = measure.marching_cubes(U12, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('KY',fontsize = 10)
#ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

ax = fig.add_subplot(234, projection='3d')
verts, faces = measure.marching_cubes(U13, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('KZ',fontsize = 10)
#ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

# Save Figure
plt.tight_layout(0.5) 
plt.savefig("figure2.eps", bbox_inches='tight')

# figure 3
"""
fig = plt.figure(3)
ax = fig.add_subplot(231, projection='3d')
verts, faces = measure.marching_cubes(U14, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('RSQ',fontsize = 10)
ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

ax = fig.add_subplot(232, projection='3d')
verts, faces = measure.marching_cubes(U15, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('THETA',fontsize = 10)
ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

ax = fig.add_subplot(233, projection='3d')
verts, faces = measure.marching_cubes(U16, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('VARPHI',fontsize = 10)
ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

ax = fig.add_subplot(234, projection='3d')
verts, faces = measure.marching_cubes(U17, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('OUTPUT_PER',fontsize = 10)
#ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

# Save Figure
plt.tight_layout(0.5) 
plt.savefig("figure3.eps", bbox_inches='tight')

"""

# File size
num_lines_1 = sum(1 for line in open(Filename_18))
num_lines_2 = sum(1 for line in open(Filename_19))
num_lines_3 = sum(1 for line in open(Filename_20))

# Initialize arrays
PHI = []
S = []
UC = []


# Build 1D arrays
i = 0
for line in U18:
  words1 = line.split(" ")
  i=i+1
  if i > 2 :       #skip rows
     PHI.append(float(words1[0]))#PHI



i = 0
for line in U19:
  words2 = line.split(" ")
  i=i+1
  if i > 2 :       #skip rows
     S.append(float(words2[0]))# S

i = 0
for line in U20:
  words3 = line.split(" ")
  i=i+1
  if i > 2 :       #skip rows
     UC.append(float(words3[0]))# U

print 'number of elemets size = ', num_lines_1, num_lines_2, num_lines_3 
print 'number of elemets size (PHI, S, UC) = ', len(PHI), len(S), len(UC)

# Reshape 1D array to 3D arrays
sqrt1 = np.cbrt(len(PHI))
PHIS = np.reshape(PHI, (int(sqrt1), int(sqrt1), int(sqrt1)))

sqrt2 = np.cbrt(len(S))
SS = np.reshape(S, (int(sqrt2), int(sqrt2), int(sqrt2)))

sqrt3 = np.cbrt(len(UC))
UCS = np.reshape(UC, (int(sqrt3), int(sqrt3), int(sqrt3)))

print "Array size, PHI(x,y,z) = ", PHIS.shape
print "Array size,   S(x,y,z) = ", SS.shape
print "Array size,  UC(x,y,z) = ", UCS.shape

# figure 4
fig = plt.figure(4)
ax = fig.add_subplot(131, projection='3d')
verts, faces = measure.marching_cubes(PHIS, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('OUTPUT_PHI',fontsize = 10)
ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

ax = fig.add_subplot(132, projection='3d')
verts, faces = measure.marching_cubes(SS, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('OUTPUT_S',fontsize = 10)
ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

ax = fig.add_subplot(133, projection='3d')
verts, faces = measure.marching_cubes(UCS, 0, spacing=(0.1, 0.1, 0.1))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('OUTPUT_UC',fontsize = 10)
ax.set_aspect('equal')
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

# Save Figure
plt.tight_layout(0.5) 
plt.savefig("figure4.eps", bbox_inches='tight')

# figure 5
fig = plt.figure(5)
ax = fig.add_subplot(111, projection='3d')
verts, faces = measure.marching_cubes(UCS, 0, spacing=(0.05, 0.05, 0.05))
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_title('OUTPUT_UC',fontsize = 10)
ax.set_aspect('equal')
ax.view_init(-0, 0) # elevation, azimuth
ax.set_xlabel('X', fontsize = 10)
ax.set_ylabel('Y', fontsize = 10)
ax.set_zlabel('Z', fontsize = 10)
ax.grid(color='k', linestyle='-', linewidth=1)

# Save Figure
plt.tight_layout(0.5) 
plt.savefig("figure5.eps", bbox_inches='tight')


# Show plots
plt.show()


