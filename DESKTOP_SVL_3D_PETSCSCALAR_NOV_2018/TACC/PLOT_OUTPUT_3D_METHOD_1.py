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

# Path Data Files
Filename_1 = np.loadtxt('/home/henry/Desktop/PETSC/Examples/SVL_3D_V_4_1_PETSC_DESKTOP/PETSCSALAR/SVL_3D_PETSCSCALAR/OUTPUT_UNIT_CELL_REAL')

print "Filename_1 = ", Filename_1
print '\n'
print 'number of elemets = ', np.size(Filename_1)
print "2D Array size, (x, y) = ", Filename_1.shape
rows = Filename_1.shape[0]
columns = Filename_1.shape[1]
print "2D Array size, (x, y) = ", (rows, columns)

# Reshape data file to 3D
U1 = Filename_1.reshape((columns, columns, columns))
print "Filename_1 = ", U1
print '\n'
print 'number of elemets size = ', np.size(U1)
print "3D Array size, U1(x,y,z) = ", U1.shape
rows_0 = U1.shape[0]
rows_1 = U1.shape[1]
rows_2 = U1.shape[2]
print "3D Array size, U1(x,y,z) = ", (rows_0, rows_1, rows_2)

# 3D Mesh 
X = np.arange(0, rows_0)
Y = np.arange(0, rows_0)
Z = np.arange(0, rows_0)
X, Y, Z = np.meshgrid(X, Y, Z)
print "3D Array size, X(x,y,z) = ", X.shape
print "3D Array size, Y(x,y,z) = ", Y.shape
print "3D Array size, Z(x,y,z) = ", Z.shape

# 3D PLOT
fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
#ax = plt.axes(projection='3d')
ax.scatter(X, Y, Z, c = U1, edgecolors='none');
#ax.scatter(X, Y, Z, c = U1, cmap='viridis', linewidth=0.5);
#ax.scatter(X, Y, Z, c = U1, cmap='cubehelix', linewidth=0.5);
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

verts, faces = measure.marching_cubes(U1, 0, spacing=(0.1, 0.1, 0.1))
ax = fig.add_subplot(122, projection='3d')
#ax = plt.axes(projection='3d')
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],cmap ='Spectral', lw=1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()


