### Load modules:
# $ module load petsc/3.7-complex
# $ module load fftw3
### Compile:
# $ make SVL_PETSC
### Submit:
# s batch mujob.sh
### Clean folder:
# $ bash SVL_Clean.sh
### Multiples submitions:
# $ cp myjob.sh OUTPUT_PETSC Clean.sh  SVL.data OUT_Node_Cores_Nx_GridSize_GridResolution_CompilationNumber/ 
### Make multiples copies of the folder OUT_Node_Cores_Nx_GridSize_GridResolution_CompilationNumber/
# $ cp -R OUT...C1 OUT..C2

CFLAGS           = -I/opt/apps/intel17/impi17_0/fftw3/3.3.6/include
FFLAGS           =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           =                                      # Working folder
EXAMPLESC        =                                      # *.c or *.cpp file names here
EXAMPLESF        =
MANSEC           =
FFTW_LIBS        = -L/opt/apps/intel17/impi17_0/fftw3/3.3.6/lib -lfftw3     -L\${TACC_FFTW3_LIB} # fftw serial libraries
FFTW_MPI         = -L/opt/apps/intel17/impi17_0/fftw3/3.3.6/lib -lfftw3_mpi -lfftw3 -L\${TACC_FFTW3_LIB} # fftw mpi libraries
#vi MPI_LIBS         = -L\${TACC_MPIP_LIB} -lmpiP     # mpip librarie no longer on stampede2
MATH_lIBS        = -lm        # math libraries
TRASH            :=  *.*~  *~  *.o
# Note: PETSC main -> SVL_PETSC_DFT_MAIN.c

SOURCE_PETSC     := SVL_PETSC_3D_DFT_MAIN_ORNL.c\
		    SVL_PETSC_3D_UNIT_CELL.c\
		    SVL_PETSC_3D_FFTW_SWAP_METHOD_1.c\
                    SVL_PETSC_3D_FFTW_SWAP_METHOD_2.c\
		    SVL_PETSC_1D_FFTW.c\
		    SVL_PETSC_3D_FFTW.c\
		    SVL_PETSC_3D_SWAP_QUADRANTS.c\
		    SVL_PETSC_3D_TRANSPOSE_1_COMPLEX.c\
		    SVL_PETSC_3D_TRANSPOSE_2_COMPLEX.c\
		    SVL_PETSC_3D_TRANSPOSE_3_COMPLEX.c\
		    SVL_PETSC_3D_TRANSPOSE_4_COMPLEX.c\
		    SVL_PETSC_3D_IFFTW.c\
		    SVL_PETSC_3D_TRUNCATED_FFTW_ARRAY.c\
 		    SVL_PETSC_3D_GRADING_VECTOR.c\
                    SVL_PETSC_3D_IMPLEMENT_IMPROVEMENTS.c\
		    SVL_PETSC_3D_ELIMINATE_GRATING_ACCORD_THEIR_AMPLITUD.c\
		    SVL_PETSC_3D_IDENTIFIED_COLLINEAR_PLANAR_GRATING.c\
                    SVL_PETSC_3D_SPHERICAL_SPATIAL_VARIANT.c\
                    SVL_PETSC_3D_SPHERICAL_ORIENTATION_FUNCTION.c\
                    SVL_PETSC_3D_SPHERICAL_LATTICE_SPACING_FUNCTION.c\
                    SVL_PETSC_3D_CYLINDRICAL_SPATIAL_VARIANT.c\
                    SVL_PETSC_3D_CYLINDRICAL_ORIENTATION_FUNCTION.c\
                    SVL_PETSC_3D_CYLINDRICAL_LATTICE_SPACING_FUNCTION.c\
                    SVL_PETSC_3D_IMPROVEMENTS.c\
		    SVL_PETSC_3D_FDDER.c\
		    SVL_PETSC_3D_LOOP.c\
		    SVL_PETSC_3D_ORIENTATION_VECTOR.c\
                    SVL_PETSC_3D_SPHERICAL_TRANSLATION.c\
		    SVL_PETSC_3D_CARTESIAN_TO_SPHERICAL.c\
		    SVL_PETSC_3D_SPHERICAL_ROTATION.c\
	            SVL_PETSC_3D_SPHERICAL_SPACING.c\
		    SVL_PETSC_3D_SPHERICAL_TO_CARTESIAN.c\
                    SVL_PETSC_3D_CYLINDRICAL_TRANSLATION.c\
		    SVL_PETSC_3D_CARTESIAN_TO_CYLINDRICAL.c\
		    SVL_PETSC_3D_CYLINDRICAL_ROTATION.c\
	            SVL_PETSC_3D_CYLINDRICAL_SPACING.c\
		    SVL_PETSC_3D_CYLINDRICAL_TO_CARTESIAN.c\
                    SVL_PETSC_3D_RHS.c\
                    SAVE_1D_TO_3D_ARRAY_REAL.c\
		    SAVE_1D_TO_3D_ARRAY_COMPLEX.c\
                                        
OBJECTS_PETSC    := $(SOURCE_PETSC:.c=.o)
EXECUTABLE_PETSC := OUTPUT_PETSC

# Version 3.5.4: These Makefiles lines must be updated every time you Update your petsc version.
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

# Version 3.7.6: These Makefiles lines must be updated every time you Update your petsc version.
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# $ icc -xMIC-AVX512  -O3 mycode.c -o myexe # will run only on KNL
SVL_PETSC: $: $(OBJECTS_PETSC) chkopts
	   -${CLINKER}  -xMIC-AVX512 -O1 -g -o $(EXECUTABLE_PETSC) $(OBJECTS_PETSC) ${FFTW_LIBS} ${MATH_LIBS} ${PETSC_LIB}  # ${MPI_LIBS} 
	    ${RM} $(OBJECTS_PETSC) $(TRASH)

