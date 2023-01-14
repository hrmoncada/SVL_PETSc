# SVL_PETSc
Spatially Variant Lattice Algorithm (SVL) with FFTW (http://www.fftw.org/) and  PETSC (https://petsc.org/release/)

We wrote a portable computer program for parallel architectures with general purpose programming language that supports structured programming. For the parallel code, we use the FFTW (Fastest Fourier Transform in the West) for handling the Fourier transform of the unit cell device and PETSc (Portable, Extensible Toolkit for Scientific Computation) for handling the numerical linear algebra operations. Using Message Passing Interface (MPI) for distributed memory helps us to improve the performance of the code that generates 2D and 3D SVL when it is executed on a parallel system. The SVL code was executed the on Stampede 2 supercomputers at the Texas Advanced Computing Center (TACC), the University of Texas at Austin.

The SVL code can be executed on a personal computer if PETSC and FFTW are properly installed. Also, we are showing the code we execute on Stampede 2 on KNL and SKX architectures

For more detail on the code implementation, please review the chapter my MS and Ph.D. thesis at the following links:
* MS-CPS Thesis:  See chapter 3 for more details on the code (https://scholarworks.utep.edu/dissertations/AAI10118248/)
* PhD-CPS Thesis: See chapter 2 for more details on the cod (https://scholarworks.utep.edu/dissertations/AAI10824303/)
