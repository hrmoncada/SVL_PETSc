#!/bin/sh
clear all

## petsc 3.5.4
#export PETSC_DIR=~/Desktop/PETSC/petsc-3.5.4
#export PETSC_ARCH=linux-gnu-complex #PetscScalar is Complex

# petsc 3.7.6
export PETSC_DIR=~/Desktop/PETSC/petsc-3.7.6
export PETSC_ARCH=linux-gnu-complex #PetscScalar is Complex

make SVL_PETSC

#octave PLOT_OUTPUT_3D_METHOD_1.m
#octave PLOT_OUTPUT_3D_METHOD_2.m
