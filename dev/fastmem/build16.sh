source /opt/Modules/default/init/bash
module use /usr/common/software/carl_modulefiles
module load intel/16.0.3.210 vtune sde impi
module load memkind

export SCALAPACK_LIBS="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64"
export LAPACK_LIBS="-mkl=parallel"
export BLAS_LIBS="-mkl=parallel"
export FFT_LIBS="-mkl=parallel"
export AR=ar
export MPIF90=mpiifort
export MPICC=mpiicc
export MPIF77=mpiifort
export MPICXX=mpiicc
export I_MPI_F77=ifort
export I_MPI_CXX=icpc
export I_MPI_F90=ifort
export I_MPI_CC=icc
export F90=ifort
export F77=ifort
export CXX=icpc

make -j 8 pw
