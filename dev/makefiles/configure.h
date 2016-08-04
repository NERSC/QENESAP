!
! Copyright (C) 2006 WanT Group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!
! contains configure infos
!

#define  __HAVE_CONFIG_INFO

#define   __CONF_HOST           "@host@"
#define   __CONF_ARCH           "x86_64"

#define   __CONF_CC             "cc"
#define   __CONF_CFLAGS         "-g"
#define   __CONF_DFLAGS         " -D__INTEL -D__FFTW3 -D__OPENMP -D__MPI -D__PARA -D__SCALAPACK"
#define   __CONF_CPP            "cpp"
#define   __CONF_CPPFLAGS       "-P -C -traditional"

#define   __CONF_F90            "ftn"
#define   __CONF_MPIF90         "ftn"
#define   __CONF_F90FLAGS       "$(FFLAGS) -openmp -nomodule -ffp -g"
#define   __CONF_F77            "ftn"
#define   __CONF_FFLAGS         "-O2 -assume byterecl -g -traceback"
#define   __CONF_FFLAGS_NOOPT   "-O0 -assume byterecl -g -traceback"
#define   __CONF_PRE_FDFLAGS    "-fpp "
#define   __CONF_FDFLAGS        "$(DFLAGS) $(MANUAL_DFLAGS)"

#define   __CONF_LD             "ftn"
#define   __CONF_LDFLAGS        "-openmp -static-intel "
#define   __CONF_IMOD           "-I"

#define   __CONF_BLAS_LIBS      "  -lmkl_intel_lp64  -lmkl_sequential -lmkl_core"
#define   __CONF_LAPACK_LIBS    ""
#define   __CONF_FFT_LIBS       ""
#define   __CONF_MASS_LIBS      ""

#define   __CONF_AR             "ar"
#define   __CONF_ARFLAGS        "ruv"
#define   __CONF_RANLIB         "ranlib"


