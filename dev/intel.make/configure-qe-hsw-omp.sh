#!/bin/bash

TARGET="-xCORE-AVX2"
OMPFLAG="-qopenmp -qoverride_limits"
#IPO="-ipo-separate"
OPTC=-O3
OPTF=-O2
PRFX=2017-

HERE=$(cd $(dirname $0); pwd -P)
export ELPAROOT="${HERE}/../elpa/${PRFX}hsw-omp"
#export MKLRTL="sequential"
export MKLRTL="intel_thread"
export OPENMP="--enable-openmp"
export LD_LIBS="-Wl,--as-needed -liomp5 -Wl,--no-as-needed"
export MPIF90=mpiifort
export CC=mpiicc
export AR=xiar
export dir=none

#LIBXSMM="-Wl,--wrap=sgemm_,--wrap=dgemm_ ${HOME}/libxsmm/lib/libxsmmext.a ${HOME}/libxsmm/lib/libxsmm.a"
export BLAS_LIBS="${LIBXSMM} -Wl,--start-group \
    ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
    ${MKLROOT}/lib/intel64/libmkl_core.a \
    ${MKLROOT}/lib/intel64/libmkl_${MKLRTL}.a \
    ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
  -Wl,--end-group"
export LAPACK_LIBS="${BLAS_LIBS}"
export SCALAPACK_LIBS="${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a"
#export SCALAPACK_LIBS="${HOME}/scalapack-2.0.2/libscalapack.a"
export FFT_LIBS="${BLAS_LIBS}"

./configure ${OPENMP} --with-elpa=${ELPAROOT} --with-scalapack=intel $*

# adjust generated configuration
sed -i \
  -e "s/-nomodule -openmp/-nomodule/" \
  -e "s/-par-report0 -vec-report0//" \
  -e "s/-D__FFTW3/-D__DFTI/" \
  make.sys
sed -i \
  -e "s/-D__FFTW/-D__DFTI/" -e "s/-D__ELPA/-D__ELPA3/" \
  -e "s/IFLAGS         = -I\.\.\/include/IFLAGS         = -I\.\.\/include -I\$(MKLROOT)\/include\/fftw/" \
  -e "s/-O3/${OPTC} ${IPO} ${TARGET} -fno-alias -ansi-alias/" \
  -e "s/-O2 -assume byterecl -g -traceback/${OPTF} -align array64byte -threads -heap-arrays 4096 ${IPO} ${TARGET} -assume byterecl/" \
  -e "s/LDFLAGS        =/LDFLAGS        = -static-intel -static-libgcc -static-libstdc++/" \
  -e "s/-openmp/${OMPFLAG}/" \
  make.sys

# Uncomment below block in case of compiler issue (ICE)
echo >> make.sys
cat configure-qe-tbbmalloc.mak >> make.sys
echo -e "init_us_1.o: init_us_1.f90\n\t\$(MPIF90) \$(F90FLAGS) -O1 -c \$<\n" >> make.sys
echo -e "new_ns.o: new_ns.f90\n\t\$(MPIF90) \$(F90FLAGS) -O1 -c \$<\n" >> make.sys
echo -e "us_exx.o: us_exx.f90\n\t\$(MPIF90) \$(F90FLAGS) ${OMPFLAG} -c \$<\n" >> make.sys

#sed -i -e "s/\$(MOD_FLAG)\.\.\/ELPA\/src/\$(MOD_FLAG)${ELPAROOT}\/src/" Modules/Makefile
sed -i -e "s/\$(MOD_FLAG)\.\.\/ELPA\/src//" Modules/Makefile

# patch source code files
patch -N PW/src/setup.f90 configure-qe-setup_pw.patch
if [ -e LAXlib/dspev_drv.f90 ]; then
  patch -N LAXlib/dspev_drv.f90 configure-qe-dspev_drv.patch
else
  patch -N Modules/dspev_drv.f90 configure-qe-dspev_drv.patch
fi
if [ -e LAXlib/zhpev_drv.f90 ]; then
  patch -N LAXlib/zhpev_drv.f90 configure-qe-zhpev_drv.patch
else
  patch -N Modules/zhpev_drv.f90 configure-qe-zhpev_drv.patch
fi
#patch -N Modules/wavefunctions.f90 configure-qe-wavefunctions.patch
#patch -N FFTXlib/fft_parallel.f90 configure-qe-fft_parallel.patch
patch -N FFTXlib/fftw.c configure-qe-fftw.patch

# reminder
echo "Ready to \"make pw\"!"

