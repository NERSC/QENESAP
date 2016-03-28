!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------!
! FFT scalar drivers Module - contains machine-dependent routines for      !
! FFTW, FFTW3, ESSL (both 3d for serial execution and 1d+2d FFTs for       !
! parallel execution; NEC ASL libraries (3d only, no parallel execution)   !
! Written by Carlo Cavazzoni, modified by P. Giannozzi, contributions      !
! by Martin Hilgemans, Guido Roma, Pascal Thibaudeau, Stephane Lefranc,    !
! Nicolas Lacorne, Filippo Spiga, Nicola Varini - Last update Jul 2015     !
!--------------------------------------------------------------------------!


#if defined(__FFTW3)

!=----------------------------------------------------------------------=!
   MODULE fft_scalar
!=----------------------------------------------------------------------=!

       USE, intrinsic ::  iso_c_binding
       
       IMPLICIT NONE
       SAVE

       PRIVATE
       PUBLIC :: cft_1z, cft_2xy, cfft3d, cfft3ds

! ...   Local Parameter

#include "fft_param.f90"

#if defined(__OPENMP)
#include "fftw3.f03"
#else
#include "fftw3.f"
#endif

!=----------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------=!

!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "z"
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cft_1z(c, nsl, nz, ldz, isign, cout, is_exx)

!     driver routine for nsl 1d complex fft's of length nz
!     ldz >= nz is the distance between sequences to be transformed
!     (ldz>nz is used on some architectures to reduce memory conflicts)
!     input  :  c(ldz*nsl)   (complex)
!     output : cout(ldz*nsl) (complex - NOTA BENE: transform is not in-place!)
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nz, nsl, ldz) are stored and re-used if available

     INTEGER, INTENT(IN) :: isign
     INTEGER, INTENT(IN) :: nsl, nz, ldz
     LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
     LOGICAL :: is_exx_

     COMPLEX (DP) :: c(:), cout(:)

     REAL (DP)  :: tscale
     INTEGER    :: i, err, idir, ip, void
     INTEGER :: zdims( 3, ndims ) = -1
     INTEGER :: icurrent = 1
     INTEGER, SAVE :: zdims_local( 3, ndims ) = -1
     INTEGER, SAVE :: icurrent_local = 1
     INTEGER, SAVE :: zdims_exx( 3, ndims ) = -1
     INTEGER, SAVE :: icurrent_exx = 1
     LOGICAL :: done
     INTEGER :: tid

#if defined(__OPENMP)
     INTEGER :: offset, ldz_t
     INTEGER :: omp_get_max_threads
     EXTERNAL :: omp_get_max_threads
#endif

     !   Pointers to the "C" structures containing FFT factors ( PLAN )
     !   C_POINTER is defined in include/fft_defs.h
     !   for 32bit executables, C_POINTER is integer(4)
     !   for 64bit executables, C_POINTER is integer(8)

     C_POINTER :: fw_planz( ndims ) = 0
     C_POINTER :: bw_planz( ndims ) = 0
     C_POINTER, SAVE :: fw_planz_local( ndims ) = 0
     C_POINTER, SAVE :: bw_planz_local( ndims ) = 0
     C_POINTER, SAVE :: fw_planz_exx( ndims ) = 0
     C_POINTER, SAVE :: bw_planz_exx( ndims ) = 0

     IF(PRESENT(is_exx))THEN
        is_exx_ = is_exx
     ELSE
        is_exx_ = .FALSE.
     END IF
     IF(is_exx_)THEN
        zdims = zdims_exx
        icurrent = icurrent_exx
        fw_planz = fw_planz_exx
        bw_planz = bw_planz_exx
     ELSE
        zdims = zdims_local
        icurrent = icurrent_local
        fw_planz = fw_planz_local
        bw_planz = bw_planz_local
     END IF

     IF( nsl < 0 ) THEN
       CALL fftx_error__(" fft_scalar: cft_1z ", " nsl out of range ", nsl)
     END IF

     !
     !   Here initialize table only if necessary
     !
     
     CALL lookup()

     IF( .NOT. done ) THEN

       !   no table exist for these parameters
       !   initialize a new one
      
       CALL init_plan()

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'cft_1z' )
#endif

     IF (isign < 0) THEN
        CALL dfftw_execute_dft( fw_planz( ip), c, cout)
        tscale = 1.0_DP / nz
        cout( 1 : ldz * nsl ) = cout( 1 : ldz * nsl ) * tscale
     ELSE IF (isign > 0) THEN
        CALL dfftw_execute_dft( bw_planz( ip), c, cout)
     END IF

     IF(is_exx_)THEN
        zdims_exx = zdims
        icurrent_exx = icurrent
        fw_planz_exx = fw_planz
        bw_planz_exx = bw_planz
     ELSE
        zdims_local = zdims
        icurrent_local = icurrent
        fw_planz_local = fw_planz
        bw_planz_local = bw_planz
     END IF

#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'cft_1z' )
#endif

     RETURN

   CONTAINS

     SUBROUTINE lookup()
        ! lookup for stored plan 
        DO ip = 1, ndims
           !   first check if there is already a table initialized
           !   for this combination of parameters
           !   The initialization in ESSL and FFTW v.3 depends on all three parameters
           done = ( nz == zdims(1,ip) )
           done = done .AND. ( nsl == zdims(2,ip) ) .AND. ( ldz == zdims(3,ip) )
           IF (done) EXIT
        END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
#if defined(__OPENMP)
       CALL dfftw_cleanup_threads() 
       void = fftw_init_threads()
       CALL dfftw_plan_with_nthreads(omp_get_max_threads())      
#endif

       IF( fw_planz( icurrent) /= 0 ) CALL dfftw_destroy_plan( fw_planz( icurrent) )
       IF( bw_planz( icurrent) /= 0 ) CALL dfftw_destroy_plan( bw_planz( icurrent) )
       idir = -1
       CALL dfftw_plan_many_dft( fw_planz( icurrent), 1, nz, nsl, c, &
            (/SIZE(c)/), 1, ldz, cout, (/SIZE(cout)/), 1, ldz, idir, FFTW_ESTIMATE)
       idir = 1
       CALL dfftw_plan_many_dft( bw_planz( icurrent), 1, nz, nsl, c, &
            (/SIZE(c)/), 1, ldz, cout, (/SIZE(cout)/), 1, ldz, idir, FFTW_ESTIMATE)

       zdims(1,icurrent) = nz; zdims(2,icurrent) = nsl; zdims(3,icurrent) = ldz;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

   END SUBROUTINE cft_1z

!
!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "x" and "y" direction
!
!
!
!=----------------------------------------------------------------------=!
!
!

   SUBROUTINE cft_2xy(r, nzl, nx, ny, ldx, ldy, isign, pl2ix, is_exx)

!     driver routine for nzl 2d complex fft's of lengths nx and ny
!     input : r(ldx*ldy)  complex, transform is in-place
!     ldx >= nx, ldy >= ny are the physical dimensions of the equivalent
!     2d array: r2d(ldx, ldy) (x first dimension, y second dimension)
!     (ldx>nx, ldy>ny used on some architectures to reduce memory conflicts)
!     pl2ix(nx) (optional) is 1 for columns along y to be transformed
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nx,ny,nzl,ldx) are stored and re-used if available

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: isign, ldx, ldy, nx, ny, nzl
     INTEGER, OPTIONAL, INTENT(IN) :: pl2ix(:)
     LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
     LOGICAL :: is_exx_
     COMPLEX (DP) :: r( : )
     INTEGER :: i, k, j, err, idir, ip, kk, void
     REAL(DP) :: tscale
     INTEGER :: icurrent = 1
     INTEGER :: dims( 4, ndims) = -1
     INTEGER, SAVE :: icurrent_local
     INTEGER, SAVE :: dims_local( 4, ndims) = -1
     INTEGER, SAVE :: icurrent_exx
     INTEGER, SAVE :: dims_exx( 4, ndims) = -1
     LOGICAL :: dofft( nfftx ), done
     INTEGER, PARAMETER  :: stdout = 6

#if defined __HPM
     INTEGER :: OMP_GET_THREAD_NUM
#endif
#if defined(__OPENMP)
     INTEGER :: offset
     INTEGER :: nx_t, ny_t, nzl_t, ldx_t, ldy_t
     INTEGER  :: itid, mytid, ntids
     INTEGER  :: omp_get_thread_num, omp_get_num_threads,omp_get_max_threads
     EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

     C_POINTER :: fw_plan( 2, ndims ) = 0
     C_POINTER :: bw_plan( 2, ndims ) = 0
     C_POINTER, SAVE :: fw_plan_local( 2, ndims ) = 0
     C_POINTER, SAVE :: bw_plan_local( 2, ndims ) = 0
     C_POINTER, SAVE :: fw_plan_exx( 2, ndims ) = 0
     C_POINTER, SAVE :: bw_plan_exx( 2, ndims ) = 0

     IF(PRESENT(is_exx))THEN
        is_exx_ = is_exx
     ELSE
        is_exx_ = .FALSE.
     END IF
     IF(is_exx_)THEN
        dims = dims_exx
        icurrent = icurrent_exx
        fw_plan = fw_plan_exx
        bw_plan = bw_plan_exx
     ELSE
        dims = dims_local
        icurrent = icurrent_local
        fw_plan = fw_plan_local
        bw_plan = bw_plan_local
     END IF

     dofft( 1 : nx ) = .TRUE.
     IF( PRESENT( pl2ix ) ) THEN
       IF( SIZE( pl2ix ) < nx ) &
         CALL fftx_error__( ' cft_2xy ', ' wrong dimension for arg no. 8 ', 1 )
       DO i = 1, nx
         IF( pl2ix(i) < 1 ) dofft( i ) = .FALSE.
       END DO
     END IF

     !
     !   Here initialize table only if necessary
     !
 
     CALL lookup()

     IF( .NOT. done ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       CALL init_plan()

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'cft_2xy' )
#endif

     IF ( ldx /= nx .OR. ldy /= ny ) THEN
        IF( isign < 0 ) THEN
           do j = 0, nzl-1
              CALL dfftw_execute_dft( fw_plan (1, ip), &
                   r(1+j*ldx*ldy:), r(1+j*ldx*ldy:))
           end do
           do i = 1, nx
              do k = 1, nzl
                 IF( dofft( i ) ) THEN
                    j = i + ldx*ldy * ( k - 1 )
                    call dfftw_execute_dft( fw_plan ( 2, ip), r(j:), r(j:))
                 END IF
              end do
           end do
           tscale = 1.0_DP / ( nx * ny )
           CALL ZDSCAL( ldx * ldy * nzl, tscale, r(1), 1)
        ELSE IF( isign > 0 ) THEN
           do i = 1, nx
              do k = 1, nzl
                 IF( dofft( i ) ) THEN
                    j = i + ldx*ldy * ( k - 1 )
                    call dfftw_execute_dft( bw_plan ( 2, ip), r(j:), r(j:))
                 END IF
              end do
           end do
           do j = 0, nzl-1
              CALL dfftw_execute_dft( bw_plan( 1, ip), &
                   r(1+j*ldx*ldy:), r(1+j*ldx*ldy:))
           end do
        END IF
     ELSE
        IF( isign < 0 ) THEN
           call dfftw_execute_dft( fw_plan( 1, ip), r(1:), r(1:))
           tscale = 1.0_DP / ( nx * ny )
           CALL ZDSCAL( ldx * ldy * nzl, tscale, r(1), 1)
        ELSE IF( isign > 0 ) THEN
           call dfftw_execute_dft( bw_plan( 1, ip), r(1:), r(1:))
        END IF
     END IF

     IF(is_exx_)THEN
        dims_exx = dims
        icurrent_exx = icurrent
        fw_plan_exx = fw_plan
        bw_plan_exx = bw_plan
     ELSE
        dims_local = dims
        icurrent_local = icurrent
        fw_plan_local = fw_plan
        bw_plan_local = bw_plan
     END IF

#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'cft_2xy' )
#endif

     RETURN

  CONTAINS

     SUBROUTINE lookup()
       DO ip = 1, ndims
         !   first check if there is already a table initialized
         !   for this combination of parameters
         done = ( ny == dims(1,ip) ) .AND. ( nx == dims(3,ip) )
         done = done .AND. ( ldx == dims(2,ip) ) .AND.  ( nzl == dims(4,ip) )
         IF (done) EXIT
       END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()

#if defined(__OPENMP)
       CALL dfftw_cleanup_threads() 
       void = fftw_init_threads()
       CALL dfftw_plan_with_nthreads(omp_get_max_threads())      
#endif

       IF ( ldx /= nx .OR. ldy /= ny ) THEN
          IF( fw_plan(2,icurrent) /= 0 )  CALL dfftw_destroy_plan( fw_plan(2,icurrent) )
          IF( bw_plan(2,icurrent) /= 0 )  CALL dfftw_destroy_plan( bw_plan(2,icurrent) )
          idir = -1
          CALL dfftw_plan_many_dft( fw_plan(2,icurrent), 1, ny, 1, r(1:), &
               (/ldx*ldy/), ldx, 1, r(1:), (/ldx*ldy/), ldx, 1, idir, &
               FFTW_ESTIMATE)
          idir =  1
          CALL dfftw_plan_many_dft( bw_plan(2,icurrent), 1, ny, 1, r(1:), &
               (/ldx*ldy/), ldx, 1, r(1:), (/ldx*ldy/), ldx, 1, idir, &
               FFTW_ESTIMATE)

          IF( fw_plan(1,icurrent) /= 0 ) CALL dfftw_destroy_plan( fw_plan(1,icurrent) )
          IF( bw_plan(1,icurrent) /= 0 ) CALL dfftw_destroy_plan( bw_plan(1,icurrent) )
          idir = -1
          CALL dfftw_plan_many_dft( fw_plan(1,icurrent), 1, nx, ny, r(1:), &
               (/ldx*ldy/), 1, ldx, r(1:), (/ldx*ldy/), 1, ldx, idir, &
               FFTW_ESTIMATE)
          idir =  1
          CALL dfftw_plan_many_dft( bw_plan(1,icurrent), 1, nx, ny, r(1:), &
               (/ldx*ldy/), 1, ldx, r(1:), (/ldx*ldy/), 1, ldx, idir, &
               FFTW_ESTIMATE)
       ELSE
          IF( fw_plan( 1, icurrent) /= 0 ) CALL dfftw_destroy_plan( fw_plan( 1, icurrent) )
          IF( bw_plan( 1, icurrent) /= 0 ) CALL dfftw_destroy_plan( bw_plan( 1, icurrent) )
          idir = -1
          CALL dfftw_plan_many_dft( fw_plan( 1, icurrent), 2, (/nx, ny/), nzl,&
               r(1:), (/nx, ny/), 1, nx*ny, r(1:), (/nx, ny/), 1, nx*ny, idir,&
               FFTW_ESTIMATE)
          idir = 1
          CALL dfftw_plan_many_dft( bw_plan( 1, icurrent), 2, (/nx, ny/), nzl,&
               r(1:), (/nx, ny/), 1, nx*ny, r(1:), (/nx, ny/), 1, nx*ny, idir,&
               FFTW_ESTIMATE)
       END IF

       dims(1,icurrent) = ny; dims(2,icurrent) = ldx;
       dims(3,icurrent) = nx; dims(4,icurrent) = nzl;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

   END SUBROUTINE cft_2xy


!
!=----------------------------------------------------------------------=!
!
!
!
!         3D scalar FFTs
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cfft3d( fftbox, nx, ny, nz, ldx, ldy, ldz, isign, is_exx )

  !     driver routine for 3d complex fft of lengths nx, ny, nz
  !     input  :  f(ldx*ldy*ldz)  complex, transform is in-place
  !     ldx >= nx, ldy >= ny, ldz >= nz are the physical dimensions
  !     of the equivalent 3d array: f3d(ldx,ldy,ldz)
  !     (ldx>nx, ldy>ny, ldz>nz may be used on some architectures
  !      to reduce memory conflicts - not implemented for FFTW)
  !     isign > 0 : f(G) => f(R)   ; isign < 0 : f(R) => f(G)
  !
  !     Up to "ndims" initializations (for different combinations of input
  !     parameters nx,ny,nz) are stored and re-used if available

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nx, ny, nz, ldx, ldy, ldz, isign
     LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
     LOGICAL :: is_exx_
     COMPLEX (DP) :: fftbox(nx,ny,nz)
     INTEGER :: i,ibeg,jbeg, k, j, err, idir, ip
     REAL(DP) :: tscale
     INTEGER :: icurrent = 1
     INTEGER :: dims(3,ndims) = -1
     INTEGER, SAVE :: icurrent_local = 1
     INTEGER, SAVE :: dims_local(3,ndims) = -1
     INTEGER, SAVE :: icurrent_exx = 1
     INTEGER, SAVE :: dims_exx(3,ndims) = -1
     integer,parameter::nblk=8

     C_POINTER :: fw_plan(ndims) = 0
     C_POINTER :: bw_plan(ndims) = 0
     C_POINTER, save :: fw_plan_local(ndims) = 0
     C_POINTER, save :: bw_plan_local(ndims) = 0
     C_POINTER, save :: fw_plan_exx(ndims) = 0
     C_POINTER, save :: bw_plan_exx(ndims) = 0
     integer(8),dimension(ndims),save:: planFMX,planF1X,planFMY,planF1Y,planFMZ,planF1Z
     integer(8),dimension(ndims),save:: planBMX,planB1X,planBMY,planB1Y,planBMZ,planB1Z

     IF(PRESENT(is_exx))THEN
        is_exx_ = is_exx
     ELSE
        is_exx_ = .FALSE.
     END IF
     IF(is_exx_)THEN
        dims = dims_exx
        icurrent = icurrent_exx
        fw_plan = fw_plan_exx
        bw_plan = bw_plan_exx
     ELSE
        dims = dims_local
        icurrent = icurrent_local
        fw_plan = fw_plan_local
        bw_plan = bw_plan_local
     END IF

     IF ( nx < 1 ) &
         call fftx_error__('cfft3d',' nx is less than 1 ', 1)
     IF ( ny < 1 ) &
         call fftx_error__('cfft3d',' ny is less than 1 ', 1)
     IF ( nz < 1 ) &
         call fftx_error__('cfft3',' nz is less than 1 ', 1)

     !
     !   Here initialize table only if necessary
     !
     CALL lookup()

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       CALL init_plan()

     END IF

     !
     !   Now perform the 3D FFT using the machine specific driver
     !

     IF( isign < 0 ) THEN
        tscale = 1.0_DP / DBLE( nx * ny * nz )
        !$omp parallel  private(k,jbeg,j,i,ibeg)
        !$omp do 
        do j=1,ny
          do ibeg=1,nx,nblk
            if(ibeg+nblk-1<=nx)then
              call dfftw_execute_dft(planFMZ(ip),fftbox(ibeg,j,1),fftbox(ibeg,j,1))
            else
              do i=ibeg,nx
                call dfftw_execute_dft(planF1Z(ip),fftbox(i,j,1),fftbox(i,j,1))
              enddo
            endif
          enddo
        enddo
        !$omp barrier  ! barrier not needed here as the next fft will be a different one
        !$omp do 
        do k=1,nz
          do ibeg=1,nx,nblk
            if(ibeg+nblk-1<=nx)then
              call dfftw_execute_dft(planFMY(ip),fftbox(ibeg,1,k),fftbox(ibeg,1,k))
            else
              do i=ibeg,nx
                call dfftw_execute_dft(planF1Y(ip),fftbox(i,1,k),fftbox(i,1,k))
              enddo
            endif
          enddo
          do jbeg=1,ny,nblk
            if(jbeg+nblk-1<=ny)then
              call dfftw_execute_dft(planFMX(ip),fftbox(1,jbeg,k),fftbox(1,jbeg,k))
            else
              do j=jbeg,ny
                call dfftw_execute_dft(planF1X(ip),fftbox(1,j,k),fftbox(1,j,k))
              enddo
            endif
          enddo
          call ZDSCAL( nx * ny , tscale, fftbox(1,1,k), 1)
        enddo
        !$omp end parallel 

     ELSE IF( isign > 0 ) THEN

        !$omp parallel  private(k,jbeg,j,i,ibeg)
        !$omp do 
        do j=1,ny
          do ibeg=1,nx,nblk
            if(ibeg+nblk-1<=nx)then
              call dfftw_execute_dft(planBMZ(ip),fftbox(ibeg,j,1),fftbox(ibeg,j,1))
            else
              do i=ibeg,nx
                call dfftw_execute_dft(planB1Z(ip),fftbox(i,j,1),fftbox(i,j,1))
              enddo
            endif
          enddo
        enddo
        !$omp barrier  ! barrier not needed here as the next fft will be a different one
        !$omp do 
        do k=1,nz
          do ibeg=1,nx,nblk
            if(ibeg+nblk-1<=nx)then
              call dfftw_execute_dft(planBMY(ip),fftbox(ibeg,1,k),fftbox(ibeg,1,k))
            else
              do i=ibeg,nx
                call dfftw_execute_dft(planB1Y(ip),fftbox(i,1,k),fftbox(i,1,k))
              enddo
            endif
          enddo
          do jbeg=1,ny,nblk
            if(jbeg+nblk-1<=ny)then
              call dfftw_execute_dft(planBMX(ip),fftbox(1,jbeg,k),fftbox(1,jbeg,k))
            else
              do j=jbeg,ny
                call dfftw_execute_dft(planB1X(ip),fftbox(1,j,k),fftbox(1,j,k))
              enddo
            endif
          enddo
        enddo
        !$omp end parallel 

     END IF

     IF(is_exx_)THEN
        dims_exx = dims
        icurrent_exx = icurrent
        fw_plan_exx = fw_plan
        bw_plan_exx = bw_plan
     ELSE
        dims_local = dims
        icurrent_local = icurrent
        fw_plan_local = fw_plan
        bw_plan_local = bw_plan
     END IF

     RETURN

   CONTAINS

     SUBROUTINE lookup()
     ip = -1
     DO i = 1, ndims
       !   first check if there is already a table initialized
       !   for this combination of parameters
       IF ( ( nx == dims(1,i) ) .and. &
            ( ny == dims(2,i) ) .and. &
            ( nz == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
       IF ( nx /= ldx .or. ny /= ldy .or. nz /= ldz ) &
            call fftx_error__('cfft3','not implemented',3)
       IF( fw_plan(icurrent) /= 0 ) CALL dfftw_destroy_plan( fw_plan(icurrent) )
       IF( bw_plan(icurrent) /= 0 ) CALL dfftw_destroy_plan( bw_plan(icurrent) )
       IF( planFMX(icurrent) /= 0 ) CALL dfftw_destroy_plan( planFMX(icurrent) )
       IF( planF1X(icurrent) /= 0 ) CALL dfftw_destroy_plan( planF1X(icurrent) )
       IF( planFMY(icurrent) /= 0 ) CALL dfftw_destroy_plan( planFMY(icurrent) )
       IF( planF1Y(icurrent) /= 0 ) CALL dfftw_destroy_plan( planF1Y(icurrent) )
       IF( planFMZ(icurrent) /= 0 ) CALL dfftw_destroy_plan( planFMZ(icurrent) )
       IF( planF1Z(icurrent) /= 0 ) CALL dfftw_destroy_plan( planF1Z(icurrent) )
       IF( planBMX(icurrent) /= 0 ) CALL dfftw_destroy_plan( planBMX(icurrent) )
       IF( planB1X(icurrent) /= 0 ) CALL dfftw_destroy_plan( planB1X(icurrent) )
       IF( planBMY(icurrent) /= 0 ) CALL dfftw_destroy_plan( planBMY(icurrent) )
       IF( planB1Y(icurrent) /= 0 ) CALL dfftw_destroy_plan( planB1Y(icurrent) )
       IF( planBMZ(icurrent) /= 0 ) CALL dfftw_destroy_plan( planBMZ(icurrent) )
       IF( planB1Z(icurrent) /= 0 ) CALL dfftw_destroy_plan( planB1Z(icurrent) )
       idir = -1
       CALL dfftw_plan_dft_3d ( fw_plan(icurrent), nx, ny, nz, fftbox(1,1,1), &
            fftbox(1,1,1), idir, FFTW_ESTIMATE)
       idir =  1
       CALL dfftw_plan_dft_3d ( bw_plan(icurrent), nx, ny, nz, fftbox(1,1,1), &
            fftbox(1,1,1), idir, FFTW_ESTIMATE)
        call dfftw_plan_many_dft(planFMX(icurrent),1,nx,nblk, &  ! plan,int rank, const int *n, int howmany,
                                fftbox,nx,      &  ! fftw_complex *in, const int *inembed,
                                1,nx,           &  ! int istride, int idist,
                                fftbox,nx,      &  ! fftw_complex *out, const int *onembed,
                                1,nx,           &  ! int ostride, int odist,
                                FFTW_FORWARD,FFTW_ESTIMATE+FFTW_UNALIGNED )  ! int sign, unsigned flags);
        call dfftw_plan_many_dft(planF1X(icurrent),1,nx,1, &  ! plan,int rank, const int *n, int howmany,
                                fftbox,nx,      &  ! fftw_complex *in, const int *inembed,
                                1,nx,           &  ! int istride, int idist,
                                fftbox,nx,      &  ! fftw_complex *out, const int *onembed,
                                1,nx,           &  ! int ostride, int odist,
                                FFTW_FORWARD,FFTW_ESTIMATE+FFTW_UNALIGNED )  ! int sign, unsigned flags);
        call dfftw_plan_many_dft(planFMY(icurrent),1,ny,nblk, &
                                fftbox,ny,      &
                                nx,1,           &
                                fftbox,ny,      &
                                nx,1,           &
                                FFTW_FORWARD,FFTW_ESTIMATE+FFTW_UNALIGNED )
        call dfftw_plan_many_dft(planF1Y(icurrent),1,ny,1, &
                                fftbox,ny,      &
                                nx,1,           &
                                fftbox,ny,      &
                                nx,1,           &
                                FFTW_FORWARD,FFTW_ESTIMATE+FFTW_UNALIGNED )
        call dfftw_plan_many_dft(planFMZ(icurrent),1,nz,nblk, &
                                fftbox,nz,      &
                                nx*ny,1,           &
                                fftbox,nz,      &
                                nx*ny,1,           &
                                FFTW_FORWARD,FFTW_ESTIMATE+FFTW_UNALIGNED )
        call dfftw_plan_many_dft(planF1Z(icurrent),1,nz,1, &
                                fftbox,nz,      &
                                nx*ny,1,           &
                                fftbox,nz,      &
                                nx*ny,1,           &
                                FFTW_FORWARD,FFTW_ESTIMATE+FFTW_UNALIGNED )

        call dfftw_plan_many_dft(planBMX(icurrent),1,nx,nblk, &  ! plan,int rank, const int *n, int howmany,
                                fftbox,nx,      &  ! fftw_complex *in, const int *inembed,
                                1,nx,           &  ! int istride, int idist,
                                fftbox,nx,      &  ! fftw_complex *out, const int *onembed,
                                1,nx,           &  ! int ostride, int odist,
                                FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_UNALIGNED )  ! int sign, unsigned flags);
        call dfftw_plan_many_dft(planB1X(icurrent),1,nx,1, &  ! plan,int rank, const int *n, int howmany,
                                fftbox,nx,      &  ! fftw_complex *in, const int *inembed,
                                1,nx,           &  ! int istride, int idist,
                                fftbox,nx,      &  ! fftw_complex *out, const int *onembed,
                                1,nx,           &  ! int ostride, int odist,
                                FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_UNALIGNED )  ! int sign, unsigned flags);
        call dfftw_plan_many_dft(planBMY(icurrent),1,ny,nblk, &
                                fftbox,ny,      &
                                nx,1,           &
                                fftbox,ny,      &
                                nx,1,           &
                                FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_UNALIGNED )
        call dfftw_plan_many_dft(planB1Y(icurrent),1,ny,1, &
                                fftbox,ny,      &
                                nx,1,           &
                                fftbox,ny,      &
                                nx,1,           &
                                FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_UNALIGNED )
        call dfftw_plan_many_dft(planBMZ(icurrent),1,nz,nblk, &
                                fftbox,nz,      &
                                nx*ny,1,           &
                                fftbox,nz,      &
                                nx*ny,1,           &
                                FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_UNALIGNED )
        call dfftw_plan_many_dft(planB1Z(icurrent),1,nz,1, &
                                fftbox,nz,      &
                                nx*ny,1,           &
                                fftbox,nz,      &
                                nx*ny,1,           &
                                FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_UNALIGNED )

       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

   END SUBROUTINE cfft3d

!
!=----------------------------------------------------------------------=!
!
!
!
!         3D scalar FFTs,  but using sticks!
!
!
!
!=----------------------------------------------------------------------=!
!

SUBROUTINE cfft3ds (f, nx, ny, nz, ldx, ldy, ldz, isign, &
     do_fft_x, do_fft_y, is_exx)
  !
  !     driver routine for 3d complex "reduced" fft - see cfft3d
  !     The 3D fft are computed only on lines and planes which have
  !     non zero elements. These lines and planes are defined by
  !     the two integer vectors do_fft_x(ldy*nz) and do_fft_y(nz)
  !     (1 = perform fft, 0 = do not perform fft)
  !     This routine is implemented only for fftw, essl, acml
  !     If not implemented, cfft3d is called instead
  !
  !----------------------------------------------------------------------
  !
  implicit none

  integer :: nx, ny, nz, ldx, ldy, ldz, isign
  !
  !   logical dimensions of the fft
  !   physical dimensions of the f array
  !   sign of the transformation

  complex(DP) :: f ( ldx * ldy * ldz )
  integer :: do_fft_x(:), do_fft_y(:)
  LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
  LOGICAL :: is_exx_
  !
  integer :: m, incx1, incx2
  INTEGER :: i, k, j, err, idir, ip,  ii, jj
  REAL(DP) :: tscale
  INTEGER :: icurrent = 1
  INTEGER :: dims(3,ndims) = -1
  INTEGER, SAVE :: icurrent_local = 1
  INTEGER, SAVE :: dims_local(3,ndims) = -1
  INTEGER, SAVE :: icurrent_exx = 1
  INTEGER, SAVE :: dims_exx(3,ndims) = -1


  C_POINTER :: fw_plan ( 3, ndims ) = 0
  C_POINTER :: bw_plan ( 3, ndims ) = 0
  C_POINTER, SAVE :: fw_plan_local( 3, ndims ) = 0
  C_POINTER, SAVE :: bw_plan_local( 3, ndims ) = 0
  C_POINTER, SAVE :: fw_plan_exx( 3, ndims ) = 0
  C_POINTER, SAVE :: bw_plan_exx( 3, ndims ) = 0

  IF(PRESENT(is_exx))THEN
     is_exx_ = is_exx
  ELSE
     is_exx_ = .FALSE.
  END IF
  IF(is_exx_)THEN
     dims = dims_exx
     icurrent = icurrent_exx
     fw_plan = fw_plan_exx
     bw_plan = bw_plan_exx
  ELSE
     dims = dims_local
     icurrent = icurrent_local
     fw_plan = fw_plan_local
     bw_plan = bw_plan_local
  END IF

  tscale = 1.0_DP

  ! WRITE( stdout, fmt="('DEBUG cfft3ds :',6I6)") nx, ny, nz, ldx, ldy, ldz
  ! WRITE( stdout, fmt="('DEBUG cfft3ds :',24I2)") do_fft_x
  ! WRITE( stdout, fmt="('DEBUG cfft3ds :',24I2)") do_fft_y


     IF( ny /= ldy ) &
       CALL fftx_error__(' cfft3ds ', ' wrong dimensions: ny /= ldy ', 1 )

     CALL lookup()

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       CALL init_plan()

     END IF


     IF ( isign > 0 ) THEN

        !
        !  i - direction ...
        !

        incx1 = 1;  incx2 = ldx;  m = 1

        !$omp parallel do private(k,j,jj,ii) schedule(dynamic,1)
        do k = 1, nz
           do j = 1, ny
              jj = j + ( k - 1 ) * ldy
              ii = 1 + ldx * ( jj - 1 )
              if ( do_fft_x( jj ) == 1 ) THEN
                call dfftw_execute_dft( bw_plan( 1, ip), f( ii: ), f( ii: ) )
              endif
           enddo
        enddo

        !
        !  ... j-direction ...
        !

        incx1 = ldx;  incx2 = 1;  m = nx

        !$omp parallel do private(k,ii) schedule(dynamic,1)
        do k = 1, nz
           ii = 1 + ldx * ldy * ( k - 1 )
           if ( do_fft_y( k ) == 1 ) then
             call dfftw_execute_dft( bw_plan( 2, ip), f( ii: ), f( ii: ) )
           endif
        enddo

        !
        !     ... k-direction
        !

        incx1 = ldx * ldy;  incx2 = 1;  m = ldx * ny

        
        !$omp parallel do private(j,ii)
        do j = 1, ny
          ii = 1 + ldx * ( j - 1 )
          call dfftw_execute_dft( bw_plan( 3, ip), f(ii:), f(ii:) )
        enddo

     ELSE

        !
        !     ... k-direction
        !

        incx1 = ldx * ny;  incx2 = 1;  m = ldx * ny

        !$omp parallel do private(j,ii)
        do j = 1, ny
          ii = 1 + ldx * ( j - 1 )
          call dfftw_execute_dft( fw_plan( 3, ip), f(ii:), f(ii:) )
        enddo

        !
        !     ... j-direction ...
        !

        incx1 = ldx;  incx2 = 1;  m = nx

        !$omp parallel do private(k,ii) schedule(dynamic,1)
        do k = 1, nz
           ii = 1 + ldx * ldy * ( k - 1 )
           if ( do_fft_y ( k ) == 1 ) then
             call dfftw_execute_dft( fw_plan( 2, ip), f( ii: ), f( ii: ) )
           endif
        enddo

        !
        !     i - direction ...
        !

        incx1 = 1;  incx2 = ldx;  m = 1

        !$omp parallel do private(k,j,jj,ii) schedule(dynamic,1)
        do k = 1, nz
           do j = 1, ny
              jj = j + ( k - 1 ) * ldy
              ii = 1 + ldx * ( jj - 1 )
              if ( do_fft_x( jj ) == 1 ) then
                call dfftw_execute_dft( fw_plan( 1, ip), f( ii: ), f( ii: ) )
              endif
           enddo
        enddo
        
!         call DSCAL (2 * ldx * ldy * nz, 1.0_DP/(nx * ny * nz), f(1), 1)
        
        !$omp parallel do private(k)
        do k = 1,nz
          call DSCAL (2 * ldx * ldy , 1.0_DP/(nx * ny * nz), f(1+(k-1)*ldx*ldy), 1)
        enddo

     END IF

     IF(is_exx_)THEN
        dims_exx = dims
        icurrent_exx = icurrent
        fw_plan_exx = fw_plan
        bw_plan_exx = bw_plan
     ELSE
        dims_local = dims
        icurrent_local = icurrent
        fw_plan_local = fw_plan
        bw_plan_local = bw_plan
     END IF

     RETURN

   CONTAINS


     SUBROUTINE lookup()
     ip = -1
     DO i = 1, ndims
       !   first check if there is already a table initialized
       !   for this combination of parameters
       IF( ( nx == dims(1,i) ) .and. ( ny == dims(2,i) ) .and. &
           ( nz == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
       IF( fw_plan( 1, icurrent) /= 0 ) &
            CALL dfftw_destroy_plan( fw_plan( 1, icurrent) )
       IF( bw_plan( 1, icurrent) /= 0 ) &
            CALL dfftw_destroy_plan( bw_plan( 1, icurrent) )
       IF( fw_plan( 2, icurrent) /= 0 ) &
            CALL dfftw_destroy_plan( fw_plan( 2, icurrent) )
       IF( bw_plan( 2, icurrent) /= 0 ) &
            CALL dfftw_destroy_plan( bw_plan( 2, icurrent) )
       IF( fw_plan( 3, icurrent) /= 0 ) &
            CALL dfftw_destroy_plan( fw_plan( 3, icurrent) )
       IF( bw_plan( 3, icurrent) /= 0 ) &
            CALL dfftw_destroy_plan( bw_plan( 3, icurrent) )
       idir = -1
       CALL dfftw_plan_many_dft( fw_plan( 1, icurrent), &
            1, nx, 1, f(1:), (/ldx, ldy, ldz/), 1, ldx, &
            f(1:), (/ldx, ldy, ldz/), 1, ldx, idir, FFTW_ESTIMATE)
       idir = 1
       CALL dfftw_plan_many_dft( bw_plan( 1, icurrent), &
            1, nx, 1, f(1:), (/ldx, ldy, ldz/), 1, ldx, &
            f(1:), (/ldx, ldy, ldz/), 1, ldx, idir, FFTW_ESTIMATE)
       idir = -1
       CALL dfftw_plan_many_dft( fw_plan( 2, icurrent), &
            1, ny, nx, f(1:), (/ldx, ldy, ldz/), ldx, 1, &
            f(1:), (/ldx, ldy, ldz/), ldx, 1, idir, FFTW_ESTIMATE)
       idir = 1
       CALL dfftw_plan_many_dft( bw_plan( 2, icurrent), &
            1, ny, nx, f(1:), (/ldx, ldy, ldz/), ldx, 1, &
            f(1:), (/ldx, ldy, ldz/), ldx, 1, idir, FFTW_ESTIMATE)
       idir = -1
       CALL dfftw_plan_many_dft( fw_plan( 3, icurrent), &
            1, nz, nx, f(1:), (/ldx, ldy, ldz/), ldx*ldy, 1, &
            f(1:), (/ldx, ldy, ldz/), ldx*ldy, 1, idir, FFTW_ESTIMATE)
       idir = 1
       CALL dfftw_plan_many_dft( bw_plan( 3, icurrent), &
            1, nz, nx, f(1:), (/ldx, ldy, ldz/), ldx*ldy, 1, &
            f(1:), (/ldx, ldy, ldz/), ldx*ldy, 1, idir, FFTW_ESTIMATE)

       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

   END SUBROUTINE cfft3ds

!=----------------------------------------------------------------------=!
   END MODULE fft_scalar
!=----------------------------------------------------------------------=!

#endif
