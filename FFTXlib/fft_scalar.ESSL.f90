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

#if defined(__ESSL) || defined (__LINUX_ESSL)

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

        !   Workspace that is statically allocated is defined here
        !   in order to avoid multiple copies of the same workspace
        !   lwork:   Dimension of the work space array (if any)

        !   ESSL IBM library: see the ESSL manual for DCFT

        INTEGER, PARAMETER :: lwork = 20000 + ( 2*nfftx + 256 ) * 64 + 3*nfftx
        REAL (DP) :: work( lwork )


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
     !INTEGER, SAVE :: zdims( 3, ndims ) = -1
     !INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: zdims_local( 3, ndims ) = -1
     INTEGER, SAVE :: icurrent_local = 1
     INTEGER, SAVE :: zdims_exx( 3, ndims ) = -1
     INTEGER, SAVE :: icurrent_exx = 1
     LOGICAL :: done

#if defined __HPM
     INTEGER :: OMP_GET_THREAD_NUM
#endif
     INTEGER :: tid

     ! ...   Machine-Dependent parameters, work arrays and tables of factors

     !   ltabl   Dimension of the tables of factors calculated at the
     !           initialization stage

#if defined(__OPENMP)
     INTEGER :: offset, ldz_t
     INTEGER :: omp_get_max_threads
     EXTERNAL :: omp_get_max_threads
#endif


     !   ESSL IBM library: see the ESSL manual for DCFT

     INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
     !REAL (DP), SAVE :: fw_tablez( ltabl, ndims )
     !REAL (DP), SAVE :: bw_tablez( ltabl, ndims )
     REAL (DP), SAVE :: fw_tablez_local( ltabl, ndims )
     REAL (DP), SAVE :: bw_tablez_local( ltabl, ndims )
     REAL (DP), SAVE :: fw_tablez_exx( ltabl, ndims )
     REAL (DP), SAVE :: bw_tablez_exx( ltabl, ndims )

     IF(PRESENT(is_exx))THEN
        is_exx_ = is_exx
     ELSE
        is_exx_ = .FALSE.
     END IF

     IF( nsl < 0 ) THEN
       CALL fftx_error__(" fft_scalar: cft_1z ", " nsl out of range ", nsl)
     END IF

     !
     !   Here initialize table only if necessary
     !

     DO ip = 1, ndims

        !   first check if there is already a table initialized
        !   for this combination of parameters

        !   The initialization in ESSL and FFTW v.3 depends on all three parameters
        IF(is_exx_)THEN
           done = ( nz == zdims_exx(1,ip) )
           
           done = done .AND. ( nsl == zdims_exx(2,ip) ) .AND. ( ldz == zdims_exx(3,ip) )
        ELSE
           done = ( nz == zdims_local(1,ip) )
           
           done = done .AND. ( nsl == zdims_local(2,ip) ) .AND. ( ldz == zdims_local(3,ip) )
        IF (done) EXIT
     END DO

     IF( .NOT. done ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       ! WRITE( stdout, fmt="('DEBUG cft_1z, reinitializing tables ', I3)" ) icurrent

       tscale = 1.0_DP / nz

       IF(is_exx_)THEN
          CALL DCFT ( 1, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl,  1, &
               tscale, fw_tablez_exx(1, icurrent_exx), ltabl, work(1), lwork)
          CALL DCFT ( 1, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, -1, &
               1.0_DP, bw_tablez_exx(1, icurrent_exx), ltabl, work(1), lwork)

          zdims_exx(1,icurrent_exx) = nz
          zdims_exx(2,icurrent_exx) = nsl
          zdims_exx(3,icurrent_exx) = ldz
          ip = icurrent_exx
          icurrent_exx = MOD( icurrent_exx, ndims ) + 1
       ELSE
          CALL DCFT ( 1, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl,  1, &
               tscale, fw_tablez_local(1, icurrent_local), ltabl, work(1), lwork)
          CALL DCFT ( 1, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, -1, &
               1.0_DP, bw_tablez_local(1, icurrent_local), ltabl, work(1), lwork)

          zdims_local(1,icurrent_local) = nz
          zdims_local(2,icurrent_local) = nsl
          zdims_local(3,icurrent_local) = ldz
          ip = icurrent_local
          icurrent_local = MOD( icurrent_local, ndims ) + 1
       END IF

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'cft_1z' )
#endif


     ! essl uses a different convention for forward/backward transforms
     ! wrt most other implementations: notice the sign of "idir"

     IF( isign < 0 ) THEN
        idir   =+1
        tscale = 1.0_DP / nz
        IF(is_exx_)THEN
           CALL DCFT (0, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, idir, &
                tscale, fw_tablez_exx(1, ip), ltabl, work, lwork)
        ELSE
           CALL DCFT (0, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, idir, &
                tscale, fw_tablez_local(1, ip), ltabl, work, lwork)
        END IF
     ELSE IF( isign > 0 ) THEN
        idir   =-1
        tscale = 1.0_DP
        IF(is_exx_)THEN
           CALL DCFT (0, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, idir, &
                tscale, bw_tablez_exx(1, ip), ltabl, work, lwork)
        ELSE
           CALL DCFT (0, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, idir, &
                tscale, bw_tablez_local(1, ip), ltabl, work, lwork)
        END IF
     END IF

#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'cft_1z' )
#endif

     RETURN

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
     !INTEGER, SAVE :: icurrent = 1
     !INTEGER, SAVE :: dims( 4, ndims) = -1
     INTEGER, SAVE :: icurrent_local = 1
     INTEGER, SAVE :: dims_local( 4, ndims) = -1
     INTEGER, SAVE :: icurrent_exx = 1
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

     INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
     !REAL (DP), SAVE :: fw_tablex( ltabl, ndims ), fw_tabley( ltabl, ndims )
     !REAL (DP), SAVE :: bw_tablex( ltabl, ndims ), bw_tabley( ltabl, ndims )
     REAL (DP), SAVE :: fw_tablex_local( ltabl, ndims ), fw_tabley_local( ltabl, ndims )
     REAL (DP), SAVE :: bw_tablex_local( ltabl, ndims ), bw_tabley_local( ltabl, ndims )
     REAL (DP), SAVE :: fw_tablex_exx( ltabl, ndims ), fw_tabley_exx( ltabl, ndims )
     REAL (DP), SAVE :: bw_tablex_exx( ltabl, ndims ), bw_tabley_exx( ltabl, ndims )

     IF(PRESENT(is_exx))THEN
        is_exx_ = is_exx
     ELSE
        is_exx_ = .FALSE.
     END IF

     dofft( 1 : nx ) = .TRUE.
     IF( PRESENT( pl2ix ) ) THEN
       IF( SIZE( pl2ix ) < nx ) &
         CALL fftx_error__( ' cft_2xy ', ' wrong dimension for arg no. 8 ', 1 )
       DO i = 1, nx
         IF( pl2ix(i) < 1 ) dofft( i ) = .FALSE.
       END DO
     END IF

     ! WRITE( stdout,*) 'DEBUG: ', COUNT( dofft )

     !
     !   Here initialize table only if necessary
     !

     DO ip = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       IF ( is_exx_ ) THEN
          done = ( ny == dims_exx(1,ip) ) .AND. ( nx == dims_exx(3,ip) )
          done = done .AND. ( ldx == dims_exx(2,ip) ) .AND.  ( nzl == dims_exx(4,ip) )
       ELSE
          done = ( ny == dims_local(1,ip) ) .AND. ( nx == dims_local(3,ip) )
          done = done .AND. ( ldx == dims_local(2,ip) ) .AND.  ( nzl == dims_local(4,ip) )
       END IF
       IF (done) EXIT

     END DO

     IF( .NOT. done ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       ! WRITE( stdout, fmt="('DEBUG cft_2xy, reinitializing tables ', I3)" ) icurrent

#if defined(__OPENMP)

       tscale = 1.0_DP / ( nx * ny )
       IF ( is_exx_ ) THEN
          CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, nx,  1, 1.0_DP, &
               fw_tabley_exx( 1, icurrent_exx), ltabl, work(1), lwork )
          CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, nx, -1, 1.0_DP, &
               bw_tabley_exx(1, icurrent_exx), ltabl, work(1), lwork )
          CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny,  1, &
               tscale, fw_tablex_exx( 1, icurrent_exx), ltabl, work(1), lwork)
          CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, -1, &
               1.0_DP, bw_tablex_exx(1, icurrent_exx), ltabl, work(1), lwork)
       ELSE
          CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, nx,  1, 1.0_DP, &
               fw_tabley_local( 1, icurrent_local), ltabl, work(1), lwork )
          CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, nx, -1, 1.0_DP, &
               bw_tabley_local(1, icurrent_local), ltabl, work(1), lwork )
          CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny,  1, &
               tscale, fw_tablex_local( 1, icurrent_local), ltabl, work(1), lwork)
          CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, -1, &
               1.0_DP, bw_tablex_local(1, icurrent_local), ltabl, work(1), lwork)
       END IF

#else

       tscale = 1.0_DP / ( nx * ny )
       IF ( is_exx_ ) THEN
          CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, 1,  1, 1.0_DP, &
               fw_tabley_exx( 1, icurrent_exx), ltabl, work(1), lwork )
          CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, 1, -1, 1.0_DP, &
               bw_tabley_exx(1, icurrent_exx), ltabl, work(1), lwork )
          CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny,  1, &
               tscale, fw_tablex_exx( 1, icurrent_exx), ltabl, work(1), lwork)
          CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, -1, &
               1.0_DP, bw_tablex_exx(1, icurrent_exx), ltabl, work(1), lwork)
       ELSE
          CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, 1,  1, 1.0_DP, &
               fw_tabley_local( 1, icurrent_local), ltabl, work(1), lwork )
          CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, 1, -1, 1.0_DP, &
               bw_tabley_local(1, icurrent_local), ltabl, work(1), lwork )
          CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny,  1, &
               tscale, fw_tablex_local( 1, icurrent_local), ltabl, work(1), lwork)
          CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, -1, &
               1.0_DP, bw_tablex_local(1, icurrent_local), ltabl, work(1), lwork)
       END IF

#endif

       IF ( is_exx_ ) THEN
          dims_exx(1,icurrent_exx) = ny; dims_exx(2,icurrent_exx) = ldx;
          dims_exx(3,icurrent_exx) = nx; dims_exx(4,icurrent_exx) = nzl;
          ip = icurrent_exx
          icurrent_exx = MOD( icurrent_exx, ndims ) + 1
       ELSE
          dims_local(1,icurrent_local) = ny; dims_local(2,icurrent_local) = ldx;
          dims_local(3,icurrent_local) = nx; dims_local(4,icurrent_local) = nzl;
          ip = icurrent_local
          icurrent_local = MOD( icurrent_local, ndims ) + 1
       END IF

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'cft_2xy' )
#endif


#if defined(__OPENMP)

   IF( isign < 0 ) THEN
      tscale = 1.0_DP / ( nx * ny )
      do k = 1, nzl
         kk = 1 + ( k - 1 ) * ldx * ldy
         IF ( is_exx_ ) THEN
            CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, &
                 1, tscale, fw_tablex_exx( 1, ip ), ltabl, work( 1 ), lwork)
            CALL DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, nx, &
                 1, 1.0_DP, fw_tabley_exx(1, ip), ltabl, work( 1 ), lwork)
         ELSE
            CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, &
                 1, tscale, fw_tablex_local( 1, ip ), ltabl, work( 1 ), lwork)
            CALL DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, nx, &
                 1, 1.0_DP, fw_tabley_local(1, ip), ltabl, work( 1 ), lwork)
         END IF
      end do
   ELSE IF( isign > 0 ) THEN
      DO k = 1, nzl
         kk = 1 + ( k - 1 ) * ldx * ldy
         IF ( is_exx_ ) THEN
            CALL DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, nx, &
                 -1, 1.0_DP, bw_tabley_local(1, ip), ltabl, work( 1 ), lwork)
            CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, &
                 -1, 1.0_DP, bw_tablex_local(1, ip), ltabl, work( 1 ), lwork)
         ELSE
            CALL DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, nx, &
                 -1, 1.0_DP, bw_tabley_local(1, ip), ltabl, work( 1 ), lwork)
            CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, &
                 -1, 1.0_DP, bw_tablex_local(1, ip), ltabl, work( 1 ), lwork)
         END IF
      END DO
   END IF

#else

   IF( isign < 0 ) THEN
      idir = 1
      tscale = 1.0_DP / ( nx * ny )
      do k = 1, nzl
         kk = 1 + ( k - 1 ) * ldx * ldy
         IF ( is_exx_ ) THEN
            CALL DCFT ( 0, r(kk), 1, ldx, r(kk), 1, ldx, nx, ny, idir, &
                 tscale, fw_tablex_exx( 1, ip ), ltabl, work( 1 ), lwork)
         ELSE
            CALL DCFT ( 0, r(kk), 1, ldx, r(kk), 1, ldx, nx, ny, idir, &
                 tscale, fw_tablex_local( 1, ip ), ltabl, work( 1 ), lwork)
         END IF
         do i = 1, nx
            IF( dofft( i ) ) THEN
               kk = i + ( k - 1 ) * ldx * ldy
               IF ( is_exx_ ) THEN
                  call DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, 1, &
                       idir, 1.0_DP, fw_tabley_exx(1, ip), ltabl, work( 1 ), lwork)
               ELSE
                  call DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, 1, &
                       idir, 1.0_DP, fw_tabley_local(1, ip), ltabl, work( 1 ), lwork)
               END IF
            END IF
         end do
      end do
   ELSE IF( isign > 0 ) THEN
      idir = -1
      DO k = 1, nzl
         do i = 1, nx
            IF( dofft( i ) ) THEN
               kk = i + ( k - 1 ) * ldx * ldy
               IF ( is_exx_ ) THEN
                  call DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, 1, &
                       idir, 1.0_DP, bw_tabley_exx(1, ip), ltabl, work( 1 ), lwork)
               ELSE
                  call DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, 1, &
                       idir, 1.0_DP, bw_tabley_local(1, ip), ltabl, work( 1 ), lwork)
               END IF
            END IF
         end do
         kk = 1 + ( k - 1 ) * ldx * ldy
         IF ( is_exx_ ) THEN
            CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, idir, &
                 1.0_DP, bw_tablex_exx(1, ip), ltabl, work( 1 ), lwork)
         ELSE
            CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, idir, &
                 1.0_DP, bw_tablex_local(1, ip), ltabl, work( 1 ), lwork)
         END IF
      END DO
   END IF
#endif

#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'cft_2xy' )
#endif

     RETURN

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

   SUBROUTINE cfft3d( f, nx, ny, nz, ldx, ldy, ldz, isign, is_exx )

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
     COMPLEX (DP) :: f(:)
     LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
     LOGICAL :: is_exx_
     INTEGER :: i, k, j, err, idir, ip
     REAL(DP) :: tscale
     !INTEGER, SAVE :: icurrent = 1
     !INTEGER, SAVE :: dims(3,ndims) = -1
     INTEGER, SAVE :: icurrent_local = 1
     INTEGER, SAVE :: dims_local(3,ndims) = -1
     INTEGER, SAVE :: icurrent_exx = 1
     INTEGER, SAVE :: dims_local(3,ndims) = -1

     IF(PRESENT(is_exx))THEN
        is_exx_ = is_exx
     ELSE
        is_exx_ = .FALSE.
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
     ip = -1
     DO i = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       IF ( is_exx_ ) THEN
          IF ( ( nx == dims_exx(1,i) ) .and. &
               ( ny == dims_exx(2,i) ) .and. &
               ( nz == dims_exx(3,i) ) ) THEN
             ip = i
             EXIT
          END IF
       ELSE
          IF ( ( nx == dims_local(1,i) ) .and. &
               ( ny == dims_local(2,i) ) .and. &
               ( nz == dims_local(3,i) ) ) THEN
             ip = i
             EXIT
          END IF
       END IF
     END DO

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       ! no initialization for 3d FFT's from ESSL

       IF ( is_exx_ ) THEN
          dims_exx(1,icurrent_exx) = nx
          dims_exx(2,icurrent_exx) = ny
          dims_exx(3,icurrent_exx) = nz
          ip = icurrent_exx
          icurrent_exx = MOD( icurrent_exx, ndims ) + 1
       ELSE
          dims_local(1,icurrent_local) = nx
          dims_local(2,icurrent_local) = ny
          dims_local(3,icurrent_local) = nz
          ip = icurrent_local
          icurrent_local = MOD( icurrent_local, ndims ) + 1
       END IF

     END IF

     !
     !   Now perform the 3D FFT using the machine specific driver
     !

     IF ( isign < 0 ) THEN
       tscale = 1.0_DP / ( nx * ny * nz )
       idir = +1
     ELSE IF( isign > 0 ) THEN
       tscale = 1.0_DP
       idir = -1
     END IF

     IF( isign /= 0 ) CALL dcft3( f(1), ldx,ldx*ldy, f(1), ldx,ldx*ldy, &
          nx,ny,nz, idir, tscale, work(1), lwork)

     RETURN
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
  LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
  LOGICAL :: is_exx_
  !
  !   logical dimensions of the fft
  !   physical dimensions of the f array
  !   sign of the transformation

  complex(DP) :: f ( ldx * ldy * ldz )
  integer :: do_fft_x(:), do_fft_y(:)
  !
  integer :: m, incx1, incx2
  INTEGER :: i, k, j, err, idir, ip,  ii, jj
  REAL(DP) :: tscale
  !INTEGER, SAVE :: icurrent = 1
  !INTEGER, SAVE :: dims(3,ndims) = -1
  INTEGER, SAVE :: icurrent_local = 1
  INTEGER, SAVE :: dims_local(3,ndims) = -1
  INTEGER, SAVE :: icurrent_exx = 1
  INTEGER, SAVE :: dims_exx(3,ndims) = -1

  INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
  !REAL (DP), SAVE :: fw_table( ltabl, 3, ndims )
  !REAL (DP), SAVE :: bw_table( ltabl, 3, ndims )
  REAL (DP), SAVE :: fw_table_local( ltabl, 3, ndims )
  REAL (DP), SAVE :: bw_table_local( ltabl, 3, ndims )
  REAL (DP), SAVE :: fw_table_exx( ltabl, 3, ndims )
  REAL (DP), SAVE :: bw_table_exx( ltabl, 3, ndims )

  IF(PRESENT(is_exx))THEN
     is_exx_ = is_exx
  ELSE
     is_exx_ = .FALSE.
  END IF

  tscale = 1.0_DP

  ! WRITE( stdout, fmt="('DEBUG cfft3ds :',6I6)") nx, ny, nz, ldx, ldy, ldz
  ! WRITE( stdout, fmt="('DEBUG cfft3ds :',24I2)") do_fft_x
  ! WRITE( stdout, fmt="('DEBUG cfft3ds :',24I2)") do_fft_y


  IF( ny /= ldy ) &
    CALL fftx_error__(' cfft3ds ', ' wrong dimensions: ny /= ldy ', 1 )

     ip = -1
     DO i = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       IF ( is_exx_ ) THEN
          IF( ( nx == dims_exx(1,i) ) .and. ( ny == dims_exx(2,i) ) .and. &
               ( nz == dims_exx(3,i) ) ) THEN
             ip = i
             EXIT
          END IF
       ELSE
          IF( ( nx == dims_local(1,i) ) .and. ( ny == dims_local(2,i) ) .and. &
               ( nz == dims_local(3,i) ) ) THEN
             ip = i
             EXIT
          END IF
       END IF

     END DO

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       !
       ! ESSL sign convention for fft's is the opposite of the "usual" one
       !
       tscale = 1.0_DP
       IF ( is_exx_ ) THEN
          !  x - direction
          incx1 = 1; incx2 = ldx; m = 1
          CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nx, m,  1, 1.0_DP, &
               fw_table_exx( 1, 1, icurrent_exx), ltabl, work(1), lwork )
          CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nx, m, -1, 1.0_DP, &
               bw_table_exx(1, 1, icurrent_exx), ltabl, work(1), lwork )
          !  y - direction
          incx1 = ldx; incx2 = 1; m = nx;
          CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, ny, m,  1, 1.0_DP, &
               fw_table_exx( 1, 2, icurrent_exx), ltabl, work(1), lwork )
          CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, ny, m, -1, 1.0_DP, &
               bw_table_exx(1, 2, icurrent_exx), ltabl, work(1), lwork )
          !  z - direction
          incx1 = ldx * ldy; incx2 = 1; m = ldx * ny
          CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nz, m,  1, 1.0_DP, &
               fw_table_exx(1, 3, icurrent_exx), ltabl, work(1), lwork )
          CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nz, m, -1, 1.0_DP, &
               bw_table_exx(1, 3, icurrent_exx), ltabl, work(1), lwork )
       ELSE
          !  x - direction
          incx1 = 1; incx2 = ldx; m = 1
          CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nx, m,  1, 1.0_DP, &
               fw_table_local( 1, 1, icurrent_local), ltabl, work(1), lwork )
          CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nx, m, -1, 1.0_DP, &
               bw_table_local(1, 1, icurrent_local), ltabl, work(1), lwork )
          !  y - direction
          incx1 = ldx; incx2 = 1; m = nx;
          CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, ny, m,  1, 1.0_DP, &
               fw_table_local( 1, 2, icurrent_local), ltabl, work(1), lwork )
          CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, ny, m, -1, 1.0_DP, &
               bw_table_local(1, 2, icurrent_local), ltabl, work(1), lwork )
          !  z - direction
          incx1 = ldx * ldy; incx2 = 1; m = ldx * ny
          CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nz, m,  1, 1.0_DP, &
               fw_table_local(1, 3, icurrent_local), ltabl, work(1), lwork )
          CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nz, m, -1, 1.0_DP, &
               bw_table_local(1, 3, icurrent_local), ltabl, work(1), lwork )
       END IF

       IF ( is_exx_ ) THEN
          dims_exx(1,icurrent_exx) = nx
          dims_exx(2,icurrent_exx) = ny
          dims_exx(3,icurrent_exx) = nz
          ip = icurrent_exx
          icurrent_exx = MOD( icurrent_exx, ndims ) + 1
       ELSE
          dims_local(1,icurrent_local) = nx
          dims_local(2,icurrent_local) = ny
          dims_local(3,icurrent_local) = nz
          ip = icurrent_local
          icurrent_local = MOD( icurrent_local, ndims ) + 1
       END IF

     END IF


     IF ( isign > 0 ) THEN

        !
        !  i - direction ...
        !

        incx1 = 1;  incx2 = ldx;  m = 1

        do k = 1, nz
           do j = 1, ny
              jj = j + ( k - 1 ) * ldy
              ii = 1 + ldx * ( jj - 1 )
              if ( do_fft_x( jj ) == 1 ) THEN
                call dcft (0, f (ii), incx1,incx2, f (ii), incx1,incx2, nx, m, &
                -isign, 1.0_DP, bw_table ( 1, 1,  ip ), ltabl, work( 1 ), lwork)
              endif
           enddo
        enddo

        !
        !  ... j-direction ...
        !

        incx1 = ldx;  incx2 = 1;  m = nx

        do k = 1, nz
           ii = 1 + ldx * ldy * ( k - 1 )
           if ( do_fft_y( k ) == 1 ) then
             call dcft (0, f (ii), incx1, incx2, f (ii), incx1, incx2, nx, m, &
               -isign, 1.0_DP, bw_table ( 1, 2,  ip ), ltabl, work( 1 ), lwork)
           endif
        enddo

        !
        !     ... k-direction
        !

        incx1 = ldx * ldy;  incx2 = 1;  m = ldx * ny

        call dcft (0, f( 1 ), incx1, incx2, f( 1 ), incx1, incx2, nz, m, &
          -isign, 1.0_DP, bw_table ( 1, 3, ip ), ltabl, work( 1 ), lwork)

     ELSE

        !
        !     ... k-direction
        !

        incx1 = ldx * ny;  incx2 = 1;  m = ldx * ny

         call dcft (0, f( 1 ), incx1, incx2, f( 1 ), incx1, incx2, nz, m, &
          -isign, 1.0_DP, fw_table ( 1, 3, ip ), ltabl, work( 1 ), lwork)

        !
        !     ... j-direction ...
        !

        incx1 = ldx;  incx2 = 1;  m = nx

        do k = 1, nz
           ii = 1 + ldx * ldy * ( k - 1 )
           if ( do_fft_y ( k ) == 1 ) then
             call dcft (0, f (ii), incx1, incx2, f (ii), incx1, incx2, ny, m, &
               -isign, 1.0_DP, fw_table ( 1, 2, ip ), ltabl, work( 1 ), lwork)
           endif
        enddo

        !
        !     i - direction ...
        !

        incx1 = 1;  incx2 = ldx;  m = 1

        do k = 1, nz
           do j = 1, ny
              jj = j + ( k - 1 ) * ldy
              ii = 1 + ldx * ( jj - 1 )
              if ( do_fft_x( jj ) == 1 ) then
                call dcft (0, f (ii), incx1,incx2, f (ii), incx1,incx2, nx, m, &
                 -isign, 1.0_DP, fw_table ( 1, 1, ip ), ltabl, work( 1 ), lwork)
              endif
           enddo
        enddo

        call DSCAL (2 * ldx * ldy * nz, 1.0_DP/(nx * ny * nz), f(1), 1)

     END IF
     RETURN
   END SUBROUTINE cfft3ds

!=----------------------------------------------------------------------=!
   END MODULE fft_scalar
!=----------------------------------------------------------------------=!

#endif
