!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!=----------------------------------------------------------------------=!
   MODULE fft_scalar_sx6
!=----------------------------------------------------------------------=!

     USE fft_param
       
     IMPLICIT NONE
     SAVE
#if defined(__SX6)
        PRIVATE
        PUBLIC :: cft_1z, cft_2xy, cfft3d, cfft3ds

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

     INTEGER :: tid

     ! ...   Machine-Dependent parameters, work arrays and tables of factors

     !   ltabl   Dimension of the tables of factors calculated at the
     !           initialization stage

#if defined(__OPENMP)
     INTEGER :: offset, ldz_t
     INTEGER :: omp_get_max_threads
     EXTERNAL :: omp_get_max_threads
#endif

     !   NEC MathKeisan

     INTEGER, PARAMETER :: ltabl = 2 * nfftx + 64
     !REAL (DP), SAVE :: tablez (ltabl, ndims)
     REAL (DP), SAVE :: tablez_local (ltabl, ndims)
     REAL (DP), SAVE :: tablez_exx (ltabl, ndims)
     REAL (DP)       :: work(4*nz*nsl)
     COMPLEX (DP)    :: DUMMY
     !INTEGER, SAVE :: isys = 1
     INTEGER, SAVE :: isys_local = 1
     INTEGER, SAVE :: isys_exx = 1

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

        IF ( is_exx_ ) THEN
           done = ( nz == zdims_exx(1,ip) )
        ELSE
           done = ( nz == zdims_local(1,ip) )
        END IF
        IF (done) EXIT
     END DO

     IF( .NOT. done ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       ! WRITE( stdout, fmt="('DEBUG cft_1z, reinitializing tables ', I3)" ) icurrent

       IF ( is_exx_ ) THEN

          CALL ZZFFTM (0, nz, 1, 1.0_DP, DUMMY, ldz, DUMMY, ldz, &
               tablez_exx (1, icurrent_exx), work, isys_exx)

          zdims_exx(1,icurrent_exx) = nz
          zdims_exx(2,icurrent_exx) = nsl
          zdims_exx(3,icurrent_exx) = ldz
          ip = icurrent_exx
          icurrent_exx = MOD( icurrent_exx, ndims ) + 1

       ELSE
          
          CALL ZZFFTM (0, nz, 1, 1.0_DP, DUMMY, ldz, DUMMY, ldz, &
               tablez_local (1, icurrent_local), work, isys_local)

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


     IF ( isign < 0 ) THEN
        idir   = -1
        tscale = 1.0_DP / nz
     ELSE IF ( isign > 0 ) THEN
        idir   = 1
        tscale = 1.0_DP
     END IF
     IF ( is_exx_ ) THEN
        IF (isign /= 0) CALL ZZFFTM (idir, nz, nsl, tscale, c(1), ldz, &
             cout(1), ldz, tablez_exx (1, ip), work, isys_exx)
     ELSE
        IF (isign /= 0) CALL ZZFFTM (idir, nz, nsl, tscale, c(1), ldz, &
             cout(1), ldz, tablez_local (1, ip), work, isys_local)
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
     INTEGER, SAVE :: icurrent_local = 1
     INTEGER, SAVE :: dims_local( 4, ndims) = -1
     INTEGER, SAVE :: icurrent_exx = 1
     INTEGER, SAVE :: dims_exx( 4, ndims) = -1
     LOGICAL :: dofft( nfftx ), done
     INTEGER, PARAMETER  :: stdout = 6

#if defined(__OPENMP)
     INTEGER :: offset
     INTEGER :: nx_t, ny_t, nzl_t, ldx_t, ldy_t
     INTEGER  :: itid, mytid, ntids
     INTEGER  :: omp_get_thread_num, omp_get_num_threads,omp_get_max_threads
     EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


     INTEGER, PARAMETER :: ltabl = 2*nfftx + 64
     REAL (DP), SAVE :: tablex_local(ltabl, ndims), tabley_local(ltabl, ndims)
     REAL (DP), SAVE :: tablex_exx(ltabl, ndims), tabley_exx(ltabl, ndims)
     REAL (DP)       :: work(4*nx*ny)
     COMPLEX (DP) :: XY(ldx*ny)
     COMPLEX (DP) :: DUMMY
     INTEGER, SAVE :: isys_local = 1
     INTEGER, SAVE :: isys_exx = 1

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

       IF ( is_exx_ )THEN
          CALL ZZFFT(0, ny, 1.0_DP, DUMMY, DUMMY, &
               tabley_exx (1, icurrent_exx), work, isys_exx)
          CALL ZZFFTM  (0, nx, 1, 1.0_DP, DUMMY, ldx, DUMMY, ldx, &
               tablex_exx(1, icurrent_exx), work, isys_exx)
          
          dims_exx(1,icurrent_exx) = ny; dims_exx(2,icurrent_exx) = ldx;
          dims_exx(3,icurrent_exx) = nx; dims_exx(4,icurrent_exx) = nzl;
          ip = icurrent_exx
          icurrent_exx = MOD( icurrent_exx, ndims ) + 1
       ELSE
          CALL ZZFFT(0, ny, 1.0_DP, DUMMY, DUMMY, &
               tabley_local (1, icurrent_local), work, isys_local)
          CALL ZZFFTM  (0, nx, 1, 1.0_DP, DUMMY, ldx, DUMMY, ldx, &
               tablex_local(1, icurrent_local), work, isys_local)
          
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


      IF( isign < 0 ) THEN

       idir = -1
       tscale = 1.0_DP / (nx * ny)
       DO k = 0, nzl-1
          kk = k * ldx * ldy
! FORWARD: ny FFTs in the X direction
          IF ( is_exx_ ) THEN
             CALL ZZFFTM ( idir, nx, ny, tscale, r(kk+1), ldx, r(kk+1), ldx, &
               tablex_exx (1, ip), work(1), isys_exx )
          ELSE
             CALL ZZFFTM ( idir, nx, ny, tscale, r(kk+1), ldx, r(kk+1), ldx, &
               tablex_local (1, ip), work(1), isys_local )
          END IF
! FORWARD: nx FFTs in the Y direction
          DO i = 1, nx
             IF ( dofft(i) ) THEN
                DO j = 0, ny-1
                   XY(j+1) = r(i + (j) * ldx + kk)
                END DO
                IF ( is_exx_ ) THEN
                   CALL ZZFFT(idir, ny, 1.0_DP, XY, XY, tabley_exx (1, ip), &
                        work(1), isys_exx)
                ELSE
                   CALL ZZFFT(idir, ny, 1.0_DP, XY, XY, tabley_local (1, ip), &
                        work(1), isys_local)
                END IF
                DO j = 0, ny-1
                   r(i + (j) * ldx + kk) = XY(j+1)
                END DO
             END IF
          END DO
       END DO

     ELSE IF ( isign > 0 ) THEN

       idir = 1
       tscale = 1.0_DP
       DO k = 0, nzl-1
! BACKWARD: nx FFTs in the Y direction
          kk = (k) * ldx * ldy
          DO i = 1, nx
             IF ( dofft(i) ) THEN
                DO j = 0, ny-1
                   XY(j+1) = r(i + (j) * ldx + kk)
                END DO
                IF ( is_exx_ ) THEN
                   CALL ZZFFT(idir, ny, 1.0_DP, XY, XY, tabley_exx (1, ip), &
                        work(1), isys_exx)
                ELSE
                   CALL ZZFFT(idir, ny, 1.0_DP, XY, XY, tabley_local (1, ip), &
                        work(1), isys_local)
                END IF
                DO j = 0, ny-1
                   r(i + (j) * ldx + kk) = XY(j+1)
                END DO
             END IF
          END DO
! BACKWARD: ny FFTs in the X direction
          IF ( is_exx_ ) THEN
             CALL ZZFFTM ( idir, nx, ny, tscale, r(kk+1), ldx, r(kk+1), ldx, &
                  tablex_exx (1, ip), work(1), isys_exx )
          ELSE
             CALL ZZFFTM ( idir, nx, ny, tscale, r(kk+1), ldx, r(kk+1), ldx, &
                  tablex_local (1, ip), work(1), isys_local )
          END IF
       END DO

     END IF

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

   SUBROUTINE cfft3d( f, nx, ny, nz, ldx, ldy, ldz, howmany, isign, is_exx )

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

     INTEGER, INTENT(IN) :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
     LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
     LOGICAL :: is_exx_
     COMPLEX (DP) :: f(:)
     INTEGER :: i, k, j, err, idir, ip
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent_local = 1
     INTEGER, SAVE :: dims_local(3,ndims) = -1
     INTEGER, SAVE :: icurrent_exx = 1
     INTEGER, SAVE :: dims_exx(3,ndims) = -1

     INTEGER, PARAMETER :: ltabl = 60
     INTEGER, PARAMETER :: lwork = 195+6*nfftx
     INTEGER, SAVE  :: iw0_local(ltabl, ndims)
     INTEGER, SAVE  :: iw0_exx(ltabl, ndims)
     INTEGER :: k_off, kj_offset
     REAL (DP), SAVE :: auxp_local (lwork, ndims)
     REAL (DP), SAVE :: auxp_exx (lwork, ndims)
     ! not sure whether auxp is work space or not
     COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: cw2
     COMPLEX (DP) :: f_out(size(f))

#if defined(ASL) && defined(MICRO)
     INTEGER :: nbtasks
     COMMON/NEC_ASL_PARA/nbtasks
#endif

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
         call fftx_error__('cfft3d',' nz is less than 1 ', 1)
     IF ( howmany /= 1 ) &
         call fftx_error__('cfft3d',' homany different from 1, not yet implemented for SX6 ', 1)

#if defined(ASL)
       ALLOCATE (cw2(ldx*ldy*ldz))
       CALL zfc3cl (f(1), nx, ny, nz, ldx, ldy, ldz, err)
#else
       ALLOCATE (cw2(6*ldx*ldy*ldz))
#endif
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

#if defined(ASL)
#if defined(MICRO)
       IF ( is_exx_ ) THEN
          CALL hfc3fb (nx,ny,nz, f(1) , ldx, ldy, ldz, 0, &
               iw0_exx(1,icurrent_exx), auxp_exx(1,icurrent_exx), cw2(1), nbtasks, err)
       ELSE
          CALL hfc3fb (nx,ny,nz, f(1) , ldx, ldy, ldz, 0, &
               iw0_local(1,icurrent_local), auxp_local(1,icurrent_local), cw2(1), nbtasks, err)
       END IF
#else
       IF ( is_exx_ ) THEN
          CALL zfc3fb (nx,ny,nz, f(1), ldx, ldy, ldz, 0, &
               iw0_exx(1,icurrent_exx), auxp_exx(1,icurrent_exx), cw2(1), err)
       ELSE
          CALL zfc3fb (nx,ny,nz, f(1), ldx, ldy, ldz, 0, &
               iw0_local(1,icurrent_local), auxp_local(1,icurrent_local), cw2(1), err)
       END IF
#endif
#else
       ! for some reason the error variable is not set by this driver on NEC SX machines
       err = 0 
       IF ( is_exx_ ) THEN
          CALL ZZFFT3D (0, nx,ny,nz, 1.0_DP, f(1), ldx, ldy, &
               &      f(1), ldx, ldy, auxp_exx(1,icurrent_exx), cw2(1), err)
       ELSE
          CALL ZZFFT3D (0, nx,ny,nz, 1.0_DP, f(1), ldx, ldy, &
               &      f(1), ldx, ldy, auxp_local(1,icurrent_local), cw2(1), err)
       END IF
#endif

       IF (err /= 0) CALL fftx_error__('cfft3d','FFT init returned an error ', err)

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

#if defined(ASL)
#if defined(MICRO)
     IF ( is_exx_ ) THEN
        CALL hfc3bf (nx,ny,nz, f(1), ldx,ldy, ldz, &
             -isign, iw0_exx(1,ip), auxp_exx(1,ip), cw2(1), nbtasks, err)
     ELSE
        CALL hfc3bf (nx,ny,nz, f(1), ldx,ldy, ldz, &
             -isign, iw0_local(1,ip), auxp_local(1,ip), cw2(1), nbtasks, err)
     END IF
#else
     IF ( is_exx_ ) THEN
        CALL zfc3bf (nx,ny,nz, f(1), ldx,ldy, ldz, &
             -isign, iw0_exx(1,ip), auxp_exx(1,ip), cw2(1), err)
     ELSE
        CALL zfc3bf (nx,ny,nz, f(1), ldx,ldy, ldz, &
             -isign, iw0_local(1,ip), auxp_local(1,ip), cw2(1), err)
     END IF
#endif
     IF ( isign < 0) THEN
        tscale = 1.0_DP / DBLE( nx * ny * nz )
        call ZDSCAL( ldx * ldy * ldz, tscale, f(1), 1)
     END IF
#else
     ! for some reason the error variable is not set by this driver on NEC SX machines
     err = 0 
     tscale = 1.0_DP
     IF ( isign < 0) THEN
        tscale = tscale / DBLE( nx * ny * nz )
     END IF
     IF ( is_exx_ ) THEN
        CALL ZZFFT3D (isign, nx,ny,nz, tscale, f(1), ldx,ldy, &
             f_out(1), ldx,ldy, auxp_exx(1,ip), cw2(1), err)
     ELSE
        CALL ZZFFT3D (isign, nx,ny,nz, tscale, f(1), ldx,ldy, &
             f_out(1), ldx,ldy, auxp_local(1,ip), cw2(1), err)
     END IF
!$omp parallel do private(j,i,k_off,kj_offset)
     do k=1,nz
        k_off = (k-1)*ldx*ldy
        do j=1,ny
           kj_offset = (j-1)*ldx + k_off
           do i=1,nx
              f(i+kj_offset) = f_out(i+kj_offset)
           end do
        end do
     end do
!$omp end parallel do
#endif

     IF (err /= 0) CALL fftx_error__('cfft3d','FFT returned an error ', err)
     DEALLOCATE(cw2)


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

SUBROUTINE cfft3ds (f, nx, ny, nz, ldx, ldy, ldz, howmany, isign, &
     do_fft_z, do_fft_y, is_exx)
  !
  !     driver routine for 3d complex "reduced" fft - see cfft3d
  !     The 3D fft are computed only on lines and planes which have
  !     non zero elements. These lines and planes are defined by
  !     the two integer vectors do_fft_y(nx) and do_fft_z(ldx*ny)
  !     (1 = perform fft, 0 = do not perform fft)
  !     This routine is implemented only for fftw, essl, acml
  !     If not implemented, cfft3d is called instead
  !
  !----------------------------------------------------------------------
  !
  implicit none

  integer :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
  LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
  LOGICAL :: is_exx_
  !
  !   logical dimensions of the fft
  !   physical dimensions of the f array
  !   sign of the transformation

  complex(DP) :: f ( ldx * ldy * ldz )
  integer :: do_fft_y(:), do_fft_z(:)
  !
  integer :: m, incx1, incx2
  INTEGER :: i, k, j, err, idir, ip,  ii, jj
  REAL(DP) :: tscale
  !INTEGER, SAVE :: icurrent = 1
  !INTEGER, SAVE :: dims(3,ndims) = -1

  IF(PRESENT(is_exx))THEN
     is_exx_ = is_exx
  ELSE
     is_exx_ = .FALSE.
  END IF

  CALL cfft3d (f, nx, ny, nz, ldx, ldy, ldz, howmany, isign, is_exx = is_exx_)
  RETURN
END SUBROUTINE cfft3ds
#endif
!=----------------------------------------------------------------------=!
END MODULE fft_scalar_sx6
!=----------------------------------------------------------------------=!
