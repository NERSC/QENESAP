!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------------=!
MODULE fft_interfaces

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: fwfft, invfft
#if defined(__USE_3D_FFT) & defined(__USE_MANY_FFT)
	PUBLIC :: fwfftm, invfftm
#endif
  
  
  INTERFACE invfft
     !! invfft is the interface to both the standard fft **invfft_x**,
     !! and to the "box-grid" version **invfft_b**, used only in CP 
     !! (the latter has an additional argument)
     
     SUBROUTINE invfft_x( grid_type, f, dfft, dtgs, howmany, is_exx )
       USE fft_types,  ONLY: fft_type_descriptor
       USE task_groups,   ONLY: task_groups_descriptor
       IMPLICIT NONE
       INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
       CHARACTER(LEN=*),  INTENT(IN) :: grid_type
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft
       TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
       INTEGER, OPTIONAL, INTENT(IN) :: howmany
       COMPLEX(DP) :: f(:)
       LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
     END SUBROUTINE invfft_x
     !
     SUBROUTINE invfft_b( f, dfft, ia, is_exx )
       USE fft_smallbox_type,  ONLY: fft_box_descriptor
       IMPLICIT NONE
       INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
       INTEGER, INTENT(IN) :: ia
       TYPE(fft_box_descriptor), INTENT(IN) :: dfft
       COMPLEX(DP) :: f(:)
       LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
     END SUBROUTINE invfft_b
  END INTERFACE

#if defined(__USE_MANY_FFT) & defined(__USE_3D_FFT)
    INTERFACE invfftm
        !many version
        SUBROUTINE invfft_xm( grid_type, f, dfft, dtgs, howmany, is_exx )
            USE fft_types,  ONLY: fft_type_descriptor
            USE task_groups,   ONLY: task_groups_descriptor
            IMPLICIT NONE
            INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
            CHARACTER(LEN=*),  INTENT(IN) :: grid_type
            TYPE(fft_type_descriptor), INTENT(IN) :: dfft
            TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
            COMPLEX(DP) :: f(:,:)
            LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
            INTEGER, OPTIONAL, INTENT(IN) :: howmany
        END SUBROUTINE invfft_xm
    END INTERFACE
#endif

  INTERFACE fwfft
     SUBROUTINE fwfft_x( grid_type, f, dfft, dtgs, howmany, is_exx )
       USE fft_types,  ONLY: fft_type_descriptor
       USE task_groups,   ONLY: task_groups_descriptor
       IMPLICIT NONE
       INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
       CHARACTER(LEN=*), INTENT(IN) :: grid_type
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft
       TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
       INTEGER, OPTIONAL, INTENT(IN) :: howmany
       COMPLEX(DP) :: f(:)
       LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
     END SUBROUTINE fwfft_x
  END INTERFACE
  
#if defined(__USE_MANY_FFT) & defined(__USE_3D_FFT)
	INTERFACE fwfftm
       SUBROUTINE fwfft_xm( grid_type, f, dfft, dtgs, howmany, is_exx )
         USE fft_types,  ONLY: fft_type_descriptor
         USE task_groups,   ONLY: task_groups_descriptor
         IMPLICIT NONE
         INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
         CHARACTER(LEN=*), INTENT(IN) :: grid_type
         TYPE(fft_type_descriptor), INTENT(IN) :: dfft
         TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
         COMPLEX(DP) :: f(:,:)
         LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
  	   	 INTEGER, OPTIONAL, INTENT(IN) :: howmany
       END SUBROUTINE fwfft_xm
	END INTERFACE
#endif

END MODULE fft_interfaces
!=---------------------------------------------------------------------------=!
!
!=---------------------------------------------------------------------------=!
SUBROUTINE invfft_x( grid_type, f, dfft, dtgs, howmany, is_exx )
  !! Compute G-space to R-space for a specific grid type
  !! 
  !! **grid_type = 'Dense'** : 
  !!   inverse fourier transform of potentials and charge density f
  !!   on the dense grid . On output, f is overwritten
  !! 
  !! **grid_type = 'Smooth'** :
  !!   inverse fourier transform of  potentials and charge density f
  !!   on the smooth grid . On output, f is overwritten
  !! 
  !! **grid_type = 'Wave'** :
  !!   inverse fourier transform of  wave functions f
  !!   on the smooth grid . On output, f is overwritten
  !! 
  !! **grid_type = 'Custom'** : 
  !!   inverse fourier transform of potentials and charge density f
  !!   on a custom grid. On output, f is overwritten
  !! 
  !! **grid_type = 'CustomWave'** :
  !!   inverse fourier transform of  wave functions f
  !!   on a custom grid. On output, f is overwritten
  !! 
  !! **dfft = FFT descriptor**, IMPORTANT NOTICE: grid is specified only by dfft.
  !!   No check is performed on the correspondence between dfft and grid_type.
  !!   grid_type is now used only to distinguish cases 'Wave' / 'CustomWave' 
  !!   from all other cases
  
  USE fft_scalar,    ONLY: cfft3d, cfft3ds
  USE fft_smallbox,  ONLY: cft_b, cft_b_omp
  USE fft_parallel,  ONLY: tg_cft3s
  USE fft_types,     ONLY: fft_type_descriptor
  USE task_groups,   ONLY: task_groups_descriptor

  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: grid_type
  COMPLEX(DP) :: f(:)
  LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
  LOGICAL :: is_exx_
  TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  INTEGER :: howmany_ = 1
  IF( present( is_exx ) ) THEN
     is_exx_ = is_exx
  ELSE
     is_exx_ = .FALSE.
  END IF

  IF(PRESENT(howmany) ) THEN
     howmany_ = howmany
  END IF
  !
  IF( grid_type == 'Dense' ) THEN
     CALL start_clock( 'fft' )
  ELSE IF( grid_type == 'Smooth' ) THEN
     CALL start_clock( 'ffts' )
  ELSE IF( grid_type == 'Wave' ) THEN
     CALL start_clock('fftw')
  ELSE IF( grid_type == 'Custom' ) THEN
     CALL start_clock('fftc')
  ELSE IF( grid_type == 'CustomWave' ) THEN
     CALL start_clock('fftcw')
  ELSE 
     CALL fftx_error__( ' invfft ', ' unknown grid: '//grid_type , 1 )
  END IF

  IF( dfft%lpara ) THEN

     IF( howmany_ /= 1 ) THEN
        CALL fftx_error__( ' invfft ', ' howmany not yet implemented for parallel driver ', 1 )
     END IF

! pcarrier@cray.com CHANGED: simplified syntax, and added is_exx syntax for USE_3D_FFT
#if defined(__MPI) && !defined(__USE_3D_FFT)
     
     IF( grid_type == 'Dense' .OR. grid_type == 'Smooth' .OR. grid_type == 'Custom' ) THEN
        CALL tg_cft3s( f, dfft, 1, is_exx=is_exx_ )
     ELSE IF( grid_type == 'Wave' .OR. grid_type == 'CustomWave' ) THEN
        CALL tg_cft3s( f, dfft, 2, dtgs, is_exx=Is_exx_ )
     END IF

#endif

#if defined(__MPI) && defined(__USE_3D_FFT)

     IF ( .NOT. is_exx_ ) THEN

        IF( grid_type == 'Dense' .OR. grid_type == 'Smooth' .OR. grid_type == 'Custom' ) THEN
           CALL tg_cft3s( f, dfft, 1, is_exx=is_exx_ )
        ELSE IF( grid_type == 'Wave' .OR. grid_type == 'CustomWave' ) THEN
           CALL tg_cft3s( f, dfft, 2, dtgs, is_exx=Is_exx_ )
        END IF

     ELSE
        CALL cfft3d( f, dfft%nr1, dfft%nr2, dfft%nr3, &
             dfft%nr1x,dfft%nr2x,dfft%nr3x, 1, 1)
     ENDIF
#endif

  ELSE

     IF( grid_type == 'Dense' .OR. grid_type == 'Smooth' .OR. &
         grid_type == 'Custom' ) THEN
        CALL cfft3d( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , 1)
     ELSE IF( grid_type == 'Wave' .OR. grid_type == 'CustomWave' ) THEN
        CALL cfft3ds( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                         dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , 1, &
                         dfft%isind, dfft%iplw )
     END IF

  END IF

  IF( grid_type == 'Dense' ) THEN
     CALL stop_clock( 'fft' )
  ELSE IF( grid_type == 'Smooth' ) THEN
     CALL stop_clock( 'ffts' )
  ELSE IF( grid_type == 'Wave' ) THEN
     CALL stop_clock('fftw')
  ELSE IF( grid_type == 'Custom' ) THEN
     CALL stop_clock('fftc')
  ELSE IF( grid_type == 'CustomWave' ) THEN
     CALL stop_clock('fftcw')
  END IF

  RETURN

END SUBROUTINE invfft_x
!=---------------------------------------------------------------------------=!
!
!=---------------------------------------------------------------------------=!
#if defined(__USE_3D_FFT) & defined(__USE_MANY_FFT)
SUBROUTINE invfft_xm( grid_type, f, dfft, dtgs, howmany, is_exx )
  !! Compute G-space to R-space for a specific grid type
  !! 
  !! **grid_type = 'Custom'** : 
  !!   inverse fourier transform of potentials and charge density f
  !!   on a custom grid. On output, f is overwritten. Only that is supported
  !! 
  !! **dfft = FFT descriptor**, IMPORTANT NOTICE: grid is specified only by dfft.
  !!   No check is performed on the correspondence between dfft and grid_type.
  !!   grid_type is now used only to distinguish cases 'Wave' / 'CustomWave' 
  !!   from all other cases
  
  USE fft_scalar,    ONLY: cfft3dm
  USE fft_smallbox,  ONLY: cft_b, cft_b_omp
  USE fft_parallel,  ONLY: tg_cft3s
  USE fft_types,     ONLY: fft_type_descriptor
  USE task_groups,   ONLY: task_groups_descriptor

  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: grid_type
  COMPLEX(DP) :: f(:)
  LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
  TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
  INTEGER, INTENT(IN) :: howmany
  LOGICAL :: is_exx_
  IF( present( is_exx ) ) THEN
     is_exx_ = is_exx
  ELSE
     is_exx_ = .FALSE.
  END IF
  !
  IF( grid_type == 'Custom' ) THEN
     CALL start_clock('fftcm')
  ELSE 
     CALL fftx_error__( ' invfft ', ' unknown grid: '//grid_type , 1 )
  END IF

     CALL cfft3dm( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                     dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany, 1)

  CALL stop_clock('fftcm')

  RETURN

END SUBROUTINE invfft_xm
#endif
!=---------------------------------------------------------------------------=!
!
!=---------------------------------------------------------------------------=!
SUBROUTINE fwfft_x( grid_type, f, dfft, dtgs, howmany, is_exx )
  !! Compute R-space to G-space for a specific grid type
  !! 
  !! **grid_type = 'Dense'**
  !!   forward fourier transform of potentials and charge density f
  !!   on the dense grid . On output, f is overwritten
  !! 
  !! **grid_type = 'Smooth'**
  !!   forward fourier transform of potentials and charge density f
  !!   on the smooth grid . On output, f is overwritten
  !! 
  !! **grid_type = 'Wave'**
  !!   forward fourier transform of  wave functions f
  !!   on the smooth grid . On output, f is overwritten
  !! 
  !! **grid_type = 'Custom'**
  !!   forward fourier transform of potentials and charge density f
  !!   on a custom grid . On output, f is overwritten
  !! 
  !! **grid_type = 'CustomWave'**
  !!   forward fourier transform of  wave functions
  !!   on a custom grid . On output, f is overwritten
  !! 
  !! **dfft = FFT descriptor**, IMPORTANT NOTICE: grid is specified only by dfft.
  !!   No check is performed on the correspondence between dfft and grid_type.
  !!   grid_type is now used only to distinguish cases 'Wave' / 'CustomWave' 
  !!   from all other cases
  
  USE fft_scalar,    ONLY: cfft3d, cfft3ds
  USE fft_parallel,  ONLY: tg_cft3s
  USE fft_types,     ONLY: fft_type_descriptor
  USE task_groups,   ONLY: task_groups_descriptor

  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: grid_type
  COMPLEX(DP) :: f(:)
  LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
  LOGICAL :: is_exx_
  TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  INTEGER :: howmany_ = 1
  IF( present(is_exx) ) THEN
     is_exx_ = is_exx
  ELSE
     is_exx_ = .FALSE.
  END IF

  IF(PRESENT(howmany) ) THEN
     howmany_ = howmany
  END IF

  IF( grid_type == 'Dense' ) THEN
     CALL start_clock( 'fft' )
  ELSE IF( grid_type == 'Smooth' ) THEN
     CALL start_clock( 'ffts' )
  ELSE IF( grid_type == 'Wave' ) THEN
     CALL start_clock( 'fftw' )
  ELSE IF( grid_type == 'Custom' ) THEN
     CALL start_clock('fftc')
  ELSE IF( grid_type == 'CustomWave' ) THEN
     CALL start_clock('fftcw')
  ELSE
     CALL fftx_error__( ' fwfft ', ' unknown grid: '//grid_type , 1 )
  END IF

  IF( dfft%lpara ) THEN

     IF( howmany_ /= 1 ) THEN
        CALL fftx_error__( ' fwfft ', ' howmany not yet implemented for parallel driver ', 1 )
     END IF

! pcarrier@cray.com CHANGED: simplified syntax, and added is_exx syntax for USE_3D_FFT
#if defined(__MPI) && !defined(__USE_3D_FFT)
     
     IF( grid_type == 'Dense' .OR. grid_type == 'Smooth' .OR. grid_type == 'Custom' ) THEN
        CALL tg_cft3s( f, dfft, -1, is_exx=is_exx_ )
     ELSE IF( grid_type == 'Wave' .OR. grid_type == 'CustomWave' ) THEN
        CALL tg_cft3s( f, dfft, -2, dtgs, is_exx=Is_exx_ )
  END IF

#endif

#if defined(__MPI) && defined(__USE_3D_FFT)

  IF ( .NOT. is_exx_ ) THEN

    IF( grid_type == 'Dense' .OR. grid_type == 'Smooth' .OR. grid_type == 'Custom' ) THEN
       CALL tg_cft3s( f, dfft, -1, is_exx=is_exx_ )
    ELSE IF( grid_type == 'Wave' .OR. grid_type == 'CustomWave' ) THEN
       CALL tg_cft3s( f, dfft, -2, dtgs, is_exx=Is_exx_ )
    END IF

 ELSE

     CALL cfft3d( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                     dfft%nr1x,dfft%nr2x,dfft%nr3x, 1, -1)
  ENDIF
#endif

  ELSE

     IF( grid_type == 'Dense' .OR. grid_type == 'Smooth' .OR. &
         grid_type == 'Custom' ) THEN
        CALL cfft3d( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , -1)
     ELSE IF( grid_type == 'Wave' .OR. grid_type == 'CustomWave' ) THEN
        CALL cfft3ds( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                         dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , -1, &
                         dfft%isind, dfft%iplw )
     END IF

  END IF

  IF( grid_type == 'Dense' ) THEN
     CALL stop_clock( 'fft' )
  ELSE IF( grid_type == 'Smooth' ) THEN
     CALL stop_clock( 'ffts' )
  ELSE IF( grid_type == 'Wave' ) THEN
     CALL stop_clock( 'fftw' )
  ELSE IF( grid_type == 'Custom' ) THEN
     CALL stop_clock('fftc')
  ELSE IF( grid_type == 'CustomWave' ) THEN
     CALL stop_clock('fftcw')
  END IF
  
  RETURN
  !
END SUBROUTINE fwfft_x
!=---------------------------------------------------------------------------=!
!
!=---------------------------------------------------------------------------=!
#if defined(__USE_3D_FFT) & defined(__USE_MANY_FFT)
SUBROUTINE fwfft_xm( grid_type, f, dfft, dtgs, howmany, is_exx )
  !! Compute R-space to G-space for a specific grid type
  !! 
  !! **grid_type = 'Custom'**
  !!   forward fourier transform of potentials and charge density f
  !!   on a custom grid . On output, f is overwritten
  !!
  !! **dfft = FFT descriptor**, IMPORTANT NOTICE: grid is specified only by dfft.
  !!   No check is performed on the correspondence between dfft and grid_type.
  !!   grid_type is now used only to distinguish cases 'Wave' / 'CustomWave' 
  !!   from all other cases
  
  USE fft_scalar,    ONLY: cfft3dm
  USE fft_parallel,  ONLY: tg_cft3s
  USE fft_types,     ONLY: fft_type_descriptor
  USE task_groups,   ONLY: task_groups_descriptor

  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: grid_type
  COMPLEX(DP) :: f(:)
  TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
  LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  LOGICAL :: is_exx_
  IF( present(is_exx) ) THEN
     is_exx_ = is_exx
  ELSE
     is_exx_ = .FALSE.
  END IF

  IF( grid_type == 'Custom' ) THEN
     CALL start_clock('fftc')
  ELSE
     CALL fftx_error__( ' fwfft ', ' unknown grid: '//grid_type , 1 )
  END IF


     CALL cfft3dm( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                     dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany, -1)

  CALL stop_clock('fftc')

  
  RETURN
  !
END SUBROUTINE fwfft_xm
#endif
!=---------------------------------------------------------------------------=!
!
!=---------------------------------------------------------------------------=!
SUBROUTINE invfft_b( f, dfft, ia, is_exx )
  !! Not-so-parallel 3d fft for box grid, implemented ONLY for sign=1
  !!
  !! ComputeG-space to R-space, $$ output = \sum_G f(G)exp(+iG*R) $$
  !! The array f (overwritten on output) is NOT distributed:
  !! a copy is present on each processor.
  !! The fft along z  is done on the entire grid.
  !! The fft along xy is done ONLY on planes that have components on the
  !! dense grid for each processor. Note that the final array will no
  !! longer be the same on all processors.
  !! 
  !! **grid_type** = 'Box' (only allowed value!)
  !! 
  !! **dfft** = fft descriptor for the box grid
  !! 
  !! **ia**   = index of the atom with a box grid. Used to find the number
  !!         of planes on this processors, contained in dfft%np3(ia)
  
  USE fft_scalar,    ONLY: cfft3d, cfft3ds
  USE fft_smallbox,  ONLY: cft_b, cft_b_omp
  USE fft_parallel,  ONLY: tg_cft3s
  USE fft_smallbox_type, ONLY: fft_box_descriptor
  
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  
  TYPE(fft_box_descriptor), INTENT(IN) :: dfft
! Removed the 'OPTIONAL' attribute. When present, the specific interfaces
! 'invfft_x' and 'invfft_b' cannot be disambiguated when the generic interface
! call is made. This is a violation the Fortran standard. The Cray compiler
! errors out on this, while the Intel only issues a warning. --rbw
! INTEGER, OPTIONAL, INTENT(IN) :: ia
  INTEGER, INTENT(IN) :: ia
  COMPLEX(DP) :: f(:)
  !
  INTEGER :: imin3, imax3, np3
  LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
  LOGICAL :: is_exx_
  IF( present(is_exx) )THEN
     is_exx_ = is_exx
  ELSE
     is_exx_ = .FALSE.
  END IF

  ! clocks called inside a parallel region do not work properly!
  ! in the future we probably need a thread safe version of the clock

  CALL start_clock( 'fftb' )

#if defined(__MPI) && !defined(__USE_3D_FFT)
     
  IF( dfft%np3( ia ) > 0 ) THEN

#if defined(__OPENMP)

     CALL cft_b_omp( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, &
                        dfft%imin3( ia ), dfft%imax3( ia ), 1 )

#else
     CALL cft_b( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                    dfft%nr1x,dfft%nr2x,dfft%nr3x, &
                    dfft%imin3( ia ), dfft%imax3( ia ), 1 )
#endif

  END IF

#else

#if defined(__OPENMP)
  CALL cft_b_omp( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                     dfft%nr1x,dfft%nr2x,dfft%nr3x, &
                     dfft%imin3( ia ), dfft%imax3( ia ), 1 )
#else
  CALL cfft3d( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                  dfft%nr1x,dfft%nr2x,dfft%nr3x, 1, 1)
#endif

#endif

  CALL stop_clock( 'fftb' )

  RETURN
END SUBROUTINE invfft_b
!=---------------------------------------------------------------------------=!
