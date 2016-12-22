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
!#if defined(__USE_3D_FFT) & defined(__USE_MANY_FFT)
!	PUBLIC :: fwfftm, invfftm
!#endif
  
  
  INTERFACE invfft
     !! invfft is the interface to both the standard fft **invfft_x**,
     !! and to the "box-grid" version **invfft_b**, used only in CP 
     !! (the latter has an additional argument)
     
     SUBROUTINE invfft_x( grid_type, f, dfft, dtgs, howmany, is_exx )
       USE fft_types,  ONLY: fft_type_descriptor
       USE task_groups,   ONLY: task_groups_descriptor
       USE fft_param,  ONLY :DP
       IMPLICIT NONE
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
       USE fft_param,  ONLY :DP
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ia
       TYPE(fft_box_descriptor), INTENT(IN) :: dfft
       COMPLEX(DP) :: f(:)
       LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
     END SUBROUTINE invfft_b
  END INTERFACE

!#if defined(__USE_MANY_FFT) & defined(__USE_3D_FFT)
!    INTERFACE invfftm
!        !many version
!        SUBROUTINE invfft_xm( grid_type, f, dfft, dtgs, howmany, is_exx )
!            USE fft_types,  ONLY: fft_type_descriptor
!            USE task_groups,   ONLY: task_groups_descriptor
!            IMPLICIT NONE
!            INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
!            CHARACTER(LEN=*),  INTENT(IN) :: grid_type
!            TYPE(fft_type_descriptor), INTENT(IN) :: dfft
!            TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
!            COMPLEX(DP) :: f(:,:)
!            LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
!            INTEGER, OPTIONAL, INTENT(IN) :: howmany
!        END SUBROUTINE invfft_xm
!    END INTERFACE
!#endif

  INTERFACE fwfft
     SUBROUTINE fwfft_x( grid_type, f, dfft, dtgs, howmany, is_exx )
       USE fft_types,  ONLY: fft_type_descriptor
       USE task_groups,   ONLY: task_groups_descriptor
       USE fft_param,  ONLY :DP
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN) :: grid_type
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft
       TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
       INTEGER, OPTIONAL, INTENT(IN) :: howmany
       COMPLEX(DP) :: f(:)
       LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
     END SUBROUTINE fwfft_x
  END INTERFACE
  
!#if defined(__USE_MANY_FFT) & defined(__USE_3D_FFT)
!	INTERFACE fwfftm
!       SUBROUTINE fwfft_xm( grid_type, f, dfft, dtgs, howmany, is_exx )
!         USE fft_types,  ONLY: fft_type_descriptor
!         USE task_groups,   ONLY: task_groups_descriptor
!         IMPLICIT NONE
!         INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
!         CHARACTER(LEN=*), INTENT(IN) :: grid_type
!         TYPE(fft_type_descriptor), INTENT(IN) :: dfft
!         TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
!         COMPLEX(DP) :: f(:,:)
!         LOGICAL, OPTIONAL, INTENT(IN) :: is_exx
!  	   	 INTEGER, OPTIONAL, INTENT(IN) :: howmany
!       END SUBROUTINE fwfft_xm
!	END INTERFACE
!#endif

END MODULE fft_interfaces
!=---------------------------------------------------------------------------=!
