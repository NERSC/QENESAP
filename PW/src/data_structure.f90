!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE data_structure( gamma_only, is_exx )
  !-----------------------------------------------------------------------
  ! this routine sets the data structure for the fft arrays
  ! (both the smooth and the dense grid)
  ! In the parallel case, it distributes columns to processes, too
  !
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_max
  USE mp_bands,   ONLY : me_bgrp, nproc_bgrp, root_bgrp, intra_bgrp_comm, &
                         ntask_groups
  USE mp_exx,   ONLY : me_egrp, nproc_egrp, root_egrp, intra_egrp_comm
  USE mp_pools,   ONLY : inter_pool_comm
  USE fft_base,   ONLY : dfftp, dffts
  USE cell_base,  ONLY : bg, tpiba
  USE klist,      ONLY : xk, nks
  USE gvect,      ONLY : gcutm, gvect_init, deallocate_gvect_exx
  USE gvecs,      ONLY : gcutms, gvecs_init, deallocate_gvecs
  USE stick_set,  ONLY : pstickset
  USE wvfct,      ONLY : ecutwfc
  USE io_global,  ONLY : stdout, ionode
  !
  IMPLICIT NONE
  LOGICAL, INTENT(in) :: gamma_only
  INTEGER, OPTIONAL, INTENT(in) :: is_exx
  INTEGER :: is_exx_
  REAL (DP) :: gkcut
  INTEGER :: ik, ngm_, ngs_, ngw_
  IF( PRESENT(is_exx) ) THEN
     is_exx_ = is_exx
  ELSE
     is_exx_ = 0
  END IF
  !
  ! ... calculate gkcut = max |k+G|^2, in (2pi/a)^2 units
  !
  IF (nks == 0) THEN
     !
     ! if k-points are automatically generated (which happens later)
     ! use max(bg)/2 as an estimate of the largest k-point
     !
     gkcut = 0.5d0 * max ( &
        sqrt (sum(bg (1:3, 1)**2) ), &
        sqrt (sum(bg (1:3, 2)**2) ), &
        sqrt (sum(bg (1:3, 3)**2) ) )
  ELSE
     gkcut = 0.0d0
     DO ik = 1, nks
        gkcut = max (gkcut, sqrt ( sum(xk (1:3, ik)**2) ) )
     ENDDO
  ENDIF
  gkcut = (sqrt (ecutwfc) / tpiba + gkcut)**2
  !
  ! ... find maximum value among all the processors
  !
  CALL mp_max (gkcut, inter_pool_comm )
  !
  ! ... set up fft descriptors, including parallel stuff: sticks, planes, etc.
  !
  IF (is_exx_.eq.1) THEN
     CALL pstickset( gamma_only, bg, gcutm, gkcut, gcutms, &
          dfftp, dffts, ngw_ , ngm_ , ngs_ , me_egrp, &
          root_egrp, nproc_egrp, intra_egrp_comm, ntask_groups, ionode, stdout )
  ELSE
     CALL pstickset( gamma_only, bg, gcutm, gkcut, gcutms, &
          dfftp, dffts, ngw_ , ngm_ , ngs_ , me_bgrp, &
          root_bgrp, nproc_bgrp, intra_bgrp_comm, ntask_groups, ionode, stdout )
  END IF
  !
  !     on output, ngm_ and ngs_ contain the local number of G-vectors
  !     for the two grids. Initialize local and global number of G-vectors
  !
  IF (is_exx_.gt.0) THEN
     call deallocate_gvect_exx()
     call deallocate_gvecs()
  END IF
  IF (is_exx_.eq.1) THEN
     call gvect_init ( ngm_ , intra_egrp_comm )
     call gvecs_init ( ngs_ , intra_egrp_comm );
  ELSE
     call gvect_init ( ngm_ , intra_bgrp_comm )
     call gvecs_init ( ngs_ , intra_bgrp_comm );
  END IF

END SUBROUTINE data_structure

