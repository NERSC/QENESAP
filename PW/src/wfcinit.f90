! 
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE wfcinit()
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes an estimate of the starting wavefunctions
  ! ... from superposition of atomic wavefunctions and/or random wavefunctions.
  ! ... It also open needed files or memory buffers
  !
  USE io_global,            ONLY : stdout
  USE basis,                ONLY : natomwfc, starting_wfc
  USE bp,                   ONLY : lelfield
  USE klist,                ONLY : xk, nks, ngk
  USE control_flags,        ONLY : io_level, lscf
  USE fixed_occ,            ONLY : one_atom_occupations
  USE ldaU,                 ONLY : lda_plus_u, U_projection, wfcU
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE io_files,             ONLY : nwordwfc, nwordwfcU, iunhub, iunwfc, iunigk
  USE buffers,              ONLY : open_buffer, get_buffer, save_buffer
  USE uspp,                 ONLY : nkb, vkb
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npw, current_k, igk
  USE wannier_new,          ONLY : use_wannier
  USE pw_restart,           ONLY : pw_readfile
  USE mp_bands,             ONLY : nbgrp, root_bgrp,inter_bgrp_comm
  USE mp,                   ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ierr
  LOGICAL :: exst_mem, exst_file
  !
  !
  CALL start_clock( 'wfcinit' )
  !
  ! ... Orthogonalized atomic functions needed for LDA+U and other cases
  !
  IF ( use_wannier .OR. one_atom_occupations ) CALL orthoatwfc ( use_wannier )
  IF ( lda_plus_u ) CALL orthoUwfc()
  !
  ! ... open files/buffer for wavefunctions (nwordwfc set in openfil)
  ! ... io_level > 1 : open file, otherwise: open buffer
  !
  CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst_file )
  !
  ! ... now the various possible wavefunction initializations
  ! ... first a check: is "tmp_dir"/"prefix".wfc found on disk?
  !
  IF ( TRIM(starting_wfc) == 'file' .AND. .NOT. exst_file) THEN
     !
     ! ... "tmp_dir"/"prefix".wfc not found on disk: try to read
     ! ... wavefunctions in "collected" format from "prefix".save/, 
     ! ... rewrite them (in pw_readfile) using the internal format
     !
     ierr = 1
     CALL pw_readfile( 'wave', ierr )
     IF ( ierr > 0 ) THEN
        WRITE( stdout, '(5X,"Cannot read wfc : file not found")' )
        starting_wfc = 'atomic+random'
     END IF
     !
     ! ... workaround: with k-point parallelization and 1 k-point per pool,
     ! ... pw_readfile does not leave evc properly initialized on all pools
     !
     IF ( nks == 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, 1 )
     !
  ELSE IF ( TRIM(starting_wfc) == 'file' .AND. exst_file) THEN
     !
     ! ... wavefunctions are read from file (or buffer) in routine 
     ! ... c_bands, but not if there is a single k-point. In such
     ! ... a case, we read wavefunctions (directly from file in 
     ! ... order to avoid a useless buffer allocation) here
     !
     IF ( nks == 1 ) CALL davcio ( evc, 2*nwordwfc, iunwfc, nks, -1 )
     !
  END IF
  !
  ! ... state what will happen
  !
  IF ( TRIM(starting_wfc) == 'file' ) THEN
     !
     WRITE( stdout, '(5X,"Starting wfc from file")' )
     !
  ELSE IF ( starting_wfc == 'atomic' ) THEN
     !
     IF ( natomwfc >= nbnd ) THEN
        WRITE( stdout, '(5X,"Starting wfc are ",I4," atomic wfcs")' ) natomwfc
     ELSE
        WRITE( stdout, '(5X,"Starting wfc are ",I4," atomic + ", &
             &           I4," random wfc")' ) natomwfc, nbnd-natomwfc
     END IF
     !
  ELSE IF ( TRIM(starting_wfc) == 'atomic+random' .AND. natomwfc > 0) THEN
     !
     IF ( natomwfc >= nbnd ) THEN
        WRITE( stdout, '(5X,"Starting wfc are ",I4," randomized atomic wfcs")')&
             natomwfc
     ELSE
        WRITE( stdout, '(5X,"Starting wfc are ",I4," randomized atomic wfcs + "&
             &          ,I4," random wfc")' ) natomwfc, nbnd-natomwfc
     END IF
     !
  ELSE
     !
     WRITE( stdout, '(5X,"Starting wfc are random")' )
     !
  END IF
  !
  ! ... exit here if starting from file or for non-scf calculations.
  ! ... In the latter case the starting wavefunctions are not 
  ! ... calculated here but just before diagonalization (to reduce I/O)
  !
  IF (  ( .NOT. lscf .AND. .NOT. lelfield ) .OR. TRIM(starting_wfc) == 'file' ) THEN
     !
     CALL stop_clock( 'wfcinit' )
     RETURN
     !
  END IF
  !
  IF ( nks > 1 ) REWIND( iunigk )
  !
  ! ... calculate and write all starting wavefunctions to file
  !
  DO ik = 1, nks
     !
     ! ... various initializations: k, spin, number of PW, indices
     !
     current_k = ik
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk (ik)
     IF ( nks > 1 ) READ( iunigk ) igk
     !
     call g2_kin (ik)
     !
     ! ... Calculate nonlocal pseudopotential projectors |beta>
     !
     IF ( nkb > 0 ) CALL init_us_2( npw, igk, xk(1,ik), vkb )
     !
     ! ... Needed for LDA+U
     !
     IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
        CALL get_buffer( wfcU, nwordwfcU, iunhub, ik )
     !
     ! ... calculate starting wavefunctions
     !
     CALL init_wfc ( ik )

!     ! BGRP:  make sure all bgrp have the same starting wfc
!     IF ( nbgrp > 1 ) CALL mp_bcast(evc,root_bgrp,inter_bgrp_comm)

     !
     ! ... write  starting wavefunctions to file
     !
     IF ( nks > 1 .OR. (io_level > 1) .OR. lelfield ) &
         CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
     !
  END DO
  !
  CALL stop_clock( 'wfcinit' )
  RETURN
  !
END SUBROUTINE wfcinit
!
!----------------------------------------------------------------------------
SUBROUTINE init_wfc ( ik )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes starting wavefunctions for k-point ik
  !
  USE kinds,                ONLY : DP
  USE bp,                   ONLY : lelfield
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   bec_type, becp
  USE constants,            ONLY : tpi
  USE cell_base,            ONLY : tpiba2
  USE basis,                ONLY : natomwfc, starting_wfc
  USE gvect,                ONLY : g, gstart
  USE klist,                ONLY : xk
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, et
  USE uspp,                 ONLY : nkb, okvan
  USE noncollin_module,     ONLY : npol
  USE wavefunctions_module, ONLY : evc
  USE random_numbers,       ONLY : randy
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm, my_bgrp_id
  USE mp,                   ONLY : mp_sum
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  !
  INTEGER :: ibnd, ig, ipol, n_starting_wfc, n_starting_atomic_wfc
  LOGICAL :: lelfield_save
  !
  REAL(DP) :: rr, arg
  REAL(DP), ALLOCATABLE :: etatom(:) ! atomic eigenvalues
  !
  COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:,:) ! atomic wfcs for initialization
  !
  !
  IF ( starting_wfc(1:6) == 'atomic' ) THEN
     !
     n_starting_wfc = MAX( natomwfc, nbnd )
     n_starting_atomic_wfc = natomwfc
     !
  ELSE IF ( starting_wfc == 'random' ) THEN
     !
     n_starting_wfc = nbnd
     n_starting_atomic_wfc = 0
     !
  ELSE
     !
     ! ...case 'file' should not be done here
     !
     CALL errore ( 'init_wfc', &
          'invalid value for startingwfc: ' // TRIM ( starting_wfc ) , 1 )
     !
  END IF
  !
  ALLOCATE( wfcatom( npwx, npol, n_starting_wfc ) )
  !
  IF ( starting_wfc(1:6) == 'atomic' ) THEN
     !
     CALL atomic_wfc( ik, wfcatom )
     !
     IF ( starting_wfc == 'atomic+random' .AND. &
         n_starting_wfc == n_starting_atomic_wfc ) THEN
         !
         ! ... in this case, introduce a small randomization of wavefunctions
         ! ... to prevent possible "loss of states"
         !
         DO ibnd = 1, n_starting_atomic_wfc
            !
            DO ipol = 1, npol
               !
               DO ig = 1, npw
                  !
                  rr  = randy()
                  arg = tpi * randy()
                  !
                  wfcatom(ig,ipol,ibnd) = wfcatom(ig,ipol,ibnd) * &
                     ( 1.0_DP + 0.05_DP * CMPLX( rr*COS(arg), rr*SIN(arg) ,kind=DP) ) 
                  !
               END DO
               !
            END DO
            !
         END DO
         !
     END IF
     !
  END IF
  !
  ! ... if not enough atomic wfc are available,
  ! ... fill missing wfcs with random numbers
  !
  DO ibnd = n_starting_atomic_wfc + 1, n_starting_wfc
     !
     DO ipol = 1, npol
        ! 
        wfcatom(:,ipol,ibnd) = (0.0_dp, 0.0_dp)
        !
        DO ig = 1, npw
           !
           rr  = randy()
           arg = tpi * randy()
           !
           wfcatom(ig,ipol,ibnd) = &
                CMPLX( rr*COS( arg ), rr*SIN( arg ) ,kind=DP) / &
                       ( ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                         ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                         ( xk(3,ik) + g(3,igk(ig)) )**2 + 1.0_DP )
        END DO
        !
     END DO
     !
  END DO
  
  ! when band paralleization is active, the first bgrp distributes the wfc to the others
  ! making sure all bgrp have the same starting wfc
  if (my_bgrp_id > 0) wfcatom(:,:,:) = (0.d0,0.d0) ; call mp_sum(wfcatom,inter_bgrp_comm)

  !
  ! ... Diagonalize the Hamiltonian on the basis of atomic wfcs
  !
  ALLOCATE( etatom( n_starting_wfc ) )
  !
  ! ... Allocate space for <beta|psi>
  !
  CALL allocate_bec_type ( nkb, n_starting_wfc, becp, intra_bgrp_comm )
  !
  ! ... the following trick is for electric fields with Berry's phase:
  ! ... by setting lelfield = .false. one prevents the calculation of
  ! ... electric enthalpy in the Hamiltonian (cannot be calculated
  ! ... at this stage: wavefunctions at previous step are missing)
  ! 
  lelfield_save = lelfield
  lelfield = .FALSE.
  !
  CALL rotate_wfc ( npwx, npw, n_starting_wfc, gstart, &
                    nbnd, wfcatom, npol, okvan, evc, etatom )
  !
  lelfield = lelfield_save
  !
  ! ... copy the first nbnd eigenvalues
  ! ... eigenvectors are already copied inside routine rotate_wfc
  !
  et(1:nbnd,ik) = etatom(1:nbnd)
  !
  CALL deallocate_bec_type ( becp )
  DEALLOCATE( etatom )
  DEALLOCATE( wfcatom )
  !
  RETURN
  !
END SUBROUTINE init_wfc
