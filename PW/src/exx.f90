! Copyright (C) 2005-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------
MODULE exx
  !--------------------------------------
  !
  USE kinds,                ONLY : DP
  USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
                                   vcut_get,  vcut_spheric_get
  USE noncollin_module,     ONLY : noncolin, npol
  USE io_global,            ONLY : ionode
  USE fft_custom,           ONLY : fft_cus
  !
  USE control_flags, ONLY : tqr
  !<<<
  USE fft_types,            ONLY : fft_dlay_descriptor
  USE exx_old,              ONLY : exx_restart_old, exxinit_old, vexx_old, &
                                   exxenergy_old, exxenergy2_old, &
                                   exx_divergence_old, exx_stress_old
  !>>>

  IMPLICIT NONE
  SAVE
  COMPLEX(DP), ALLOCATABLE :: psi_exx(:,:), hpsi_exx(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc_exx(:,:), vkb_exx(:,:), psic_exx(:)
  INTEGER :: lda_original, n_original
  !
  ! general purpose vars
  !
  REAL(DP):: exxalfa=0._dp                ! 1 if exx, 0 elsewhere
  !
  ! variables defining the auxiliary k-point grid 
  ! used in X BZ integration
  !
  INTEGER :: nq1=1, nq2=1, nq3=1         ! integers defining the X integration mesh
  INTEGER :: nqs=1                       ! number of points in the q-grid
  INTEGER :: nkqs                        ! total number of different k+q
  !
  REAL(DP),    ALLOCATABLE :: xkq_collect(:,:)  ! xkq(3,nkqs) the auxiliary k+q set
  REAL(DP),    ALLOCATABLE :: x_occupation(:,:)           
                                         ! x_occupation(nbnd,nkstot) the weight of 
                                         ! auxiliary functions in the density matrix

  INTEGER :: x_nbnd_occ                  ! number of bands of auxiliary functions with
                                         ! at least some x_occupation > eps_occ

  INTEGER :: ibnd_start = 0              ! starting band index used in bgrp parallelization
  INTEGER :: ibnd_end = 0                ! ending band index used in bgrp parallelization

  COMPLEX(DP), ALLOCATABLE :: exxbuff(:,:,:)
                                         ! temporary buffer for wfc storage
  !
  !
  ! let xk(:,ik) + xq(:,iq) = xkq(:,ikq) = S(isym)*xk(ik') + G
  ! 
  !     index_xkq(ik,iq) = ikq
  !     index_xk(ikq)    = ik'
  !     index_sym(ikq)   = isym
  !
  INTEGER, ALLOCATABLE :: index_xkq(:,:) ! index_xkq(nks,nqs) 
  INTEGER, ALLOCATABLE :: index_xk(:)    ! index_xk(nkqs)  
  INTEGER, ALLOCATABLE :: index_sym(:)   ! index_sym(nkqs)
  INTEGER, ALLOCATABLE :: rir(:,:)       ! rotations to take k to q
!
!  Used for k points pool parallelization. All pools need these quantities.
!  They are allocated only IF needed.
!
  REAL(DP),    ALLOCATABLE :: xk_collect(:,:)
  REAL(DP),    ALLOCATABLE :: wk_collect(:)
  REAL(DP),    ALLOCATABLE :: wg_collect(:,:)
  !
  ! Internal:
  LOGICAL :: exx_grid_initialized = .false.
  !
  ! variables to deal with Coulomb divergence
  ! and related issues
  !
  REAL(DP)         :: eps  = 1.d-6
  REAL(DP)         :: eps_qdiv = 1.d-8 ! |q| > eps_qdiv
  REAL(DP), PARAMETER :: eps_occ  = 1.d-8 ! skip band with occupation < eps_occ
  REAL(DP)         :: exxdiv = 0._dp
  CHARACTER(32)    :: exxdiv_treatment  = ''
  !
  ! x_gamma_extrapolation
  LOGICAL           :: x_gamma_extrapolation =.TRUE.
  LOGICAl           :: on_double_grid =.FALSE.
  REAL(DP)          :: grid_factor = 1.d0 !8.d0/7.d0 
  !
  ! Gygi-Baldereschi 
  LOGICAL           :: use_regularization = .TRUE.
  !
  ! yukawa method
  REAL(DP)          :: yukawa = 0._dp
  !
  ! erfc screening
  REAL(DP)          :: erfc_scrlen = 0._dp
  !
! Debug Zhenfei Liu 9/25/2015
  REAL(DP)          :: beta_in_rsh = 0._dp
! DONE Debug.
  ! erf screening
  REAL(DP)          :: erf_scrlen = 0._dp
  !
  ! gau-pbe screening
  REAL (DP)         :: gau_scrlen = 0.d0
  !
  ! cutoff techniques
  LOGICAL           :: use_coulomb_vcut_ws = .FALSE.
  LOGICAL           :: use_coulomb_vcut_spheric = .FALSE.
  REAL(DP)          :: ecutvcut
  TYPE(vcut_type)   :: vcut

  !
  ! energy related variables
  !
  REAL(DP) :: fock0 = 0.0_DP, & !   sum <phi|Vx(phi)|phi>
              fock1 = 0.0_DP, & !   sum <psi|vx(phi)|psi>
              fock2 = 0.0_DP, & !   sum <psi|vx(psi)|psi>
              dexx  = 0.0_DP    !   fock1  - 0.5*(fock2+fock0)
  !
  ! custom fft grids
  !
  TYPE(fft_cus) exx_fft         ! Custom grid for Vx*psi calculation
  REAL(DP)  :: ecutfock         ! energy cutoff for custom grid
  !
  ! mapping for the data structure conversion
  !
  TYPE comm_packet
     INTEGER :: size
     INTEGER, ALLOCATABLE :: indices(:)
     COMPLEX(DP), ALLOCATABLE :: msg(:,:)
     COMPLEX(DP), ALLOCATABLE :: msg_evc(:,:)
     COMPLEX(DP), ALLOCATABLE :: msg_vkb(:,:)
  END TYPE comm_packet
  TYPE(comm_packet), ALLOCATABLE :: comm_recv(:,:), comm_send(:,:)
  TYPE(comm_packet), ALLOCATABLE :: comm_recv_reverse(:,:)
  TYPE(comm_packet), ALLOCATABLE :: comm_send_reverse(:,:)
  INTEGER, ALLOCATABLE :: lda_local(:,:)
  INTEGER, ALLOCATABLE :: lda_exx(:,:)
  INTEGER, ALLOCATABLE :: ngk_local(:), ngk_exx(:)
  INTEGER, ALLOCATABLE :: igk_exx(:)
  INTEGER :: npwx_local = 0
  INTEGER :: npwx_exx = 0
  INTEGER :: npw_local = 0
  INTEGER :: npw_exx = 0
  INTEGER :: nwordwfc_exx
  LOGICAL :: first_data_structure_change = .TRUE.

  INTEGER :: ngm_loc, ngm_g_loc, gstart_loc
  INTEGER, ALLOCATABLE :: ig_l2g_loc(:)
  REAL(DP), ALLOCATABLE :: g_loc(:,:), gg_loc(:)
  INTEGER, ALLOCATABLE :: mill_loc(:,:), nl_loc(:)
  INTEGER :: ngms_loc, ngms_g_loc
  INTEGER, ALLOCATABLE :: nls_loc(:)

  !gcutm
  INTEGER :: ngm_exx, ngm_g_exx, gstart_exx
  INTEGER, ALLOCATABLE :: ig_l2g_exx(:)
  REAL(DP), ALLOCATABLE :: g_exx(:,:), gg_exx(:)
  INTEGER, ALLOCATABLE :: mill_exx(:,:), nl_exx(:)
  INTEGER :: ngms_exx, ngms_g_exx
  INTEGER, ALLOCATABLE :: nls_exx(:)
  !gcutms

  !the coulomb factor is reused between iterations
  !REAL(DP), ALLOCATABLE :: coulomb_fac(:,:,:)
  !list of which coulomb factors have been calculated already
  !LOGICAL, ALLOCATABLE :: coulomb_done(:,:)

  TYPE(fft_dlay_descriptor) :: dfftp_loc, dffts_loc
  TYPE(fft_dlay_descriptor) :: dfftp_exx, dffts_exx
  INTEGER :: ngw_loc, ngs_loc
  INTEGER :: ngw_exx, ngs_exx

 CONTAINS
#define _CX(A)  CMPLX(A,0._dp,kind=DP)
#define _CY(A)  CMPLX(0._dp,-A,kind=DP)
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_fft_create ()
    USE wvfct,        ONLY : ecutwfc, npw
    USE gvect,        ONLY : ecutrho, ig_l2g
    USE uspp,         ONLY : okvan
    USE paw_variables,ONLY : okpaw
    USE control_flags,ONLY : gamma_only
    USE klist,        ONLY : qnorm
    USE cell_base,    ONLY : at, bg, tpiba2
    USE fft_custom,   ONLY : set_custom_grid, ggent
    USE grid_subroutines,   ONLY : realspace_grid_init

    IMPLICIT NONE

    IF( exx_fft%initialized) RETURN

    ! Initialise the custom grid that allows us to put the wavefunction
    ! onto the new (smaller) grid for rho (and vice versa)
    !
    exx_fft%ecutt=ecutwfc
    ! with k-points the following instructions guarantees that the sphere in 
    ! G space contains k+G points - needed if ecutfock \simeq ecutwfc
    IF ( gamma_only ) THEN
       exx_fft%dual_t = ecutfock/ecutwfc
    ELSE
       exx_fft%dual_t = MAX(ecutfock,(sqrt(ecutwfc)+qnorm)**2)/ecutwfc
    END IF
    !
    exx_fft%gcutmt = exx_fft%dual_t*exx_fft%ecutt / tpiba2
    CALL realspace_grid_init(exx_fft%dfftt, at, bg, exx_fft%gcutmt)
    CALL data_structure_custom(exx_fft, gamma_only)
    CALL ggent(exx_fft)
    exx_fft%initialized = .true.
    exx_fft%dfftt%have_task_groups = .FALSE.

    IF (gamma_only .AND. MAXVAL(ABS(ig_l2g(1:npw)-exx_fft%ig_l2gt(1:npw)))/=0) &
       CALL errore('exx_fft_create', ' exx fft grid not compatible with' &
                   //' the smooth fft grid ', 1 )

    RETURN
    !------------------------------------------------------------------------
  END SUBROUTINE exx_fft_create
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE deallocate_exx ()
    !------------------------------------------------------------------------
    !
    USE becmod, ONLY : deallocate_bec_type, is_allocated_bec_type, bec_type
    USE us_exx, ONLY : becxx
    USE fft_custom,  ONLY : deallocate_fft_custom
    !
    IMPLICIT NONE
    INTEGER :: ikq
    !
    IF ( ALLOCATED(index_xkq) ) DEALLOCATE(index_xkq)
    IF ( ALLOCATED(index_xk ) ) DEALLOCATE(index_xk )
    IF ( ALLOCATED(index_sym) ) DEALLOCATE(index_sym)
    IF ( ALLOCATED(rir)       ) DEALLOCATE(rir)
    IF ( ALLOCATED(x_occupation) ) DEALLOCATE(x_occupation)
    IF ( ALLOCATED(xkq_collect) )  DEALLOCATE(xkq_collect)
    IF ( ALLOCATED(exxbuff) )      DEALLOCATE(exxbuff)
    !
    IF(ALLOCATED(becxx)) THEN
      DO ikq = 1, nkqs
        IF(is_allocated_bec_type(becxx(ikq))) CALL deallocate_bec_type(becxx(ikq))
      ENDDO
      DEALLOCATE(becxx)
    ENDIF
    !
    CALL deallocate_fft_custom(exx_fft)
    !
    !  Pool variables deallocation
    !
    IF ( allocated (xk_collect) )  DEALLOCATE( xk_collect )
    IF ( allocated (wk_collect) )  DEALLOCATE( wk_collect )
    IF ( allocated (wg_collect) )  DEALLOCATE( wg_collect )
    !
    !
    !------------------------------------------------------------------------
  END SUBROUTINE deallocate_exx
  !------------------------------------------------------------------------
  !
  SUBROUTINE exx_grid_reinit()
    IMPLICIT NONE
    DEALLOCATE(xkq_collect,index_xk,index_sym)
    exx_grid_initialized = .false.
    CALL exx_grid_init()
  END SUBROUTINE exx_grid_reinit
  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_init()
    !------------------------------------------------------------------------
    !
    USE symm_base,  ONLY : nsym, s
    USE cell_base,  ONLY : bg, at
    USE spin_orb,   ONLY : domag
    USE noncollin_module, ONLY : nspin_lsda
    USE klist,      ONLY : xk, wk, nkstot, nks, qnorm
    USE wvfct,      ONLY : nbnd
    USE io_global,  ONLY : stdout
    USE start_k,    ONLY : nk1,nk2,nk3
    USE mp_pools,   ONLY : npool
    !
    IMPLICIT NONE
    !
    CHARACTER(13) :: sub_name='exx_grid_init'
    INTEGER       :: iq1, iq2, iq3, isym, ik, ikq, iq, max_nk, temp_nkqs
    INTEGER, allocatable :: temp_index_xk(:), temp_index_sym(:)
    INTEGER, allocatable :: temp_index_ikq(:), new_ikq(:)
    REAL(DP),allocatable :: temp_xkq(:,:)
    LOGICAL      :: xk_not_found
    REAL(DP)     :: sxk(3), dxk(3), xk_cryst(3)
    REAL(DP)     :: dq1, dq2, dq3
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    !
    CALL start_clock ('exx_grid')
    !
    IF(nq1<=0) nq1 = nk1
    IF(nq2<=0) nq2 = nk2
    IF(nq3<=0) nq3 = nk3
    IF(nkstot==nspin_lsda) THEN
      nq1=1; nq2=1; nq3=1
    ENDIF
     
    IF(ANY((/nq1,nq2,nq3/)<=0)) CALL errore('exx_grid_init',"wrong EXX q grid", 1)
    !
    IF(exx_grid_initialized) CALL errore('exx_grid_init', "grid already initialized",1)
    exx_grid_initialized = .true.
    !
    ! definitions and checks
    !
    grid_factor = 1._dp
    IF (x_gamma_extrapolation) &
        grid_factor = 8.d0/7.d0
    !
    nqs = nq1 * nq2 * nq3
    !
    ! all processors need to have access to all k+q points
    !
    IF ( .NOT.allocated (xk_collect) )  ALLOCATE(xk_collect(3,nkstot))
    IF ( .NOT.allocated (wk_collect) )  ALLOCATE(wk_collect(nkstot))
    ! the next if/then if probably not necessary, as xk_wk collect can
    ! deal with npool==1, leaving it for clarity.
    IF ( npool > 1 ) THEN
      CALL xk_wk_collect(xk_collect, wk_collect, xk, wk, nkstot, nks)
    ELSE
      xk_collect(:,1:nks) = xk(:,1:nks)
      wk_collect(1:nks) = wk(1:nks)
    ENDIF
    !
    ! set a safe limit as the maximum number of auxiliary points we may need
    ! and allocate auxiliary arrays
    max_nk = nkstot * min(48, 2 * nsym)
    ALLOCATE( temp_index_xk(max_nk), temp_index_sym(max_nk) )
    ALLOCATE( temp_index_ikq(max_nk), new_ikq(max_nk) )
    ALLOCATE( temp_xkq(3,max_nk) )
    !
    ! find all k-points equivalent by symmetry to the points in the k-list
    !
    temp_nkqs = 0
    DO isym=1,nsym
      DO ik =1, nkstot
          ! go to crystalline coordinates
          xk_cryst(:) = xk_collect(:,ik)
          CALL cryst_to_cart(1, xk_cryst, at, -1)
          ! rotate with this sym.op.
          sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                   s(:,2,isym)*xk_cryst(2) + &
                   s(:,3,isym)*xk_cryst(3)
          ! add sxk to the auxiliary list IF it is not already present
          xk_not_found = .true.
          ! *** do-loop skipped the first time because temp_nksq == 0
          DO ikq=1, temp_nkqs
            IF (xk_not_found ) THEN
                dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
                IF ( abs(dxk(1)).le.eps .and. &
                     abs(dxk(2)).le.eps .and. &
                     abs(dxk(3)).le.eps ) xk_not_found = .false.
            ENDIF
          ENDDO
          IF (xk_not_found) THEN
            temp_nkqs                 = temp_nkqs + 1
            temp_xkq(:,temp_nkqs)     = sxk(:)
            temp_index_xk(temp_nkqs)  = ik
            temp_index_sym(temp_nkqs) = isym 
          ENDIF

          sxk(:) = - sxk(:)
          xk_not_found = .true.
          DO ikq=1, temp_nkqs
            IF (xk_not_found ) THEN
                dxk(:) = sxk(:) - temp_xkq(:,ikq) - nint(sxk(:) - temp_xkq(:,ikq))
                IF ( abs(dxk(1)).le.eps .and. &
                     abs(dxk(2)).le.eps .and. &
                     abs(dxk(3)).le.eps ) xk_not_found = .false.
            ENDIF
          ENDDO
          IF (xk_not_found .and. .not. (noncolin.and.domag) ) THEN
            temp_nkqs                 = temp_nkqs + 1
            temp_xkq(:,temp_nkqs)     = sxk(:)
            temp_index_xk(temp_nkqs)  = ik
            temp_index_sym(temp_nkqs) =-isym 
          ENDIF

      ENDDO
    ENDDO

    !
    ! define the q-mesh step-sizes
    !
    dq1= 1._dp/DBLE(nq1)
    dq2= 1._dp/DBLE(nq2)
    dq3= 1._dp/DBLE(nq3)
    !
    ! allocate and fill the array index_xkq(nkstot,nqs)
    !
    if(.not.allocated(index_xkq))    ALLOCATE( index_xkq(nkstot,nqs) )
    if(.not.allocated(x_occupation)) ALLOCATE( x_occupation(nbnd,nkstot) )
    nkqs = 0
    new_ikq(:) = 0
    DO ik=1,nkstot
      ! go to crystalline coordinates
      xk_cryst(:) = xk_collect(:,ik)
      CALL cryst_to_cart(1, xk_cryst, at, -1)
      !
      iq = 0
      DO iq1=1, nq1
        sxk(1) = xk_cryst(1) + (iq1-1) * dq1
        DO iq2 =1, nq2
          sxk(2) = xk_cryst(2) + (iq2-1) * dq2
          DO iq3 =1, nq3
              sxk(3) = xk_cryst(3) + (iq3-1) * dq3
              iq = iq + 1

              xk_not_found = .true.
              !
              DO ikq=1, temp_nkqs
                IF ( xk_not_found ) THEN
                    dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
                    IF ( ALL(abs(dxk) < eps ) ) THEN
                        xk_not_found = .false.
                        IF ( new_ikq(ikq) == 0) THEN
                            nkqs = nkqs + 1
                            temp_index_ikq(nkqs) = ikq
                            new_ikq(ikq) = nkqs
                        ENDIF
                        index_xkq(ik,iq) = new_ikq(ikq)
                    ENDIF
                ENDIF
              ENDDO ! ikq
              !
              IF (xk_not_found) THEN
                write (*,*) ik, iq, temp_nkqs
                write (*,*) sxk(:)
                CALL errore(sub_name, ' k + q is not an S*k ', (ik-1) * nqs + iq )
              ENDIF

          ENDDO
        ENDDO
      ENDDO

    ENDDO
    !
    ! allocate and fill the arrays xkq(3,nkqs), index_xk(nkqs) and index_sym(nkqs)
    !
    ALLOCATE( xkq_collect(3,nspin_lsda*nkqs), index_xk(nspin_lsda*nkqs),  &
              index_sym(nspin_lsda*nkqs) )

    DO ik =1, nkqs
      ikq               = temp_index_ikq(ik)
      xkq_collect(:,ik) = temp_xkq(:,ikq)
      index_xk(ik)      = temp_index_xk(ikq)
      index_sym(ik)     = temp_index_sym(ikq)
    ENDDO
    CALL cryst_to_cart(nkqs, xkq_collect, bg, +1)

    IF( nkqs > 1) THEN
      WRITE(stdout, '(5x,3a)') "EXX: setup a grid of "//TRIM(int_to_char(nkqs))&
                           //" q-points centered on each k-point"
      WRITE( stdout, '(5x,a)' ) '(k+q)-points:'
      do ik = 1, nkqs
          WRITE( stdout, '(3f12.7,5x,2i5)') (xkq_collect (ikq, ik) , ikq = 1, 3) , &
                 index_xk(ik), index_sym(ik)
      enddo
    ELSE
      WRITE(stdout, '("EXX: grid of k+q points same as grid of k-points")')
    ENDIF
    
    ! if nspin == 2, the kpoints are repeated in couples (spin up, spin down)
    IF (nspin_lsda == 2) THEN
      DO ik = 1, nkstot/2
          DO iq =1, nqs
            index_xkq(nkstot/2+ik,iq) = index_xkq(ik,iq) + nkqs
          END DO
      ENDDO
      DO ikq=1,nkqs
          xkq_collect(:,ikq + nkqs) = xkq_collect(:,ikq)
          index_xk(ikq + nkqs)  = index_xk(ikq) + nkstot/2
          index_sym(ikq + nkqs) = index_sym(ikq)
      ENDDO
      nkqs = 2 * nkqs
    ENDIF
    !
    ! clean up
    DEALLOCATE(temp_index_xk, temp_index_sym, temp_index_ikq, new_ikq, temp_xkq)
    !
    ! check that everything is what it should be
    CALL exx_grid_check () 
    !
    ! qnorm = max |k+q|, useful for reduced-cutoff calculations with k-points 
    !
    qnorm = 0.0_dp
    DO iq = 1,nkqs
       DO ik = 1,nks
          qnorm = MAX(qnorm, SQRT( SUM((xk(:,ik)-xkq_collect(:,iq))**2) ))
       ENDDO
    ENDDO
    !
    CALL stop_clock ('exx_grid')
    !
    RETURN
    !------------------------------------------------------------------------
  END SUBROUTINE exx_grid_init
  !------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE exx_n_plane_waves(ecutwfc, tpiba2, g, ngm, npwx)
    !-----------------------------------------------------------------------
    !
    ! Find maximum number of plane waves npwx among the entire grid of k and
    ! of k+q points - should be called after a previous call to n_plane_waves 
    ! (for k-points only), providing in input the value of npwx found by
    ! n_plane_waves, to ensure that the final npwx is the largest of the two
    !
    USE kinds, ONLY : DP
    USE funct, ONLY : dft_is_hybrid
    USE uspp,  ONLY : okvan
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)  :: ngm
    REAL(DP),INTENT(in)  :: ecutwfc, tpiba2, g (3, ngm)
    INTEGER, INTENT(inout) :: npwx
    INTEGER, ALLOCATABLE :: ngkq(:)
    INTEGER :: npwx_
    !
    IF( .NOT. okvan .OR. .NOT.dft_is_hybrid() ) RETURN
    IF( .NOT.exx_grid_initialized) &
        CALL errore("exx_n_plane_waves","you must initialize the grid first",1)
    ALLOCATE(ngkq(nkqs))
    CALL n_plane_waves (ecutwfc, tpiba2, nkqs, xkq_collect, g, ngm, npwx_, ngkq)
    DEALLOCATE(ngkq)
    npwx = MAX (npwx, npwx_)
    !
    RETURN
    !------------------------------------------------------------------------
  END SUBROUTINE exx_n_plane_waves
  !------------------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_div_check()
    !------------------------------------------------------------------------
    !
    USE cell_base,  ONLY : at, alat
    USE io_global,  ONLY : stdout
    USE funct,      ONLY : get_screening_parameter, get_rsh_beta
! Debug Zhenfei Liu 9/26/2015, added get_rsh_beta above
    !
    IMPLICIT NONE
    !
    REAL(DP)     :: atws(3,3)
    CHARACTER(13) :: sub_name='exx_div_check'

    !
    ! EXX singularity treatment
    !
    SELECT CASE ( TRIM(exxdiv_treatment) ) 
    CASE ( "gygi-baldereschi", "gygi-bald", "g-b", "gb" )
      !
      use_regularization = .TRUE.
      !
      !
    CASE ( "vcut_ws" )
      !
      use_coulomb_vcut_ws = .TRUE.
      IF ( x_gamma_extrapolation ) &
            CALL errore(sub_name,'cannot USE x_gamm_extrap and vcut_ws', 1)
      !
    CASE ( "vcut_spherical" ) 
      !
      use_coulomb_vcut_spheric = .TRUE.
      IF ( x_gamma_extrapolation ) &
            CALL errore(sub_name,'cannot USE x_gamm_extrap and vcut_spherical', 1)
      !
    CASE ( "none" )
      use_regularization = .FALSE.
      !
    CASE DEFAULT
      CALL errore(sub_name,'invalid exxdiv_treatment: '//TRIM(exxdiv_treatment), 1)
    END SELECT
    !
    ! Set variables for Coulomb vcut
    ! NOTE: some memory is allocated inside this routine (in the var vcut)
    !       and should be deallocated somewehre, at the end of the run
    !
    IF ( use_coulomb_vcut_ws .OR. use_coulomb_vcut_spheric ) THEN
        !
        ! build the superperiodicity direct lattice
        !
        atws = alat * at
        !
        atws(:,1) = atws(:,1) * nq1
        atws(:,2) = atws(:,2) * nq2
        atws(:,3) = atws(:,3) * nq3
        !
        CALL vcut_init( vcut, atws, ecutvcut )
        !
        IF ( ionode ) CALL vcut_info( stdout, vcut )
        !          
    ENDIF
    RETURN
  !------------------------------------------------------------------------
  END SUBROUTINE exx_div_check 
  !------------------------------------------------------------------------


  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_check ( )
    !------------------------------------------------------------------------
    USE symm_base,  ONLY : s
    USE cell_base,  ONLY : at
    USE klist,      ONLY : nkstot, xk
    USE mp_pools,   ONLY : npool
    IMPLICIT NONE
    REAL(DP) :: sxk(3), dxk(3), xk_cryst(3), xkk_cryst(3)
    INTEGER :: iq1, iq2, iq3, isym, ik, ikk, ikq, iq
    REAL(DP) :: dq1, dq2, dq3
    dq1= 1._dp/DBLE(nq1)
    dq2= 1._dp/DBLE(nq2)
    dq3= 1._dp/DBLE(nq3)

    DO ik =1, nkstot
      xk_cryst(:) = xk_collect(:,ik)
      CALL cryst_to_cart(1, xk_cryst, at, -1)
      !
      iq = 0
      DO iq1=1, nq1
        sxk(1) = xk_cryst(1) + (iq1-1) * dq1
        DO iq2 =1, nq2
          sxk(2) = xk_cryst(2) + (iq2-1) * dq2
          DO iq3 =1, nq3
              sxk(3) = xk_cryst(3) + (iq3-1) * dq3
              iq = iq + 1
              
              ikq  = index_xkq(ik,iq) 
              ikk  = index_xk(ikq)
              isym = index_sym(ikq)

              IF (npool>1) THEN
                xkk_cryst(:) = at(1,:)*xk_collect(1,ikk)+at(2,:)*xk_collect(2,ikk)+at(3,:)*xk_collect(3,ikk)
              ELSE
                xkk_cryst(:) = at(1,:)*xk(1,ikk)+at(2,:)*xk(2,ikk)+at(3,:)*xk(3,ikk)
              ENDIF

              IF (isym < 0 ) xkk_cryst(:) = - xkk_cryst(:)
              isym = abs (isym)
              dxk(:) = s(:,1,isym)*xkk_cryst(1) + &
                       s(:,2,isym)*xkk_cryst(2) + &
                       s(:,3,isym)*xkk_cryst(3) - sxk(:)
              dxk(:) = dxk(:) - nint(dxk(:))
              IF ( .not. ( abs(dxk(1)).le.eps .and. &
                           abs(dxk(2)).le.eps .and. &
                           abs(dxk(3)).le.eps )   ) THEN
                  write(*,*) ik,iq
                  write(*,*) ikq,ikk,isym
                  write(*,*) dxk(:)
                  CALL errore('exx_grid_check', 'something wrong', 1 )
              ENDIF

          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    return

    !------------------------------------------------------------------------
  END SUBROUTINE exx_grid_check
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_restart(l_exx_was_active)
     !------------------------------------------------------------------------
     !This SUBROUTINE is called when restarting an exx calculation
     USE funct,                ONLY : get_exx_fraction, start_exx, &
! Debug Zhenfei Liu 9/25/2015
!                                      exx_is_active, get_screening_parameter
                                      exx_is_active, get_screening_parameter, &
                                      get_rsh_beta
! DONE Debug.
     USE mp_exx,               ONLY : use_old_exx

     IMPLICIT NONE
     LOGICAL, INTENT(IN) :: l_exx_was_active

     !<<<
     IF (use_old_exx) THEN
        CALL exx_restart_old(l_exx_was_active)
        RETURN
     END IF
     !>>>
     IF (.not. l_exx_was_active ) return ! nothing had happened yet
     !
     erfc_scrlen = get_screening_parameter()
! Debug Zhenfei Liu 9/25/2015
     beta_in_rsh = get_rsh_beta()
! DONE Debug.
     exxalfa = get_exx_fraction()
     exxdiv = exx_divergence() 
! Debug Zhenfei Liu 10/8/2015 - I think it's important to get exxalfa
! first, then execute the divergence. This is important as explicitly tested.

     CALL start_exx()
     CALL weights()
     CALL exxinit()
     fock0 = exxenergy2()
     RETURN
     !------------------------------------------------------------------------
  END SUBROUTINE exx_restart
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE exxinit()
  !------------------------------------------------------------------------

    ! This SUBROUTINE is run before the first H_psi() of each iteration.
    ! It saves the wavefunctions for the right density matrix, in real space
    !
    USE wavefunctions_module, ONLY : evc, psic
    USE io_files,             ONLY : nwordwfc, iunwfc_exx, iunigk_exx
    USE buffers,              ONLY : get_buffer
    USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg
    USE control_flags,        ONLY : gamma_only
    USE klist,                ONLY : ngk, nks, nkstot
    USE symm_base,            ONLY : nsym, s, sr, ftau
    USE mp_pools,             ONLY : npool, nproc_pool, me_pool, inter_pool_comm
    USE mp_exx,               ONLY : me_egrp, set_egrp_indices, negrp, &
                                     init_index_over_band, &
                                     inter_egrp_comm, intra_egrp_comm, &
                                     iexx_start, iexx_end, use_old_exx
    USE mp,                   ONLY : mp_sum
    USE funct,                ONLY : get_exx_fraction, start_exx,exx_is_active,&
                                     get_screening_parameter, &
                                     get_gau_parameter, get_rsh_beta
! Debug Zhenfei Liu 9/25/2015, added get_rsh_beta above
    USE scatter_mod,          ONLY : gather_grid, scatter_grid
    USE fft_interfaces,       ONLY : invfft
    USE becmod,               ONLY : allocate_bec_type, bec_type
    USE uspp,                 ONLY : nkb, okvan
    USE us_exx,               ONLY : becxx
    USE paw_variables,        ONLY : okpaw
    USE paw_exx,              ONLY : PAW_init_keeq

    IMPLICIT NONE
    INTEGER :: ik,ibnd, i, j, k, ir, isym, ikq, ig
    INTEGER :: h_ibnd
    INTEGER :: ibnd_loop_start, ibnd_buff_start, ibnd_buff_end
    INTEGER :: ipol, jpol
    COMPLEX(DP),ALLOCATABLE :: temppsic(:)
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:), psic_nc(:,:)
    INTEGER :: nxxs, nrxxs
#ifdef __MPI
    COMPLEX(DP),allocatable  :: temppsic_all(:),      psic_all(:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_all_nc(:,:), psic_all_nc(:,:)
#endif
    COMPLEX(DP) :: d_spin(2,2,48)
    INTEGER :: current_ik
    integer       :: find_current_k
    INTEGER :: ibnd_start_new, ibnd_end_new
    INTEGER :: ibnd_exx
    INTEGER, ALLOCATABLE :: test_igk(:)
    !<<<
    IF (use_old_exx) THEN
       CALL exxinit_old()
       RETURN
    END IF
    !>>>
    
    CALL start_clock ('exxinit')
    CALL convert_evc(1) !!!ERROR, IK NOT ASSIGNED YET
    !
    !  prepare the symmetry matrices for the spin part
    !
    IF (noncolin) THEN
       DO isym=1,nsym
          CALL find_u(sr(:,:,isym), d_spin(:,:,isym))
       ENDDO
    ENDIF

    CALL exx_fft_create()

    ! Note that nxxs is not the same as nrxxs in parallel case
    nxxs = exx_fft%dfftt%nr1x *exx_fft%dfftt%nr2x *exx_fft%dfftt%nr3x 
    nrxxs= exx_fft%dfftt%nnr

    !allocate psic_exx
    IF(.not.allocated(psic_exx))THEN
       ALLOCATE(psic_exx(nrxxs))
       psic_exx = 0.0_DP
    END IF
#ifdef __MPI
    IF (noncolin) THEN
       ALLOCATE(psic_all_nc(nxxs,npol), temppsic_all_nc(nxxs,npol) )
    ELSE IF ( .NOT. gamma_only ) THEN
       ALLOCATE(psic_all(nxxs), temppsic_all(nxxs) )
    ENDIF
#endif
    IF (noncolin) THEN
       ALLOCATE(temppsic_nc(nrxxs, npol), psic_nc(nrxxs, npol))
    ELSE IF ( .NOT. gamma_only ) THEN
       ALLOCATE(temppsic(nrxxs))
    ENDIF
    !
    IF (.not.exx_is_active()) THEN 
       !
       erfc_scrlen = get_screening_parameter()
! Debug Zhenfei Liu 9/25/2015
       beta_in_rsh = get_rsh_beta()
! DONE Debug.
       gau_scrlen = get_gau_parameter()
       exxalfa = get_exx_fraction()
       exxdiv  = exx_divergence() 
! Debug Zhenfei Liu 10/8/2015 - I think it is important to get exxalfa
! first, then execute the divergence. This is important as explicitly tested.
       !
       CALL start_exx()
    ENDIF

    IF (.NOT.allocated (wg_collect)) ALLOCATE(wg_collect(nbnd,nkstot))
    IF (npool>1) THEN
      CALL wg_all(wg_collect, wg, nkstot, nks)
    ELSE
      wg_collect = wg
    ENDIF

    IF ( .NOT. gamma_only ) CALL exx_set_symm ( )

    ! set appropriately the x_occupation and get an upperbound to the number of 
    ! bands with non zero occupation. used to distribute bands among band groups
    x_nbnd_occ = 0
    DO ik =1,nkstot
       IF(ABS(wk_collect(ik)) > eps_occ ) THEN
          x_occupation(1:nbnd,ik) = wg_collect (1:nbnd, ik) / wk_collect(ik)
          do ibnd = max(1,x_nbnd_occ), nbnd
             if (abs(x_occupation(ibnd,ik)) > eps_occ ) x_nbnd_occ = ibnd
          end do
       ELSE
          x_occupation(1:nbnd,ik) = 0._dp
       ENDIF
    ENDDO

    CALL set_egrp_indices(x_nbnd_occ,ibnd_start,ibnd_end)
    CALL init_index_over_band(inter_egrp_comm,nbnd,nbnd)

    !this will cause exxbuff to be calculated for every band
    ibnd_start_new = 1
    ibnd_end_new = nbnd

    IF ( gamma_only ) THEN
        ibnd_buff_start = ibnd_start_new/2
        IF(MOD(ibnd_start_new,2)==0) ibnd_buff_start = ibnd_buff_start -1
        !
        ibnd_buff_end = ibnd_end_new/2
        IF(MOD(ibnd_end_new,2)==1) ibnd_buff_end = ibnd_buff_end +1
    ELSE
        ibnd_buff_start = ibnd_start_new
        ibnd_buff_end   = ibnd_end_new
    ENDIF
    !<<<<<<<<<
    ibnd_start_new = iexx_start
    ibnd_end_new = iexx_end
    !>>>>>>>>>
    !
    IF (.NOT. allocated(exxbuff)) &
        ALLOCATE( exxbuff(nrxxs*npol, ibnd_buff_start:ibnd_buff_end, nkqs))
    exxbuff=(0.0_DP,0.0_DP)

    !
    !   This is parallelized over pool. Each pool computes only its k-points
    !
    IF ( nks > 1 ) REWIND( iunigk_exx )
    KPOINTS_LOOP : &
    DO ik = 1, nks
       npw = ngk (ik)
       IF ( nks > 1 ) THEN
          READ( iunigk_exx ) igk_exx
          CALL get_buffer(evc_exx, nwordwfc_exx, iunwfc_exx, ik)
       ENDIF
       !
       ! only useful for npool>1, but always work
       current_ik=find_current_k(ik, nkstot, nks)
       !
       IF_GAMMA_ONLY : & 
       IF (gamma_only) THEN
          !
          h_ibnd = ibnd_start_new/2
          !
          IF(MOD(ibnd_start_new,2)==0) THEN
             h_ibnd=h_ibnd-1
             ibnd_loop_start=ibnd_start_new-1
          ELSE
             ibnd_loop_start=ibnd_start_new
          ENDIF

          DO ibnd = ibnd_loop_start, ibnd_end_new, 2
             ibnd_exx = ibnd
             h_ibnd = h_ibnd + 1
             !
             psic(:) = ( 0._dp, 0._dp )
             !
             if ( ibnd < ibnd_end_new ) then
                DO ig=1,exx_fft%npwt
                   psic_exx(exx_fft%nlt(ig))  = evc_exx(ig,ibnd_exx)  &
                        + ( 0._dp, 1._dp ) * evc_exx(ig,ibnd_exx+1)
                   psic_exx(exx_fft%nltm(ig)) = CONJG( evc_exx(ig,ibnd_exx) ) &
                        + ( 0._dp, 1._dp ) * CONJG( evc_exx(ig,ibnd_exx+1) )
                END DO
             else
                DO ig=1,exx_fft%npwt
                   psic_exx(exx_fft%nlt (ig)) = evc_exx(ig,ibnd_exx)
                   psic_exx(exx_fft%nltm(ig)) = CONJG( evc_exx(ig,ibnd_exx) )
                END DO
             end if

             CALL invfft ('CustomWave', psic_exx, exx_fft%dfftt, is_exx=.TRUE.)

             exxbuff(1:nrxxs,h_ibnd,ik)=psic_exx(1:nrxxs)
             
          END DO
          !
       ELSE IF_GAMMA_ONLY 
          !
          IBND_LOOP_K : &
          DO ibnd = ibnd_start_new, ibnd_end_new
             !
             ibnd_exx = ibnd
             IF (noncolin) THEN
                temppsic_nc(:,:) = ( 0._dp, 0._dp )
                temppsic_nc(exx_fft%nlt(igk_exx(1:npw)),1) = evc_exx(1:npw,ibnd_exx)
                CALL invfft ('CustomWave', temppsic_nc(:,1), exx_fft%dfftt, is_exx=.TRUE.)
                temppsic_nc(exx_fft%nlt(igk_exx(1:npw)),2) = evc_exx(npwx+1:npwx+npw,ibnd_exx)
                CALL invfft ('CustomWave', temppsic_nc(:,2), exx_fft%dfftt, is_exx=.TRUE.)
             ELSE
                temppsic(:) = ( 0._dp, 0._dp )
                temppsic(exx_fft%nlt(igk_exx(1:npw))) = evc_exx(1:npw,ibnd_exx)
                CALL invfft ('CustomWave', temppsic, exx_fft%dfftt, is_exx=.TRUE.)
             ENDIF
             !
             DO ikq=1,nkqs
                !
                IF (index_xk(ikq) /= current_ik) CYCLE
                isym = abs(index_sym(ikq) )
                !
                IF (noncolin) THEN ! noncolinear
#ifdef __MPI
                   DO ipol=1,npol
                      CALL gather_grid(exx_fft%dfftt, temppsic_nc(:,ipol), temppsic_all_nc(:,ipol))
                   ENDDO
                   IF ( me_egrp == 0 ) THEN
                      psic_all_nc(:,:) = (0.0_DP, 0.0_DP)
                      DO ipol=1,npol
                         DO jpol=1,npol
                            psic_all_nc(:,ipol)=psic_all_nc(:,ipol) &
                              +  CONJG(d_spin(jpol,ipol,isym))* &
                                 temppsic_all_nc(rir(:,isym),jpol)
                         ENDDO
                      ENDDO
                   ENDIF
                   DO ipol=1,npol
                      CALL scatter_grid(exx_fft%dfftt,psic_all_nc(:,ipol), psic_nc(:,ipol))
                   ENDDO
#else
                   psic_nc(:,:) = (0._dp, 0._dp)
                   DO ipol=1,npol
                      DO jpol=1,npol
                         psic_nc(:,ipol) = psic_nc(:,ipol) + &
                              CONJG(d_spin(jpol,ipol,isym))* &
                                        temppsic_nc(rir(:,isym),jpol)
                      END DO
                   END DO
#endif
                   exxbuff(      1:  nrxxs,ibnd,ikq)=psic_nc(:,1)
                   exxbuff(nrxxs+1:2*nrxxs,ibnd,ikq)=psic_nc(:,2)
                ELSE ! noncolinear
#ifdef __MPI
                  CALL gather_grid(exx_fft%dfftt,temppsic,temppsic_all)
                  IF ( me_egrp == 0 ) &
                    psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
                  CALL scatter_grid(exx_fft%dfftt,psic_all,psic_exx)
#else
                  psic_exx(1:nrxxs) = temppsic(rir(1:nrxxs,isym))
#endif
                  IF (index_sym(ikq) < 0 ) psic_exx(1:nrxxs) = CONJG(psic_exx(1:nrxxs))
                  exxbuff(1:nrxxs,ibnd,ikq)=psic_exx(1:nrxxs)
                  !
                ENDIF ! noncolinear

             ENDDO
             !
          ENDDO &
          IBND_LOOP_K 
          !
       ENDIF & 
       IF_GAMMA_ONLY
    ENDDO &
    KPOINTS_LOOP
    !
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc, psic_nc)
#ifdef __MPI
       DEALLOCATE(temppsic_all_nc, psic_all_nc)
#endif 
    ELSE IF ( .NOT. gamma_only ) THEN
       DEALLOCATE(temppsic)
#ifdef __MPI
       DEALLOCATE(temppsic_all, psic_all)
#endif 
    ENDIF
    !
    !<<<<<<<<<<<
    !   All band groups must have the complete set of wavefunctions
    IF (negrp>1) CALL mp_sum(exxbuff, inter_egrp_comm)
    !IF (negrp>1) THEN
    !   CALL MPI_ALLGATHER( psi_gather, &
    !        lda_max_local*m, MPI_DOUBLE_COMPLEX, &
    !        psi_work, &
    !        lda_max_local*m, MPI_DOUBLE_COMPLEX, &
    !        inter_egrp_comm, ierr )
    !END IF
    !>>>>>>>>>>>
    !   All pools must have the complete set of wavefunctions (i.e. from every kpoint)
    IF (npool>1) CALL mp_sum(exxbuff, inter_pool_comm)
    !
    ! For US/PAW only: prepare space for <beta_I|phi_j> scalar products
    !
    IF(.not. allocated(becxx) .and. okvan) THEN 
        ALLOCATE(becxx(nkqs))
        DO ikq = 1,nkqs
            CALL allocate_bec_type( nkb, nbnd, becxx(ikq))
        ENDDO
    ENDIF
    ! compute <beta_I|psi_j,k+q> for the entire de-symmetrized k+q grid
    !
    CALL compute_becxx()
    !
    ! CHECKME: probably it's enough that each pool computes its own bec
    !          and then I sum them like exxbuff, but check it. In this case this
    !          call should only act when index_xk(ikq) = current_ik
    !
    ! Initialize 4-wavefunctions one-center Fock integrals
    !    \int \psi_a(r)\phi_a(r)\phi_b(r')\psi_b(r')/|r-r'|
    !
    IF(okpaw) CALL PAW_init_keeq()
    !
    CALL change_data_structure(.FALSE.)
    CALL stop_clock ('exxinit')  
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE exxinit
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE exx_set_symm ( )
    !-----------------------------------------------------------------------
    !
    ! Uses nkqs and index_sym from module exx, computes rir
    !
    USE fft_custom,           ONLY : fft_cus
    USE symm_base,            ONLY : nsym, s, sr, ftau
    !
    IMPLICIT NONE
    !
    INTEGER :: nxxs, nr1,nr2,nr3, nr1x,nr2x,nr3x
    INTEGER :: ikq, isym, i,j,k, ri,rj,rk, ir
    LOGICAL :: ispresent(nsym)
    !
    nr1 = exx_fft%dfftt%nr1
    nr2 = exx_fft%dfftt%nr2
    nr3 = exx_fft%dfftt%nr3
    nr1x= exx_fft%dfftt%nr1x
    nr2x= exx_fft%dfftt%nr2x
    nr3x= exx_fft%dfftt%nr3x
    nxxs = nr1x*nr2x*nr3x
    IF(.NOT. ALLOCATED(rir)) ALLOCATE(rir(nxxs,nsym))
    rir = 0
    ispresent(1:nsym) = .false.

    DO ikq =1,nkqs
       isym = abs(index_sym(ikq))
       IF (.not. ispresent(isym) ) THEN
          ispresent(isym) = .true.
          IF ( mod(s(2, 1, isym) * nr1, nr2) /= 0 .or. &
               mod(s(3, 1, isym) * nr1, nr3) /= 0 .or. &
               mod(s(1, 2, isym) * nr2, nr1) /= 0 .or. &
               mod(s(3, 2, isym) * nr2, nr3) /= 0 .or. &
               mod(s(1, 3, isym) * nr3, nr1) /= 0 .or. &
               mod(s(2, 3, isym) * nr3, nr2) /= 0 ) THEN
             CALL errore ('exxinit',' EXX smooth grid is not compatible with &
                                    & symmetry: change ecutfock',isym)
          ENDIF
          DO ir=1, nxxs
             rir(ir,isym) = ir
          ENDDO
          DO k = 1, nr3
             DO j = 1, nr2
                DO i = 1, nr1
                   CALL ruotaijk (s(1,1,isym), ftau(1,isym), i,j,k, nr1,nr2,nr3, ri,rj,rk)
                   ir =   i + ( j-1)*nr1x + ( k-1)*nr1x*nr2x
                   rir(ir,isym) = ri + (rj-1)*nr1x + (rk-1)*nr1x*nr2x
                ENDDO
             ENDDO
          ENDDO

       ENDIF
    ENDDO
  END SUBROUTINE exx_set_symm
  !
  !-----------------------------------------------------------------------
  SUBROUTINE compute_becxx ( )
    !-----------------------------------------------------------------------
    !
    ! prepare the necessary quantities, then call calbec to compute <beta_I|phi_j,k+q>
    ! and store it becxx(ikq). This must be called AFTER exxbuff and xkq_collected are done
    ! (i.e. at the end of exxinit)
    !
    USE kinds,                ONLY : DP
    USE wvfct,                ONLY : g2kin, npwx, ecutwfc, nbnd
    USE gvect,                ONLY : g, ngm
    USE gvecs,                ONLY : nls, nlsm
    USE cell_base,            ONLY : tpiba2
    USE uspp,                 ONLY : nkb, okvan
    USE becmod,               ONLY : calbec
    USE fft_base,             ONLY : dffts
    USE fft_interfaces,       ONLY : fwfft
    USE control_flags,        ONLY : gamma_only
    USE us_exx,               ONLY : becxx

    IMPLICIT NONE
    !
    INTEGER  :: npwq, ibnd, i, ikq, j, h_ibnd, ibnd_loop_start
    REAL(DP) :: gcutwfc
    INTEGER,ALLOCATABLE     :: igkq(:)   ! order of wavefunctions at k+q[+G]
    COMPLEX(DP),ALLOCATABLE :: vkbq(:,:) ! |beta_I> 
    COMPLEX(DP),ALLOCATABLE :: evcq(:,:) ! |psi_j,k> in g-space
    COMPLEX(DP),ALLOCATABLE :: phi(:)    ! aux space for fwfft
    COMPLEX(DP) :: fp, fm
    !
    ! NOTE: I do not want to use vkb from uspp, as you never know if it is going to be used again or not,
    !       this way we are wasting some memory, but the fault is with uspp that should not use global
    !       variables for temporary data (lp-2012-10-03)
    !
    IF(.not. okvan) RETURN
    !
    CALL start_clock('becxx')
    !
    gcutwfc = ecutwfc / tpiba2
    ALLOCATE(igkq(npwx))
    ALLOCATE(vkbq(npwx,nkb))
    ALLOCATE(phi(dffts%nnr))
    ALLOCATE(evcq(npwx,nbnd))
    !
    DO ikq = 1,nkqs
      ! each pool only does its own k-points, then it calls mp_sum (to be tested)
      ! bands count is reset at each k-point
      !
      ! prepare the g-vectors mapping
      CALL gk_sort(xkq_collect(:, ikq), ngm, g, gcutwfc, npwq, igkq, g2kin )
      ! prepare the |beta> function at k+q
      CALL init_us_2(npwq, igkq, xkq_collect(:, ikq), vkbq)
      !
      ! take rotated phi to G space
      IF (gamma_only) THEN
         !
         h_ibnd=ibnd_start/2
         !
         IF(MOD(ibnd_start,2)==0) THEN
            h_ibnd=h_ibnd-1
            ibnd_loop_start=ibnd_start-1
         ELSE
            ibnd_loop_start=ibnd_start
         ENDIF

         DO ibnd = ibnd_loop_start,ibnd_end,2
            h_ibnd = h_ibnd + 1
            phi(:) = exxbuff(:,h_ibnd,ikq)
            CALL fwfft ('Wave', phi, dffts, is_exx=.TRUE.)
            IF (ibnd < ibnd_end) THEN
               ! two ffts at the same time
               DO j = 1, npwq
                  fp = (phi (nls(igkq(j))) + phi (nlsm(igkq(j))))*0.5d0
                  fm = (phi (nls(igkq(j))) - phi (nlsm(igkq(j))))*0.5d0
                  evcq( j, ibnd)   = CMPLX( DBLE(fp), AIMAG(fm),kind=DP)
                  evcq( j, ibnd+1) = CMPLX(AIMAG(fp),- DBLE(fm),kind=DP)
               ENDDO
            ELSE
               DO j = 1, npwq
                  evcq(j, ibnd)   =  phi(nls(igkq(j)))
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO ibnd = ibnd_start,ibnd_end
            phi(:) = exxbuff(:,ibnd,ikq)
            CALL fwfft ('Wave', phi, dffts, is_exx=.TRUE.)
            FORALL(i=1:npwq) evcq(i,ibnd) = phi(nls(igkq(i)))
         ENDDO
      ENDIF
      !
      ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
      CALL calbec(npwq, vkbq, evcq, becxx(ikq), nbnd)
      !
    ENDDO
    !
    ! only work for k (only to be called once...):
    ! CALL mp_sum(becxx%k, inter_pool_comm)
    !
    DEALLOCATE(igkq, vkbq, phi, evcq)
    !
    CALL stop_clock('becxx')
    !-----------------------------------------------------------------------
  END SUBROUTINE compute_becxx
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... Wrapper routine computing V_x\psi, V_x = exchange potential
    ! ... Calls generic version vexx_k or Gamma-specific one vexx_gamma
    !
    ! ... input:
    ! ...    lda   leading dimension of arrays psi and hpsi
    ! ...    n     true dimension of psi and hpsi
    ! ...    m     number of states psi
    ! ...    psi   m wavefunctions 
    ! ..     becpsi <beta|psi>, optional but needed for US and PAW case
    !
    ! ... output:
    ! ...    hpsi  V_x*psi
    !
    USE becmod,         ONLY : bec_type
    USE control_flags,  ONLY : gamma_only
    USE uspp,           ONLY : okvan
    USE paw_variables,  ONLY : okpaw
    USE mp_exx,         ONLY : negrp, inter_egrp_comm, init_index_over_band, &
                               use_old_exx
    USE wvfct,          ONLY : nbnd
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m) 
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi
    INTEGER :: i
    !<<<
    IF ( use_old_exx ) THEN
       CALL vexx_old(lda, n, m, psi, hpsi, becpsi)
       RETURN
    END IF
    !>>>
    !
    IF ( (okvan.OR.okpaw) .AND. .NOT. PRESENT(becpsi)) &
       CALL errore('vexx','becpsi needed for US/PAW case',1)
    CALL start_clock ('vexx')
    !
    IF(gamma_only) THEN
       CALL vexx_gamma(lda, n, m, psi, hpsi, becpsi) 
    ELSE
       IF(negrp.eq.1)THEN
          CALL vexx_k_pairs(lda, n, m, psi, hpsi, becpsi)
       ELSE
          CALL init_index_over_band(inter_egrp_comm,nbnd,m)
          CALL start_exx_parallelization(lda,n,m,psi,hpsi,becpsi)
          !CALL vexx_k(lda, n, m, psi_exx, hpsi_exx, becpsi)
          CALL vexx_k_pairs(lda, n, m, psi_exx, hpsi_exx, becpsi)
          CALL end_exx_parallelization(lda,n,m,psi,hpsi)
       END IF
    ENDIF
    !
    CALL stop_clock ('vexx')
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_gamma(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... Gamma-specific version of vexx
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, npw, igk, current_k, ecutwfc
    USE control_flags,  ONLY : gamma_only
    USE klist,          ONLY : xk, nks, nkstot
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_exx,       ONLY : inter_egrp_comm, intra_egrp_comm, my_egrp_id, negrp
    USE mp,             ONLY : mp_sum, mp_barrier, mp_bcast
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    !
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m) 
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: result(:), result_g(:)
    REAL(DP),ALLOCATABLE :: temppsic_dble (:)
    REAL(DP),ALLOCATABLE :: temppsic_aimag(:)
    !
    COMPLEX(DP),ALLOCATABLE :: rhoc(:), vc(:), deexx(:)
    REAL(DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ikq, iq, ipol
    INTEGER          :: ir, ig
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    INTEGER          :: h_ibnd, nrxxs
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xkq(3)
    ! <LMS> temp array for vcut_spheric
    INTEGER  :: find_current_k
    LOGICAL :: l_fft_doubleband
    LOGICAL :: l_fft_singleband
    !
    ALLOCATE( fac(exx_fft%ngmt) )
    nrxxs= exx_fft%dfftt%nnr
    !
    ALLOCATE( result(nrxxs), temppsic_dble(nrxxs), temppsic_aimag(nrxxs) )
    ALLOCATE( result_g(n) )
    !
    ALLOCATE(rhoc(nrxxs), vc(nrxxs))
    IF(okvan) ALLOCATE(deexx(nkb))
    !
    current_ik=find_current_k(current_k,nkstot,nks)
    xkp = xk_collect(:,current_ik)
    !
    ! This is to stop numerical inconsistencies creeping in through the band parallelization.
    !
!    !IF(my_bgrp_id>0) THEN
!    IF(my_egrp_id>0) THEN
!       hpsi=0.0_DP
!       psi =0.0_DP
!    ENDIF
!    !IF (nbgrp>1) THEN
!    IF (negrp>1) THEN
!       !CALL mp_bcast(hpsi,0,inter_bgrp_comm)
!       !CALL mp_bcast(psi,0,inter_bgrp_comm)
!       CALL mp_bcast(hpsi,0,inter_egrp_comm)
!       CALL mp_bcast(psi,0,inter_egrp_comm)
!    ENDIF
    !
    ! Here the loops start
    !
    INTERNAL_LOOP_ON_Q : &
    DO iq=1,nqs
       !
       ikq  = index_xkq(current_ik,iq)
       ik   = index_xk(ikq)
       xkq  = xkq_collect(:,ikq)
       !
       ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
       CALL g2_convolution(exx_fft%ngmt, exx_fft%gt, xk(:,current_k), xkq, fac) 
       IF ( okvan .AND..NOT.tqr ) CALL qvan_init (xkq, xkp)
       !
       LOOP_ON_PSI_BANDS : &
       DO im = 1,m !for each band of psi (the k cycle is outside band)
          IF(okvan) deexx(:) = 0.0_DP
          !
          result = 0.0_DP
          !
          l_fft_doubleband = .FALSE.
          l_fft_singleband = .FALSE.
          !
          IF ( MOD(im,2)==1 .AND. (im+1)<=m ) l_fft_doubleband = .TRUE.
          IF ( MOD(im,2)==1 .AND. im==m )     l_fft_singleband = .TRUE.
          !
          IF( l_fft_doubleband ) THEN 
!$omp parallel do  default(shared), private(ig)
             DO ig = 1, exx_fft%npwt
                result( exx_fft%nlt(ig) )  =       psi(ig, im) + (0._DP,1._DP) * psi(ig, im+1)
                result( exx_fft%nltm(ig) ) = CONJG(psi(ig, im) - (0._DP,1._DP) * psi(ig, im+1))
             ENDDO
!$omp end parallel do
          ENDIF
          !
          IF( l_fft_singleband ) THEN 
!$omp parallel do  default(shared), private(ig)
             DO ig = 1, exx_fft%npwt
                result( exx_fft%nlt(ig) )  =       psi(ig,im) 
                result( exx_fft%nltm(ig) ) = CONJG(psi(ig,im))
             ENDDO
!$omp end parallel do
          ENDIF
          !
          IF( l_fft_doubleband.OR.l_fft_singleband) THEN
             CALL invfft ('CustomWave', result, exx_fft%dfftt, is_exx=.TRUE.)
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                temppsic_dble(ir)  = DBLE ( result(ir) )
                temppsic_aimag(ir) = AIMAG( result(ir) )
             ENDDO
!$omp end parallel do
          ENDIF
          !
          result = 0.0_DP
          !
          h_ibnd = ibnd_start/2
          IF(MOD(ibnd_start,2)==0) THEN
             h_ibnd=h_ibnd-1
             ibnd_loop_start=ibnd_start-1
          ELSE
             ibnd_loop_start=ibnd_start
          ENDIF
          !
          IBND_LOOP_GAM : &
          DO ibnd=ibnd_loop_start,ibnd_end, 2 !for each band of psi
             !
             h_ibnd = h_ibnd + 1
             IF( ibnd < ibnd_start ) THEN
                x1 = 0.0_DP
             ELSE
                x1 = x_occupation(ibnd,  ik)
             ENDIF
             IF( ibnd == ibnd_end) THEN
                x2 = 0.0_DP
             ELSE
                x2 = x_occupation(ibnd+1,  ik)
             ENDIF
             IF ( ABS(x1) < eps_occ .AND. ABS(x2) < eps_occ ) CYCLE
             !
             ! calculate rho in real space. Gamma tricks are used. 
             ! temppsic is real; tempphic contains one band in the real part, 
             ! another one in the imaginary part; the same applies to rhoc
             !
             IF( MOD(im,2) == 0 ) THEN 
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir) = exxbuff(ir,h_ibnd,ikq) * temppsic_aimag(ir) / omega 
                ENDDO
!$omp end parallel do
             ELSE
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir) = exxbuff(ir,h_ibnd,ikq) * temppsic_dble(ir) / omega 
                ENDDO
!$omp end parallel do
             ENDIF 
             !
             ! bring rho to G-space
             !
             !   >>>> add augmentation in REAL SPACE here
             IF(okvan .AND. tqr) THEN
                IF(ibnd>=ibnd_start) &
                CALL addusxx_r(rhoc, _CX(becxx(ikq)%r(:,ibnd)), _CX(becpsi%r(:,im)))
                IF(ibnd<ibnd_end) &
                CALL addusxx_r(rhoc,_CY(becxx(ikq)%r(:,ibnd+1)),_CX(becpsi%r(:,im)))
             ENDIF
             !
             CALL fwfft ('Custom', rhoc, exx_fft%dfftt, is_exx=.TRUE.)
             !   >>>> add augmentation in G SPACE here
             IF(okvan .AND. .NOT. TQR) THEN
                ! contribution from one band added to real (in real space) part of rhoc
                IF(ibnd>=ibnd_start) &
                   CALL addusxx_g(rhoc, xkq,  xkp, 'r', &
                   becphi_r=becxx(ikq)%r(:,ibnd), becpsi_r=becpsi%r(:,im) )
                ! contribution from following band added to imaginary (in real space) part of rhoc
                IF(ibnd<ibnd_end) &
                   CALL addusxx_g(rhoc, xkq,  xkp, 'i', &
                   becphi_r=becxx(ikq)%r(:,ibnd+1), becpsi_r=becpsi%r(:,im) )
             ENDIF
             !   >>>> charge density done
             !
             vc = 0._DP
             !
!$omp parallel do default(shared), private(ig)
             DO ig = 1, exx_fft%ngmt
                ! 
                vc(exx_fft%nlt(ig))  = fac(ig) * rhoc(exx_fft%nlt(ig)) 
                vc(exx_fft%nltm(ig)) = fac(ig) * rhoc(exx_fft%nltm(ig)) 
                !                 
             ENDDO
!$omp end parallel do
             !
             !   >>>>  compute <psi|H_fock G SPACE here
             IF(okvan .and. .not. TQR) THEN
                IF(ibnd>=ibnd_start) &
                CALL newdxx_g(vc, xkq, xkp, 'r', deexx, becphi_r=x1*becxx(ikq)%r(:,ibnd))
                IF(ibnd<ibnd_end) &
                CALL newdxx_g(vc, xkq, xkp, 'i', deexx,becphi_r=x2*becxx(ikq)%r(:,ibnd+1))
             ENDIF
             !
             !brings back v in real space
             CALL invfft ('Custom', vc, exx_fft%dfftt, is_exx=.TRUE.) 
             !
             !   >>>>  compute <psi|H_fock REAL SPACE here
             IF(okvan .and. tqr) THEN
                IF(ibnd>=ibnd_start) &
                CALL newdxx_r(vc, CMPLX(x1*becxx(ikq)%r(:,ibnd), 0.0_DP, KIND=DP), deexx)
                IF(ibnd<ibnd_end) &
                CALL newdxx_r(vc, CMPLX(0.0_DP,-x2*becxx(ikq)%r(:,ibnd+1), KIND=DP), deexx)
             ENDIF
             !
             IF(okpaw) THEN
                IF(ibnd>=ibnd_start) &
                CALL PAW_newdxx(x1/nqs, _CX(becxx(ikq)%r(:,ibnd)), _CX(becpsi%r(:,im)), deexx)
                IF(ibnd<ibnd_end) &
                CALL PAW_newdxx(x2/nqs, _CX(becxx(ikq)%r(:,ibnd+1)), _CX(becpsi%r(:,im)), deexx)
             ENDIF
             !
             ! accumulates over bands and k points
             !
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                result(ir) = result(ir)+x1* DBLE(vc(ir))* DBLE(exxbuff(ir,h_ibnd,ikq))&
                                       +x2*AIMAG(vc(ir))*AIMAG(exxbuff(ir,h_ibnd,ikq))
             ENDDO
!$omp end parallel do
             !
          ENDDO &
          IBND_LOOP_GAM
          !
          !
          IF(okvan) THEN
             CALL mp_sum(deexx,intra_egrp_comm)
             CALL mp_sum(deexx,inter_egrp_comm)
          ENDIF
          !
          !
          ! brings back result in G-space
          !
          CALL fwfft( 'CustomWave' , result, exx_fft%dfftt, is_exx=.TRUE. )
          !communicate result
          DO ig = 1, n
             result_g(ig) = result(exx_fft%nlt(igk_exx(ig)))
          END DO
          CALL mp_sum( result_g(1:n), inter_egrp_comm)
          !
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
! Debug Zhenfei Liu 9/26/2015
!             hpsi(ig,im)=hpsi(ig,im) - exxalfa*result_g(ig)
!             hpsi(ig,im)=hpsi(ig,im) - exxalfa*result(exx_fft%nlt(ig))
             hpsi(ig,im)=hpsi(ig,im) - result_g(ig)
! DONE Debug.

          ENDDO
!$omp end parallel do
          ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
          IF(okvan) CALL add_nlxx_pot (lda, hpsi(:,im), xkp, npw, igk_exx, &
                                       deexx, eps_occ, exxalfa)
! Debug Zhenfei Liu 9/26/2015: I am not implementing USPP (okvan) now.
! so that the above line will not run.
       ENDDO &
       LOOP_ON_PSI_BANDS
       IF ( okvan .AND..NOT.tqr ) CALL qvan_clean ()
       !
    ENDDO &
    INTERNAL_LOOP_ON_Q
    !  
    DEALLOCATE( result, temppsic_dble, temppsic_aimag) 
    !
    DEALLOCATE(rhoc, vc, fac )
    !
    IF(okvan) DEALLOCATE( deexx )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx_gamma
  !-----------------------------------------------------------------------
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_k_pairs(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... generic, k-point version of vexx
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, npw, igk, current_k, ecutwfc
    USE control_flags,  ONLY : gamma_only
    USE klist,          ONLY : xk, nks, nkstot
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_exx,       ONLY : inter_egrp_comm, intra_egrp_comm, my_egrp_id, &
         negrp, max_pairs, egrp_pairs
    USE mp,             ONLY : mp_sum, mp_barrier, mp_bcast
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    !
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m) 
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: temppsic(:), result(:)
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:),result_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: result_g(:)
    COMPLEX(DP),ALLOCATABLE :: result_nc_g(:,:)
    INTEGER          :: request_send, request_recv
    !
    COMPLEX(DP),ALLOCATABLE :: rhoc(:), vc(:), deexx(:)
    REAL(DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ikq, iq, ipol
    INTEGER          :: ir, ig
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    INTEGER          :: h_ibnd, nrxxs
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xkq(3)
    ! <LMS> temp array for vcut_spheric
    INTEGER  :: find_current_k
    DOUBLE PRECISION :: max, tempx
    COMPLEX(DP), ALLOCATABLE :: big_result(:,:)
    INTEGER :: ir_out, ipair, jbnd, old_ibnd
    !
    CALL start_clock ('vexx_init')
    !
    ALLOCATE( fac(exx_fft%ngmt) )
    nrxxs= exx_fft%dfftt%nnr
    !
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs,npol), result_nc(nrxxs,npol) )
       ALLOCATE( result_nc_g(n,2) )
    ELSE
       ALLOCATE( temppsic(nrxxs), result(nrxxs) )
       ALLOCATE( result_g(n) )
    ENDIF
    !
    ALLOCATE(rhoc(nrxxs), vc(nrxxs))
    IF(okvan) ALLOCATE(deexx(nkb))
    !
    current_ik=find_current_k(current_k,nkstot,nks)
    xkp = xk_collect(:,current_ik)
    !
    ! This is to stop numerical inconsistencies creeping in through the band parallelization.
    !

    allocate(big_result(n,m))
    big_result = 0.0_DP
    old_ibnd = 0

    CALL stop_clock ('vexx_init')

    DO ipair=1, max_pairs

       ibnd = egrp_pairs(1,ipair,my_egrp_id+1)
       jbnd = egrp_pairs(2,ipair,my_egrp_id+1)

       !CALL start_clock ('vexx_bcom')
       !CALL communicate_exxbuff(ipair, request_send, request_recv)
       !CALL stop_clock ('vexx_bcom')
       
       IF (ibnd.eq.0.or.ibnd.gt.m) CYCLE
       
       IF (ibnd.ne.old_ibnd) THEN

          CALL start_clock ('vexx_out1')
          
          IF(okvan) deexx = 0.0_DP
          
          IF (noncolin) THEN
             temppsic_nc = 0._DP
          ELSE
             temppsic    = 0.0_DP
          END IF

          IF (noncolin) THEN
             !
!$omp parallel do  default(shared), private(ig)
             DO ig = 1, n
                temppsic_nc(exx_fft%nlt(igk_exx(ig)),1) = psi(ig,ibnd)
             ENDDO
!$omp end parallel do
!$omp parallel do  default(shared), private(ig)
             DO ig = 1, n
                temppsic_nc(exx_fft%nlt(igk_exx(ig)),2) = psi(npwx+ig,ibnd)
             ENDDO
!$omp end parallel do
             !
             CALL invfft ('CustomWave', temppsic_nc(:,1), exx_fft%dfftt, &
                  is_exx=.TRUE.)
             CALL invfft ('CustomWave', temppsic_nc(:,2), exx_fft%dfftt, &
                  is_exx=.TRUE.)
             !
          ELSE
             !
!$omp parallel do  default(shared), private(ig)
             DO ig = 1, n
                temppsic( exx_fft%nlt(igk_exx(ig)) ) = psi(ig,ibnd)
             ENDDO
!$omp end parallel do
             !
             CALL invfft ('CustomWave', temppsic, exx_fft%dfftt, is_exx=.TRUE.)
             !
          END IF

          result = 0.0_DP
          
          old_ibnd = ibnd

          CALL stop_clock ('vexx_out1')

       END IF
       
       
       
       !----------------------------------------------------------------------!
       !INNER LOOP START
       !----------------------------------------------------------------------!
       DO iq=1, nqs
          
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          xkq  = xkq_collect(:,ikq)

          !
          ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
          CALL start_clock ('vexx_g2')
          CALL g2_convolution(exx_fft%ngmt, exx_fft%gt, xk(:,current_k), xkq, fac)
          CALL stop_clock ('vexx_g2')
          IF ( okvan .AND..NOT.tqr ) CALL qvan_init (xkq, xkp)
          !
          !
          IF ( ABS(x_occupation(jbnd,ik)) < eps_occ) CYCLE
          !
          !loads the phi from file
          !
          CALL start_clock ('vexx_rho')
          IF (noncolin) THEN
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                rhoc(ir) = ( CONJG(exxbuff(ir,jbnd,ikq))*temppsic_nc(ir,1) + &
                     CONJG(exxbuff(nrxxs+ir,jbnd,ikq))*temppsic_nc(ir,2) )/omega
             ENDDO
!$omp end parallel do
          ELSE
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                rhoc(ir)=CONJG(exxbuff(ir,jbnd,ikq))*temppsic(ir) / omega
             ENDDO
!$omp end parallel do
          ENDIF
          CALL stop_clock ('vexx_rho')

          !   >>>> add augmentation in REAL space HERE
          CALL start_clock ('vexx_augr')
          IF(okvan .AND. tqr) THEN ! augment the "charge" in real space
             CALL addusxx_r(rhoc, becxx(ikq)%k(:,jbnd), becpsi%k(:,im))
          ENDIF
          CALL stop_clock ('vexx_augr')
          !
          !   >>>> brings it to G-space
          CALL start_clock ('vexx_ffft')
          CALL fwfft('Custom', rhoc, exx_fft%dfftt, is_exx=.TRUE.)
          CALL stop_clock ('vexx_ffft')
          !
          !   >>>> add augmentation in G space HERE
          CALL start_clock ('vexx_augg')
          IF(okvan .AND. .NOT. tqr) THEN
             CALL addusxx_g(rhoc, xkq, xkp, 'c', &
                  becphi_c=becxx(ikq)%k(:,jbnd),becpsi_c=becpsi%k(:,im))
          ENDIF
          CALL stop_clock ('vexx_augg')
          !   >>>> charge done
          !
          CALL start_clock ('vexx_vc')
          vc = 0._DP
          !
!$omp parallel do default(shared), private(ig)
          DO ig = 1, exx_fft%ngmt
             vc(exx_fft%nlt(ig)) = fac(ig) * &
                  rhoc(exx_fft%nlt(ig)) * x_occupation(jbnd,ik) / nqs
          ENDDO
!$omp end parallel do
          CALL stop_clock ('vexx_vc')
          !
          ! Add ultrasoft contribution (RECIPROCAL SPACE)
          ! compute alpha_I,j,k+q = \sum_J \int <beta_J|phi_j,k+q> V_i,j,k,q Q_I,J(r) d3r
          CALL start_clock ('vexx_ultr')
          IF(okvan .AND. .NOT. tqr) THEN
             CALL newdxx_g(vc, xkq, xkp, 'c', deexx, &
                  becphi_c=becxx(ikq)%k(:,jbnd))
          ENDIF
          CALL stop_clock ('vexx_ultr')
          !
          !brings back v in real space
          CALL start_clock ('vexx_ifft')
          CALL invfft ('Custom', vc, exx_fft%dfftt, is_exx=.TRUE.)
          CALL stop_clock ('vexx_ifft')
          !
          ! Add ultrasoft contribution (REAL SPACE)
          CALL start_clock ('vexx_ultg')
          IF(okvan .AND. TQR) CALL newdxx_r(vc, becxx(ikq)%k(:,jbnd),deexx)
          CALL stop_clock ('vexx_ultg')
          !
          ! Add PAW one-center contribution
          CALL start_clock ('vexx_paw')
          IF(okpaw) THEN
             CALL PAW_newdxx(x_occupation(jbnd,ik)/nqs, becxx(ikq)%k(:,jbnd), &
                  becpsi%k(:,im), deexx)
          ENDIF
          CALL stop_clock ('vexx_paw')
          !
          !accumulates over bands and k points
          !
          CALL start_clock ('vexx_res')
          IF (noncolin) THEN
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                result_nc(ir,1)= result_nc(ir,1) + vc(ir) * exxbuff(ir,jbnd,ikq)
             ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                result_nc(ir,2)= result_nc(ir,2) + vc(ir) * exxbuff(ir+nrxxs,jbnd,ikq)
             ENDDO
!$omp end parallel do
          ELSE
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                result(ir) = result(ir) + vc(ir)*exxbuff(ir,jbnd,ikq)
             ENDDO
!$omp end parallel do
          ENDIF
          CALL stop_clock ('vexx_res')
          !
          CALL start_clock ('vexx_qcln')
          IF ( okvan .AND..NOT.tqr ) CALL qvan_clean ()
          CALL stop_clock ('vexx_qcln')
          !
       END DO
       !----------------------------------------------------------------------!
       !INNER LOOP END
       !----------------------------------------------------------------------!



       IF (ipair.eq.max_pairs.or.egrp_pairs(1,min(ipair+1,max_pairs),my_egrp_id+1).ne.ibnd) THEN
          CALL start_clock ('vexx_out2')
          !
          IF(okvan) THEN
             CALL mp_sum(deexx,intra_egrp_comm)
             CALL mp_sum(deexx,inter_egrp_comm)
          ENDIF
          !
          IF (noncolin) THEN
             !brings back result in G-space
             CALL fwfft ('CustomWave', result_nc(:,1), exx_fft%dfftt, &
                  is_exx=.TRUE.)
             CALL fwfft ('CustomWave', result_nc(:,2), exx_fft%dfftt, &
                  is_exx=.TRUE.)
          ELSE
             !
             CALL fwfft ('CustomWave', result, exx_fft%dfftt, is_exx=.TRUE.)
             DO ig = 1, n
                big_result(ig,ibnd) = big_result(ig,ibnd) + result(exx_fft%nlt(igk_exx(ig)))
             ENDDO
          ENDIF
          CALL stop_clock ('vexx_out2')
       END IF


    END DO



    !sum result
    CALL start_clock ('vexx_sum')
    CALL result_sum(n, m, big_result)
    CALL stop_clock ('vexx_sum')
    DO im=1, m
       CALL start_clock ('vexx_hpsi')
!$omp parallel do default(shared), private(ig)
       DO ig = 1, n
! Debug Zhenfei Liu 9/26/2015
!          hpsi(ig,im)=hpsi(ig,im) - exxalfa*big_result(ig,im)
!          hpsi(ig,im)=hpsi(ig,im) - exxalfa*result(exx_fft%nlt(igk(ig)))
          hpsi(ig,im)=hpsi(ig,im) - big_result(ig,im)
! DONE Debug.
       ENDDO
!$omp end parallel do
       CALL stop_clock ('vexx_hpsi')
       ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       CALL start_clock ('vexx_nloc')
       IF(okvan) CALL add_nlxx_pot (lda, hpsi(:,im), xkp, npw, igk_exx, &
            deexx, eps_occ, exxalfa)
! Debug Zhenfei Liu 9/26/2015: I am not implementing USPP (okvan) now.
       CALL stop_clock ('vexx_nloc')
    END DO

    !
    CALL start_clock ('vexx_deal')
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc, result_nc)
       DEALLOCATE(result_nc_g)
    ELSE
       DEALLOCATE(temppsic, result)
       DEALLOCATE(result_g)
    END IF
    DEALLOCATE(big_result)
    !
    DEALLOCATE(rhoc, vc, fac )
    !
    IF(okvan) DEALLOCATE( deexx)
    CALL stop_clock ('vexx_deal')

  END SUBROUTINE vexx_k_pairs









  !-----------------------------------------------------------------------
  SUBROUTINE vexx_k(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... generic, k-point version of vexx
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, npw, igk, current_k, ecutwfc
    USE control_flags,  ONLY : gamma_only
    USE klist,          ONLY : xk, nks, nkstot
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_exx,       ONLY : inter_egrp_comm, intra_egrp_comm, my_egrp_id, negrp
    USE mp,             ONLY : mp_sum, mp_barrier, mp_bcast
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    !
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m) 
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: temppsic(:), result(:)
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:),result_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: result_g(:)
    COMPLEX(DP),ALLOCATABLE :: result_nc_g(:,:)
    INTEGER          :: request_send, request_recv
    !
    COMPLEX(DP),ALLOCATABLE :: rhoc(:), vc(:), deexx(:)
    REAL(DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ikq, iq, ipol
    INTEGER          :: ir, ig
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    INTEGER          :: h_ibnd, nrxxs
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xkq(3)
    ! <LMS> temp array for vcut_spheric
    INTEGER  :: find_current_k
    DOUBLE PRECISION :: max, tempx
    COMPLEX(DP), ALLOCATABLE :: big_result(:,:)
    INTEGER :: ir_out, ipair, jbnd, old_ibnd
    !
    ALLOCATE( fac(exx_fft%ngmt) )
    nrxxs= exx_fft%dfftt%nnr
    !
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs,npol), result_nc(nrxxs,npol) )
       ALLOCATE( result_nc_g(n,2) )
    ELSE
       ALLOCATE( temppsic(nrxxs), result(nrxxs) )
       ALLOCATE( result_g(n) )
    ENDIF
    !
    ALLOCATE(rhoc(nrxxs), vc(nrxxs))
    IF(okvan) ALLOCATE(deexx(nkb))
    !
    current_ik=find_current_k(current_k,nkstot,nks)
    xkp = xk_collect(:,current_ik)
    !
    ! This is to stop numerical inconsistencies creeping in through the band parallelization.
    !
    !
    LOOP_ON_PSI_BANDS : &
    DO im = 1,m !for each band of psi (the k cycle is outside band)
       IF(okvan) deexx = 0.0_DP
       !
       IF (noncolin) THEN
          temppsic_nc = 0._DP
       ELSE
!$omp parallel do  default(shared), private(ig)
          do ig=1,nrxxs
             temppsic(ig)    = 0.0_DP
          enddo
!$omp end parallel do
       ENDIF
       !
       IF (noncolin) THEN
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic_nc(exx_fft%nlt(igk_exx(ig)),1) = psi(ig,im)
          ENDDO
!$omp end parallel do
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic_nc(exx_fft%nlt(igk_exx(ig)),2) = psi(npwx+ig,im)
          ENDDO
!$omp end parallel do
          !
          CALL invfft ('CustomWave', temppsic_nc(:,1), exx_fft%dfftt, &
               is_exx=.TRUE.)
          CALL invfft ('CustomWave', temppsic_nc(:,2), exx_fft%dfftt, &
               is_exx=.TRUE.)
          !
       ELSE
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic( exx_fft%nlt(igk_exx(ig)) ) = psi(ig,im)
          ENDDO
!$omp end parallel do
          CALL invfft ('CustomWave', temppsic, exx_fft%dfftt, is_exx=.TRUE.)
          !
       ENDIF
       !
       IF (noncolin) THEN
          result_nc = 0.0_DP
       ELSE
!$omp parallel do  default(shared), private(ig)
          do ig=1,nrxxs
             result(ig)    = 0.0_DP
          enddo
!$omp end parallel do
       ENDIF
       !
       INTERNAL_LOOP_ON_Q : &
       DO iq=1,nqs
          !
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          xkq  = xkq_collect(:,ikq)
          !
          ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
          CALL g2_convolution(exx_fft%ngmt, exx_fft%gt, xk(:,current_k), xkq, fac)
          IF ( okvan .AND..NOT.tqr ) CALL qvan_init (xkq, xkp)
          !
          IBND_LOOP_K : &
          DO ibnd=ibnd_start,ibnd_end !for each band of psi
             !
             IF ( ABS(x_occupation(ibnd,ik)) < eps_occ) CYCLE IBND_LOOP_K
             !
             !loads the phi from file
             !
             !   >>>> calculate rho in real space
             IF (noncolin) THEN
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir) = ( CONJG(exxbuff(ir,ibnd,ikq))*temppsic_nc(ir,1) + &
                                 CONJG(exxbuff(nrxxs+ir,ibnd,ikq))*temppsic_nc(ir,2) )/omega
                ENDDO
!$omp end parallel do
             ELSE
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir)=CONJG(exxbuff(ir,ibnd,ikq))*temppsic(ir) / omega
                ENDDO
!$omp end parallel do
             ENDIF
             !   >>>> add augmentation in REAL space HERE
             IF(okvan .AND. tqr) THEN ! augment the "charge" in real space
                CALL addusxx_r(rhoc, becxx(ikq)%k(:,ibnd), becpsi%k(:,im))
             ENDIF
             !
             !   >>>> brings it to G-space
             !CALL fwfft('Custom', rhoc, exx_fft%dfftt)
             CALL fwfft('Custom', rhoc, exx_fft%dfftt, is_exx=.TRUE.)
             !
             !   >>>> add augmentation in G space HERE
             IF(okvan .AND. .NOT. tqr) THEN
                CALL addusxx_g(rhoc, xkq, xkp, 'c', &
                   becphi_c=becxx(ikq)%k(:,ibnd),becpsi_c=becpsi%k(:,im))
             ENDIF
             !   >>>> charge done
             !
!$omp parallel do default(shared), private(ir)
             do ir=1,nrxxs
                vc(ir) = 0._DP
             enddo
!$omp end parallel do
             !
!$omp parallel do default(shared), private(ig)
             DO ig = 1, exx_fft%ngmt
                vc(exx_fft%nlt(ig)) = fac(ig) * &
                     rhoc(exx_fft%nlt(ig)) * x_occupation(ibnd,ik) / nqs
             ENDDO
!$omp end parallel do
             !
             ! Add ultrasoft contribution (RECIPROCAL SPACE)
             ! compute alpha_I,j,k+q = \sum_J \int <beta_J|phi_j,k+q> V_i,j,k,q Q_I,J(r) d3r
             IF(okvan .AND. .NOT. tqr) THEN
                CALL newdxx_g(vc, xkq, xkp, 'c', deexx, becphi_c=becxx(ikq)%k(:,ibnd))
             ENDIF
             !
             !brings back v in real space
             CALL invfft ('Custom', vc, exx_fft%dfftt, is_exx=.TRUE.)
             !
             ! Add ultrasoft contribution (REAL SPACE)
             IF(okvan .AND. TQR) CALL newdxx_r(vc, becxx(ikq)%k(:,ibnd),deexx)
             !
             ! Add PAW one-center contribution
             IF(okpaw) THEN
                CALL PAW_newdxx(x_occupation(ibnd,ik)/nqs, becxx(ikq)%k(:,ibnd), becpsi%k(:,im), deexx)
             ENDIF
             !
             !accumulates over bands and k points
             !
             IF (noncolin) THEN
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   result_nc(ir,1)= result_nc(ir,1) + vc(ir) * exxbuff(ir,ibnd,ikq)
                ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   result_nc(ir,2)= result_nc(ir,2) + vc(ir) * exxbuff(ir+nrxxs,ibnd,ikq)
                ENDDO
!$omp end parallel do
             ELSE
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   result(ir) = result(ir) + vc(ir)*exxbuff(ir,ibnd,ikq)
                ENDDO
!$omp end parallel do
             ENDIF
             !
          ENDDO &
          IBND_LOOP_K
          IF ( okvan .AND..NOT.tqr ) CALL qvan_clean ()
          !
       ENDDO &
       INTERNAL_LOOP_ON_Q
       !
       IF(okvan) THEN
         CALL mp_sum(deexx,intra_egrp_comm)
         CALL mp_sum(deexx,inter_egrp_comm)
       ENDIF
       !
       !brings back result in G-space
       !
       IF (noncolin) THEN
          !brings back result in G-space
          CALL fwfft ('CustomWave', result_nc(:,1), exx_fft%dfftt, is_exx=.TRUE.)
          CALL fwfft ('CustomWave', result_nc(:,2), exx_fft%dfftt, is_exx=.TRUE.)
          !communicate result
          DO ig = 1, n
             result_nc_g(ig,1) = result_nc(exx_fft%nlt(igk_exx(ig)),1)
          END DO
          DO ig = 1, n
             result_nc_g(ig,2) = result_nc(exx_fft%nlt(igk_exx(ig)),2)
          END DO
          CALL mp_sum( result_nc_g(1:n,1:2), inter_egrp_comm)
          !
          !adds it to hpsi
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(ig,im)    = hpsi(ig,im)     - exxalfa*result_nc_g(ig,1)
          ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(lda+ig,im)= hpsi(lda+ig,im) - exxalfa*result_nc_g(ig,2)
          ENDDO
!$omp end parallel do
          !
       ELSE
          !
          CALL fwfft ('CustomWave', result, exx_fft%dfftt, is_exx=.TRUE.)
          !
          !communicate result
          DO ig = 1, n
             result_g(ig) = result(exx_fft%nlt(igk_exx(ig)))
          END DO
          CALL mp_sum( result_g(1:n), inter_egrp_comm)
          !
          !adds it to hpsi
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(ig,im)=hpsi(ig,im) - exxalfa*result_g(ig)
          ENDDO
!$omp end parallel do
       ENDIF
       !
       ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF(okvan) CALL add_nlxx_pot (lda, hpsi(:,im), xkp, npw, igk_exx, &
                                       deexx, eps_occ, exxalfa)
       !
    ENDDO &
    LOOP_ON_PSI_BANDS
    !
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc, result_nc) 
    ELSE
       DEALLOCATE(temppsic, result) 
    END IF
    !
    DEALLOCATE(rhoc, vc, fac )
    !
    IF(okvan) DEALLOCATE( deexx)
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx_k
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE g2_convolution(ngm, g, xk, xkq, fac)
  !-----------------------------------------------------------------------
    ! This routine calculates the 1/|r-r'| part of the exact exchange 
    ! expression in reciprocal space (the G^-2 factor).
    ! It then regularizes it according to the specified recipe
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : tpiba, at, tpiba2
    USE constants, ONLY : fpi, e2, pi
    ! 
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN)    :: ngm   ! Number of G vectors
    REAL(DP), INTENT(IN)    :: g(3,ngm) ! Cartesian components of G vectors
    REAL(DP), INTENT(IN)    :: xk(3) ! current k vector
    REAL(DP), INTENT(IN)    :: xkq(3) ! current q vector
    !
    REAL(DP), INTENT(INOUT) :: fac(ngm) ! Calculated convolution
    !
    !Local variables
    INTEGER :: ig !Counters 
    REAL(DP) :: q(3), qq, x
    REAL(DP) :: grid_factor_track(ngm), qq_track(ngm)
    REAL(DP) :: nqhalf_dble(3)
    LOGICAL :: odg(3)
    !
    ! First the types of Coulomb potential that need q(3) and an external call
    !
    IF( use_coulomb_vcut_ws ) THEN 
       DO ig = 1, ngm 
          q(:)= ( xk(:) - xkq(:) + g(:,ig) ) * tpiba
          fac(ig) = vcut_get(vcut,q)
       ENDDO
       RETURN
    ENDIF
    !
    IF ( use_coulomb_vcut_spheric ) THEN
       DO ig = 1, ngm 
          q(:)= ( xk(:) - xkq(:) + g(:,ig) ) * tpiba
          fac(ig) = vcut_spheric_get(vcut,q)
       ENDDO
       RETURN
    ENDIF
    !
    ! Now the Coulomb potential that are computed on the fly
    !
    nqhalf_dble(1:3) = (/ DBLE(nq1)*0.5_DP, DBLE(nq2)*0.5_DP, DBLE(nq3)*0.5_DP /) 
    !
    ! Set the grid_factor_track and qq_track
    !
    IF( x_gamma_extrapolation ) THEN 
!$omp parallel do default(shared), private(ig,q,x,odg)
       DO ig = 1, ngm 
          q(:)= xk(:) - xkq(:) + g(:,ig) 
          qq_track(ig) = SUM(q(:)**2) * tpiba2
          x = (q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nqhalf_dble(1)
          odg(1) = ABS(x-NINT(x))<eps
          x = (q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nqhalf_dble(2)
          odg(2) = ABS(x-NINT(x))<eps
          x = (q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nqhalf_dble(3)
          odg(3) = ABS(x-NINT(x))<eps
          IF( ALL ( odg(:) ) ) THEN
             grid_factor_track(ig) = 0._DP ! on double grid
          ELSE
             grid_factor_track(ig) = grid_factor ! not on double grid
          ENDIF
       ENDDO
!$omp end parallel do
    ELSE
!$omp parallel do default(shared), private(ig,q)
       DO ig = 1, ngm 
          q(:)= xk(:) - xkq(:) + g(:,ig) 
          qq_track(ig) = SUM(q(:)**2) * tpiba2
       ENDDO
!$omp end parallel do
       grid_factor_track = 1._DP
    ENDIF
    !
    ! The big loop
    !
!$omp parallel do default(shared), private(ig,qq)
    DO ig=1,ngm
      !
      qq = qq_track(ig) 
      !
      IF(gau_scrlen > 0) THEN
         fac(ig)=e2*((pi/gau_scrlen)**(1.5_DP))*EXP(-qq/4._DP/gau_scrlen) * grid_factor_track(ig)
         !
      ELSE IF (qq > eps_qdiv) THEN
         !
         IF ( erfc_scrlen > 0  ) THEN
            fac(ig)=e2*fpi/qq*(1._DP-EXP(-qq/4._DP/erfc_scrlen**2)) * grid_factor_track(ig)
         ELSEIF( erf_scrlen > 0 ) THEN
            fac(ig)=e2*fpi/qq*(EXP(-qq/4._DP/erf_scrlen**2)) * grid_factor_track(ig)
         ELSE
            fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor_track(ig) ! as HARTREE
         ENDIF
         !
      ELSE
         !
         fac(ig)= - exxdiv ! or rather something ELSE (see F.Gygi)
         !
         IF ( yukawa > 0._DP.AND. .NOT. x_gamma_extrapolation ) fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
         IF( erfc_scrlen > 0._DP.AND. .NOT. x_gamma_extrapolation ) fac(ig) = fac(ig) + e2*pi/(erfc_scrlen**2)
         !
      ENDIF
      !
    ENDDO
!$omp end parallel do
  END SUBROUTINE g2_convolution
  !-----------------------------------------------------------------------
  !


  !-----------------------------------------------------------------------
  FUNCTION exxenergy ()
    !-----------------------------------------------------------------------
    ! 
    ! NB: This function is meant to give the SAME RESULT as exxenergy2.
    !     It is worth keeping it in the repository because in spite of being 
    !     slower it is a simple driver using vexx potential routine so it is 
    !     good, from time to time, to replace exxenergy2 with it to check that 
    !     everything is ok and energy and potential are consistent as they should.
    !
    USE io_files,               ONLY : iunigk_exx, iunwfc_exx, nwordwfc
    USE buffers,                ONLY : get_buffer
    USE wvfct,                  ONLY : nbnd, npwx, npw, igk, wg, current_k
    USE control_flags,          ONLY : gamma_only
    USE gvect,                  ONLY : gstart
    USE wavefunctions_module,   ONLY : evc
    USE lsda_mod,               ONLY : lsda, current_spin, isk
    USE klist,                  ONLY : ngk, nks, xk
    USE mp_pools,               ONLY : inter_pool_comm
    !<<<
    !USE mp_bands,               ONLY : intra_bgrp_comm, intra_bgrp_comm, nbgrp
    USE mp_exx,                 ONLY : intra_egrp_comm, intra_egrp_comm, negrp, &
                                       use_old_exx
    !>>>
    USE mp,                     ONLY : mp_sum
    USE becmod,                 ONLY : bec_type, allocate_bec_type, deallocate_bec_type, calbec
    USE uspp,                   ONLY : okvan,nkb,vkb

    IMPLICIT NONE

    TYPE(bec_type) :: becpsi
    REAL(DP)       :: exxenergy,  energy
    INTEGER        :: ibnd, ik
    COMPLEX(DP)    :: vxpsi ( npwx*npol, nbnd ), psi(npwx*npol,nbnd)
    COMPLEX(DP),EXTERNAL :: ZDOTC
    !<<<
    IF (use_old_exx) THEN
       exxenergy = exxenergy_old()
       RETURN
    END IF
    !>>>
    !
    exxenergy=0._dp
    
    CALL start_clock ('exxenergy')

    IF(okvan) CALL allocate_bec_type( nkb, nbnd, becpsi)
    energy = 0._dp
    
    IF ( nks > 1 ) REWIND( iunigk_exx )
    DO ik=1,nks
       current_k = ik
       IF ( lsda ) current_spin = isk(ik)
       npw = ngk (ik)
       IF ( nks > 1 ) THEN
          READ( iunigk_exx ) igk_exx
          CALL get_buffer(psi, nwordwfc_exx, iunwfc_exx, ik)
       ELSE
          psi(1:npwx*npol,1:nbnd) = evc(1:npwx*npol,1:nbnd)
       END IF
       !
       IF(okvan)THEN
          ! prepare the |beta> function at k+q
          CALL init_us_2(npw, igk_exx, xk(:,ik), vkb)
          ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
          CALL calbec(npw, vkb, psi, becpsi, nbnd)
       ENDIF
       !
       vxpsi(:,:) = (0._dp, 0._dp)
       CALL vexx(npwx,npw,nbnd,psi,vxpsi,becpsi)
       !
       DO ibnd=1,nbnd
          energy = energy + DBLE(wg(ibnd,ik) * ZDOTC(npw,psi(1,ibnd),1,vxpsi(1,ibnd),1))
          IF (noncolin) energy = energy + &
                            DBLE(wg(ibnd,ik) * ZDOTC(npw,psi(npwx+1,ibnd),1,vxpsi(npwx+1,ibnd),1))
          !
       ENDDO
       IF (gamma_only .and. gstart == 2) THEN
           DO ibnd=1,nbnd
              energy = energy - &
                       DBLE(0.5_dp * wg(ibnd,ik) * CONJG(psi(1,ibnd)) * vxpsi(1,ibnd))
           ENDDO
       ENDIF
    END DO
    !
    IF (gamma_only) energy = 2 * energy

    !<<<
    !CALL mp_sum( energy, intra_bgrp_comm)
    CALL mp_sum( energy, intra_egrp_comm)
    !>>>
    CALL mp_sum( energy, inter_pool_comm )
    IF(okvan)  CALL deallocate_bec_type(becpsi)
    ! 
    exxenergy = energy
    !
    CALL stop_clock ('exxenergy')
    !-----------------------------------------------------------------------
  END FUNCTION exxenergy
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  FUNCTION exxenergy2()
    !-----------------------------------------------------------------------
    !
    USE control_flags,           ONLY : gamma_only
    USE mp_exx,                  ONLY : use_old_exx
    !
    IMPLICIT NONE
    !
    REAL(DP)   :: exxenergy2
    !<<<
    IF( use_old_exx ) THEN
       exxenergy2 = exxenergy2_old()
       RETURN
    END IF
    !>>>
    !
    CALL start_clock ('exxenergy')
    !
    IF( gamma_only ) THEN 
       exxenergy2 = exxenergy2_gamma() 
    ELSE
       exxenergy2 = exxenergy2_k() 
    ENDIF
    !
    CALL stop_clock ('exxenergy')
    !
    !-----------------------------------------------------------------------
  END FUNCTION  exxenergy2
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION exxenergy2_gamma()
    !-----------------------------------------------------------------------
    !
    USE constants,               ONLY : fpi, e2, pi
    USE io_files,                ONLY : iunigk_exx, iunwfc_exx, nwordwfc
    USE buffers,                 ONLY : get_buffer
    USE cell_base,               ONLY : alat, omega, bg, at, tpiba
    USE symm_base,               ONLY : nsym, s
    USE gvect,                   ONLY : ngm, gstart, g, nl
    USE wvfct,                   ONLY : nbnd, npwx, npw, igk, wg, ecutwfc
    USE control_flags,           ONLY : gamma_only
    USE wavefunctions_module,    ONLY : evc
    USE klist,                   ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,                ONLY : lsda, current_spin, isk
    USE mp_pools,                ONLY : inter_pool_comm
    !<<<
    !USE mp_bands,                ONLY : inter_bgrp_comm, intra_bgrp_comm, nbgrp
    USE mp_exx,                ONLY : inter_egrp_comm, intra_egrp_comm, negrp
    !>>>
    USE mp,                      ONLY : mp_sum
    USE fft_interfaces,          ONLY : fwfft, invfft
    USE gvect,                   ONLY : ecutrho
    USE klist,                   ONLY : wk
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE becmod,                  ONLY : bec_type, allocate_bec_type, deallocate_bec_type, calbec
    USE paw_variables,           ONLY : okpaw
    USE paw_exx,                 ONLY : PAW_xx_energy
    USE us_exx,                  ONLY : bexg_merge, becxx, addusxx_g, &
                                        addusxx_r, qvan_init, qvan_clean
    !
    IMPLICIT NONE
    !
    REAL(DP)   :: exxenergy2_gamma
    !
    ! local variables
    REAL(DP) :: energy 
    COMPLEX(DP), ALLOCATABLE :: temppsic(:)
    COMPLEX(DP), ALLOCATABLE :: rhoc(:)
    REAL(DP),    ALLOCATABLE :: fac(:)
    INTEGER  :: jbnd, ibnd, ik, ikk, ig, ikq, iq, ir
    INTEGER  :: h_ibnd, nrxxs, current_ik, ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: xkq(3), xkp(3), vc
    ! temp array for vcut_spheric
    INTEGER,        EXTERNAL :: find_current_k
    !
    TYPE(bec_type) :: becpsi
    COMPLEX(DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    INTEGER,     ALLOCATABLE :: igkt(:)
    REAL(DP),ALLOCATABLE :: temppsic_dble (:)
    REAL(DP),ALLOCATABLE :: temppsic_aimag(:)
    LOGICAL :: l_fft_doubleband
    LOGICAL :: l_fft_singleband
    INTEGER :: jmax
    CALL convert_evc(1) !!!ERROR, IK NOT ASSIGNED
    !
    nrxxs= exx_fft%dfftt%nnr
    ALLOCATE( fac(exx_fft%ngmt) )
    !
    ALLOCATE(temppsic(nrxxs), temppsic_dble(nrxxs),temppsic_aimag(nrxxs)) 
    ALLOCATE( rhoc(nrxxs) )
    !
    energy=0.0_DP
    !
    CALL allocate_bec_type( nkb, nbnd, becpsi)
    !
    IF ( nks > 1 ) REWIND( iunigk_exx )
    !
    IKK_LOOP : &
    DO ikk=1,nks
       current_ik=find_current_k(ikk,nkstot,nks)
       xkp = xk_collect(:,current_ik)
       !
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) THEN
          READ( iunigk_exx ) igk_exx
          CALL get_buffer (evc_exx, nwordwfc_exx, iunwfc_exx, ikk)
       END IF
       !
       ! prepare the |beta> function at k+q
       !<<<
       !CALL init_us_2(npw, igk, xk(:,ikk), vkb)
       CALL init_us_2(npw, igk_exx, xk(:,ikk), vkb_exx)
       !>>>
       ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
       !<<<
       !CALL calbec(npw, vkb, evc, becpsi, nbnd)
       CALL calbec(npw, vkb_exx, evc_exx, becpsi, nbnd)
       !>>>
       !
       IQ_LOOP : &
       DO iq = 1,nqs
          !
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          !
          xkq = xkq_collect(:,ikq)
          !
          CALL g2_convolution(exx_fft%ngmt, exx_fft%gt, xk(:,current_ik), xkq, fac) 
          fac(exx_fft%gstart_t:) = 2 * fac(exx_fft%gstart_t:)
          IF ( okvan .AND..NOT.tqr ) CALL qvan_init (xkq, xkp)
          !
          jmax = nbnd 
          DO jbnd = nbnd,1, -1
             IF ( ABS(wg(jbnd,ikk)) < eps_occ) CYCLE
             jmax = jbnd 
             EXIT
          ENDDO
          !
          JBND_LOOP : &
          DO jbnd = 1, jmax     !for each band of psi (the k cycle is outside band)
             !
             temppsic = 0._DP
             !
             l_fft_doubleband = .FALSE.
             l_fft_singleband = .FALSE.
             !
             IF ( MOD(jbnd,2)==1 .AND. (jbnd+1)<=jmax ) l_fft_doubleband = .TRUE.
             IF ( MOD(jbnd,2)==1 .AND. jbnd==jmax )     l_fft_singleband = .TRUE.
             !
             IF( l_fft_doubleband ) THEN 
!$omp parallel do  default(shared), private(ig)
                DO ig = 1, exx_fft%npwt
                   !<<<
                   !temppsic( exx_fft%nlt(ig) )  =       evc(ig,jbnd) + (0._DP,1._DP) * evc(ig,jbnd+1)
                   !temppsic( exx_fft%nltm(ig) ) = CONJG(evc(ig,jbnd) - (0._DP,1._DP) * evc(ig,jbnd+1))
                   temppsic( exx_fft%nlt(ig) )  =       evc_exx(ig,jbnd) + (0._DP,1._DP) * evc_exx(ig,jbnd+1)
                   temppsic( exx_fft%nltm(ig) ) = CONJG(evc_exx(ig,jbnd) - (0._DP,1._DP) * evc_exx(ig,jbnd+1))
                   !>>>
                ENDDO
!$omp end parallel do
             ENDIF
             !
             IF( l_fft_singleband ) THEN 
!$omp parallel do  default(shared), private(ig)
                DO ig = 1, exx_fft%npwt
                   !<<<
                   !temppsic( exx_fft%nlt(ig) )  =       evc(ig,jbnd) 
                   !temppsic( exx_fft%nltm(ig) ) = CONJG(evc(ig,jbnd))
                   temppsic( exx_fft%nlt(ig) )  =       evc_exx(ig,jbnd) 
                   temppsic( exx_fft%nltm(ig) ) = CONJG(evc_exx(ig,jbnd))
                   !>>>
                ENDDO
!$omp end parallel do
             ENDIF
             !
             IF( l_fft_doubleband.OR.l_fft_singleband) THEN
                !<<<
                !CALL invfft ('CustomWave', temppsic, exx_fft%dfftt)
                CALL invfft ('CustomWave', temppsic, exx_fft%dfftt, is_exx=.TRUE.)
                !>>>
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   temppsic_dble(ir)  = DBLE ( temppsic(ir) )
                   temppsic_aimag(ir) = AIMAG( temppsic(ir) )
                ENDDO
!$omp end parallel do
             ENDIF
             !
             h_ibnd = ibnd_start/2
             IF(MOD(ibnd_start,2)==0) THEN
                h_ibnd=h_ibnd-1
                ibnd_loop_start=ibnd_start-1
             ELSE
                ibnd_loop_start=ibnd_start
             ENDIF
             !
             IBND_LOOP_GAM : &
             DO ibnd = ibnd_loop_start, ibnd_end, 2       !for each band of psi
                !
                h_ibnd = h_ibnd + 1
                !
                IF ( ibnd < ibnd_start ) THEN
                   x1 = 0.0_DP
                ELSE
                   x1 = x_occupation(ibnd,ik)
                ENDIF
                !
                IF ( ibnd < ibnd_end ) THEN
                   x2 = x_occupation(ibnd+1,ik)
                ELSE
                   x2 = 0.0_DP
                ENDIF
                IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE IBND_LOOP_GAM
                ! calculate rho in real space. Gamma tricks are used. 
                ! temppsic is real; tempphic contains band 1 in the real part, 
                ! band 2 in the imaginary part; the same applies to rhoc
                !
                IF( MOD(jbnd,2) == 0 ) THEN
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      rhoc(ir) = exxbuff(ir,h_ibnd,ikq) * temppsic_aimag(ir) / omega
                   ENDDO
!$omp end parallel do
                ELSE
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      rhoc(ir) = exxbuff(ir,h_ibnd,ikq) * temppsic_dble(ir) / omega
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                IF(okvan .and.tqr) THEN
                   IF(ibnd>=ibnd_start) &
                   CALL addusxx_r(rhoc, _CX(becxx(ikq)%r(:,ibnd)), _CX(becpsi%r(:,jbnd)))
                   IF(ibnd<ibnd_end) &
                   CALL addusxx_r(rhoc,_CY(becxx(ikq)%r(:,ibnd+1)),_CX(becpsi%r(:,jbnd)))
                ENDIF
                !
                ! bring rhoc to G-space
                !<<<
                !CALL fwfft ('Custom', rhoc, exx_fft%dfftt)
                CALL fwfft ('Custom', rhoc, exx_fft%dfftt, is_exx=.TRUE.)
                !>>>
                !
                IF(okvan .and..not.tqr) THEN
                   IF(ibnd>=ibnd_start ) &
                      CALL addusxx_g( rhoc, xkq, xkp, 'r', &
                      becphi_r=becxx(ikq)%r(:,ibnd), becpsi_r=becpsi%r(:,jbnd) )
                   IF(ibnd<ibnd_end) &
                      CALL addusxx_g( rhoc, xkq, xkp, 'i', &
                      becphi_r=becxx(ikq)%r(:,ibnd+1), becpsi_r=becpsi%r(:,jbnd) )
                ENDIF
                !
                vc = 0.0_DP
!$omp parallel do  default(shared), private(ig),  reduction(+:vc)
                DO ig = 1,exx_fft%ngmt
                   !
                   ! The real part of rhoc contains the contribution from band ibnd
                   ! The imaginary part    contains the contribution from band ibnd+1
                   !
                   vc = vc + fac(ig) * ( x1 * &
                        ABS( rhoc(exx_fft%nlt(ig)) + CONJG(rhoc(exx_fft%nltm(ig))) )**2 &
                                        +x2 * &
                        ABS( rhoc(exx_fft%nlt(ig)) - CONJG(rhoc(exx_fft%nltm(ig))) )**2 )
                ENDDO
!$omp end parallel do
                !
                vc = vc * omega * 0.25_DP / nqs
! Debug Zhenfei Liu 9/26/2015
!                energy = energy - exxalfa * vc * wg(jbnd,ikk)
                energy = energy - vc * wg(jbnd,ikk)
! DONE Debug.
                !
                IF(okpaw) THEN
                   IF(ibnd>=ibnd_start) &
! Debug Zhenfei Liu 9/26/2015
!                   energy = energy +exxalfa*wg(jbnd,ikk)*&
                   energy = energy + wg(jbnd,ikk)*&
! DONE Debug.
                         x1 * PAW_xx_energy(_CX(becxx(ikq)%r(:,ibnd)),_CX(becpsi%r(:,jbnd)) )
                   IF(ibnd<ibnd_end) &
! Debug Zhenfei Liu 9/26/2015
!                   energy = energy +exxalfa*wg(jbnd,ikk)*&
                   energy = energy + wg(jbnd,ikk)*&
! DONE Debug.
                         x2 * PAW_xx_energy(_CX(becxx(ikq)%r(:,ibnd+1)), _CX(becpsi%r(:,jbnd)) ) 
                ENDIF
                !
             ENDDO &
             IBND_LOOP_GAM
          ENDDO &
          JBND_LOOP
          IF ( okvan .AND..NOT.tqr ) CALL qvan_clean ( )
          !
       ENDDO &
       IQ_LOOP
    ENDDO &
    IKK_LOOP
    !
    DEALLOCATE(temppsic,temppsic_dble,temppsic_aimag) 
    !
    DEALLOCATE(rhoc, fac )
    CALL deallocate_bec_type(becpsi)
    !
    !<<<
    !CALL mp_sum( energy, inter_bgrp_comm )
    !CALL mp_sum( energy, intra_bgrp_comm )
    CALL mp_sum( energy, inter_egrp_comm )
    CALL mp_sum( energy, intra_egrp_comm )
    !>>>
    CALL mp_sum( energy, inter_pool_comm )
    !
    exxenergy2_gamma = energy
    CALL change_data_structure(.FALSE.)
    !
    !-----------------------------------------------------------------------
  END FUNCTION  exxenergy2_gamma
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION exxenergy2_k()
    !-----------------------------------------------------------------------
    !
    USE constants,               ONLY : fpi, e2, pi
    USE io_files,                ONLY : iunigk_exx, iunwfc_exx, nwordwfc
    USE buffers,                 ONLY : get_buffer
    USE cell_base,               ONLY : alat, omega, bg, at, tpiba
    USE symm_base,               ONLY : nsym, s
    USE gvect,                   ONLY : ngm, gstart, g, nl
    USE wvfct,                   ONLY : nbnd, npwx, npw, igk, wg, ecutwfc
    USE control_flags,           ONLY : gamma_only
    USE wavefunctions_module,    ONLY : evc
    USE klist,                   ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,                ONLY : lsda, current_spin, isk
    USE mp_pools,                ONLY : inter_pool_comm
    !<<<
    !USE mp_bands,                ONLY : inter_bgrp_comm, intra_bgrp_comm, nbgrp
    USE mp_exx,                ONLY : inter_egrp_comm, intra_egrp_comm, negrp
    !>>>
    USE mp,                      ONLY : mp_sum
    USE fft_interfaces,          ONLY : fwfft, invfft
    USE gvect,                   ONLY : ecutrho
    USE klist,                   ONLY : wk
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE becmod,                  ONLY : bec_type, allocate_bec_type, deallocate_bec_type, calbec
    USE paw_variables,           ONLY : okpaw
    USE paw_exx,                 ONLY : PAW_xx_energy
    USE us_exx,                  ONLY : bexg_merge, becxx, addusxx_g, &
                                        addusxx_r, qvan_init, qvan_clean
    !
    IMPLICIT NONE
    !
    REAL(DP)   :: exxenergy2_k
    !
    ! local variables
    REAL(DP) :: energy 
    COMPLEX(DP), ALLOCATABLE :: temppsic(:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_nc(:,:)
    COMPLEX(DP), ALLOCATABLE :: rhoc(:)
    REAL(DP),    ALLOCATABLE :: fac(:)
    INTEGER  :: jbnd, ibnd, ik, ikk, ig, ikq, iq, ir
    INTEGER  :: h_ibnd, nrxxs, current_ik, ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: xkq(3), xkp(3), vc
    ! temp array for vcut_spheric
    INTEGER,        EXTERNAL :: find_current_k
    !
    TYPE(bec_type) :: becpsi
    COMPLEX(DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    INTEGER,     ALLOCATABLE :: igkt(:)
    INTEGER :: i
    CALL convert_evc(1) !!!ERROR, IK NOT ASSIGNED
    !
    nrxxs = exx_fft%dfftt%nnr
    ALLOCATE( fac(exx_fft%ngmt) )
    !
    IF (noncolin) THEN
       ALLOCATE(temppsic_nc(nrxxs,npol))
    ELSE
       ALLOCATE(temppsic(nrxxs)) 
    ENDIF
    ALLOCATE( rhoc(nrxxs) )
    !
    energy=0.0_DP
    !
    CALL allocate_bec_type( nkb, nbnd, becpsi)
    !
    IF ( nks > 1 ) REWIND( iunigk_exx )
    !
    IKK_LOOP : &
    DO ikk=1,nks
       current_ik=find_current_k(ikk,nkstot,nks)
       xkp = xk_collect(:,current_ik)
       !
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) THEN
          READ( iunigk_exx ) igk_exx
          CALL get_buffer (evc_exx, nwordwfc_exx, iunwfc_exx, ikk)
       END IF
       !
       ! prepare the |beta> function at k+q
       CALL init_us_2(npw, igk_exx, xk(:,ikk), vkb_exx)
       ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
       CALL calbec(npw, vkb_exx, evc_exx, becpsi, nbnd)
       !
       JBND_LOOP : &
       DO jbnd = 1, nbnd     !for each band of psi (the k cycle is outside band)
          !
          IF ( ABS(wg(jbnd,ikk)) < eps_occ) CYCLE
          !
          IF (noncolin) THEN
             temppsic_nc = 0.0_DP
          ELSE
             temppsic    = 0.0_DP
          ENDIF
          !
          IF (noncolin) THEN
             !
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic_nc(exx_fft%nlt(igk_exx(ig)),1) = evc_exx(ig,jbnd)
             ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic_nc(exx_fft%nlt(igk_exx(ig)),2) = evc_exx(npwx+ig,jbnd)
             ENDDO
!$omp end parallel do
             !
             CALL invfft ('CustomWave', temppsic_nc(:,1), exx_fft%dfftt, is_exx=.TRUE.)
             CALL invfft ('CustomWave', temppsic_nc(:,2), exx_fft%dfftt, is_exx=.TRUE.)
             !
          ELSE
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic(exx_fft%nlt(igk_exx(ig))) = evc_exx(ig,jbnd)
             ENDDO
!$omp end parallel do
             !
             CALL invfft ('CustomWave', temppsic, exx_fft%dfftt, is_exx=.TRUE.)
             !
          ENDIF
          ! 
          IQ_LOOP : &
          DO iq = 1,nqs
             !
             ikq  = index_xkq(current_ik,iq)
             ik   = index_xk(ikq)
             !
             xkq = xkq_collect(:,ikq)
             !
             CALL g2_convolution(exx_fft%ngmt, exx_fft%gt, xk(:,current_ik), xkq, fac)
             IF ( okvan .AND..NOT.tqr ) CALL qvan_init (xkq, xkp)
             !
             IBND_LOOP_K : &
             DO ibnd = ibnd_start, ibnd_end
                !
                IF ( ABS(x_occupation(ibnd,ik)) < eps_occ) CYCLE
                !
                ! load the phi at this k+q and band
                IF (noncolin) THEN
                   !
!$omp parallel do  default(shared), private(ir) 
                   DO ir = 1, nrxxs
                      rhoc(ir)=(CONJG(exxbuff(ir      ,ibnd,ikq))*temppsic_nc(ir,1) + &
                                CONJG(exxbuff(ir+nrxxs,ibnd,ikq))*temppsic_nc(ir,2) )/omega
                   ENDDO
!$omp end parallel do
                ELSE
                   !calculate rho in real space
!$omp parallel do  default(shared), private(ir)
                   DO ir = 1, nrxxs
                      rhoc(ir)=CONJG(exxbuff(ir,ibnd,ikq))*temppsic(ir) / omega
                   ENDDO
!$omp end parallel do
                ENDIF
                ! augment the "charge" in real space
                IF(okvan .AND. tqr) CALL addusxx_r(rhoc, becxx(ikq)%k(:,ibnd), becpsi%k(:,jbnd))
                !
                ! bring rhoc to G-space
                CALL fwfft ('Custom', rhoc, exx_fft%dfftt, is_exx=.TRUE.)
                ! augment the "charge" in G space
                IF(okvan .AND. .NOT. tqr) & 
                   CALL addusxx_g(rhoc, xkq, xkp, 'c', &
                   becphi_c=becxx(ikq)%k(:,ibnd),becpsi_c=becpsi%k(:,jbnd))
                !
                vc = 0.0_DP
!$omp parallel do  default(shared), private(ig), reduction(+:vc)
                DO ig=1,exx_fft%ngmt
                   vc = vc + fac(ig) * DBLE(rhoc(exx_fft%nlt(ig)) * &
                                      CONJG(rhoc(exx_fft%nlt(ig))))
                ENDDO
!$omp end parallel do
                vc = vc * omega * x_occupation(ibnd,ik) / nqs
                ! 
! Debug Zhenfei Liu 9/26/2015
!                energy = energy - exxalfa * vc * wg(jbnd,ikk)
                energy = energy - vc * wg(jbnd,ikk)
! DONE Debug.
                !
                IF(okpaw) THEN
! Debug Zhenfei Liu 9/26/2015
!                   energy = energy +exxalfa*x_occupation(ibnd,ik)/nqs*wg(jbnd,ikk) &
                   energy = energy + x_occupation(ibnd,ik)/nqs*wg(jbnd,ikk) &
! DONE Debug.
                              *PAW_xx_energy(becxx(ikq)%k(:,ibnd), becpsi%k(:,jbnd))
                ENDIF
                !
             ENDDO &
             IBND_LOOP_K 
             IF ( okvan .AND..NOT.tqr ) CALL qvan_clean ( )
          ENDDO &
          IQ_LOOP
       ENDDO &
       JBND_LOOP
    ENDDO &
    IKK_LOOP
    !
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc) 
    ELSE
       DEALLOCATE(temppsic) 
    ENDIF
    !
    DEALLOCATE(rhoc, fac )
    CALL deallocate_bec_type(becpsi)
    !
    CALL mp_sum( energy, inter_egrp_comm )
    CALL mp_sum( energy, intra_egrp_comm )
    CALL mp_sum( energy, inter_pool_comm )
    !
    exxenergy2_k = energy
    CALL change_data_structure(.FALSE.)
    !
    !-----------------------------------------------------------------------
  END FUNCTION  exxenergy2_k
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION exx_divergence ()
    !-----------------------------------------------------------------------
     USE constants,      ONLY : fpi, e2, pi
     USE cell_base,      ONLY : bg, at, alat, omega
     USE gvect,          ONLY : ngm, g
     USE wvfct,          ONLY : ecutwfc
     USE io_global,      ONLY : stdout
     USE control_flags,  ONLY : gamma_only
     USE mp_exx,         ONLY : intra_egrp_comm, use_old_exx
     USE mp,             ONLY : mp_sum

     IMPLICIT NONE
     REAL(DP) :: exx_divergence

     ! local variables
     INTEGER :: iq1,iq2,iq3, ig
     REAL(DP) :: div, dq1, dq2, dq3, xq(3), q_, qq, tpiba2, alpha, x, q(3)

     INTEGER :: nqq, iq
     REAL(DP) :: aa, dq

     !<<<
     IF( use_old_exx ) THEN
        exx_divergence = exx_divergence_old()
        RETURN
     END IF
     !>>>

     CALL start_clock ('exx_div')

     tpiba2 = (fpi / 2.d0 / alat) **2

     alpha  = 10._dp * tpiba2 / ecutwfc

     IF ( .NOT. use_regularization ) THEN
        exx_divergence = 0._dp
        RETURN
     END IF

     dq1= 1._dp/DBLE(nq1)
     dq2= 1._dp/DBLE(nq2) 
     dq3= 1._dp/DBLE(nq3) 

     div = 0._dp
     DO iq1=1,nq1
        DO iq2=1,nq2
           DO iq3=1,nq3
              xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                      bg(:,2) * (iq2-1) * dq2 + &
                      bg(:,3) * (iq3-1) * dq3 
              DO ig=1,ngm
                 q(1)= xq(1) + g(1,ig)
                 q(2)= xq(2) + g(2,ig)
                 q(3)= xq(3) + g(3,ig)
                 qq = ( q(1)**2 + q(2)**2 + q(3)**2 ) 
                 IF (x_gamma_extrapolation) THEN
                    on_double_grid = .true.
                    x= 0.5d0*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                 ENDIF
                 IF (.not.on_double_grid) THEN
                    IF ( qq > 1.d-8 ) THEN
! Debug Zhenfei Liu 9/27/2015
! this equals: (alpha+beta)*1 - beta*(1-exp[])
                          div = div + exp( -alpha * qq) / qq * &
                                ( exxalfa + beta_in_rsh * &
                                     exp(-qq*tpiba2/4.d0/erfc_scrlen**2) ) * &
                                grid_factor
!                       IF ( erfc_scrlen > 0 ) THEN
!                          div = div + exp( -alpha * qq) / qq * &
!                                (1._dp-exp(-qq*tpiba2/4.d0/erfc_scrlen**2)) * grid_factor
!                       ELSEIF ( erf_scrlen >0 ) THEN
!                          div = div + exp( -alpha * qq) / qq * &
!                                (exp(-qq*tpiba2/4.d0/erf_scrlen**2)) * grid_factor
!                       ELSE
!                          div = div + exp( -alpha * qq) / (qq + yukawa/tpiba2) &
!                                                     * grid_factor
!                       ENDIF
! DONE Debug.
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     CALL mp_sum(  div, intra_egrp_comm )
     IF (gamma_only) THEN
        div = 2.d0 * div
     ENDIF
! Debug Zhenfei Liu 09/26/2015: will enforce x_gamma_extrapolation (default)
     IF ( .not. x_gamma_extrapolation ) THEN
        IF ( yukawa > 0._dp) THEN
           div = div + tpiba2/yukawa
        ELSEIF( erfc_scrlen > 0._dp ) THEN
           div = div + tpiba2/4.d0/erfc_scrlen**2
        ELSE
           div = div - alpha
        ENDIF
     ENDIF

     div = div * e2 * fpi / tpiba2 / nqs

     alpha = alpha / tpiba2

     nqq = 100000
     dq = 5.0d0 / sqrt(alpha) /nqq
     aa = 0._dp
     DO iq=0,  nqq
        q_ = dq * (iq+0.5d0)
        qq = q_ * q_
! Debug Zhenfei Liu 9/28/2015
! this equals: (alpha+beta)*0 - beta*[]
           aa = aa + beta_in_rsh * &
                     exp(-alpha*qq)*exp(-qq/4.d0/erfc_scrlen**2) * dq
!        IF ( erfc_scrlen > 0 ) THEN
!           aa = aa  -exp( -alpha * qq) * exp(-qq/4.d0/erfc_scrlen**2) * dq
!        ELSEIF ( erf_scrlen > 0 ) THEN
!           aa = 0._dp
!        ELSE
!           aa = aa - exp( -alpha * qq) * yukawa / (qq + yukawa) * dq
!        ENDIF
! DONE Debug.
     ENDDO
     aa = aa * 8.d0 /fpi
! Debug Zhenfei Liu 10/8/2015: added exxalfa below
!     aa = aa + 1._dp/sqrt(alpha*0.25d0*fpi)
     aa = aa + 1._dp/sqrt(alpha*0.25d0*fpi) * exxalfa
! DONE Debug.
     if( erf_scrlen > 0) aa = 1._dp/sqrt((alpha+1._dp/4.d0/erf_scrlen**2)*0.25d0*fpi)
     div = div - e2*omega * aa

     exx_divergence = div * nqs
     CALL stop_clock ('exx_div')

     return
    !-----------------------------------------------------------------------
  END FUNCTION exx_divergence 
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION exx_stress()
    !-----------------------------------------------------------------------
    !
    ! This is Eq.(10) of PRB 73, 125120 (2006).
    !
    USE constants,            ONLY : fpi, e2, pi, tpi
    USE io_files,             ONLY : iunigk_exx, iunwfc_exx, nwordwfc
    USE buffers,              ONLY : get_buffer
    USE cell_base,            ONLY : alat, omega, bg, at, tpiba
    USE symm_base,            ONLY : nsym, s
    USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, current_k
    USE control_flags,        ONLY : gamma_only
    USE wavefunctions_module, ONLY : evc
    USE klist,                ONLY : xk, ngk, nks
    USE lsda_mod,             ONLY : lsda, current_spin, isk
    USE gvect,                ONLY : g, nl
    USE mp_pools,             ONLY : npool, inter_pool_comm
    USE mp_exx,               ONLY : inter_egrp_comm, intra_egrp_comm, &
                                     use_old_exx
    USE mp,                   ONLY : mp_sum 
    USE fft_base,             ONLY : dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE uspp,                 ONLY : okvan
    !
    ! ---- local variables -------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! local variables
    REAL(DP)   :: exx_stress(3,3), exx_stress_(3,3)
    !
    COMPLEX(DP),ALLOCATABLE :: tempphic(:), temppsic(:), result(:)
    COMPLEX(DP),ALLOCATABLE :: tempphic_nc(:,:), temppsic_nc(:,:), &
                               result_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: rhoc(:)
    REAL(DP),    allocatable :: fac(:), fac_tens(:,:,:), fac_stress(:)
    INTEGER  :: jbnd, ibnd, ik, ikk, ig, ir, ikq, iq, isym
    INTEGER  :: h_ibnd, nqi, iqi, beta, nrxxs, ngm
    INTEGER  :: ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), vc(3,3), x, q(3)
    ! temp array for vcut_spheric
    REAL(DP) :: delta(3,3)
    
    !<<<
    IF ( use_old_exx ) THEN
       exx_stress = exx_stress_old()
       RETURN
    END IF
    !>>>
  
    CALL start_clock ('exx_stress')

    IF (npool>1) CALL errore('exx_stress','stress not available with pools',1)
    IF (noncolin) CALL errore('exx_stress','noncolinear stress not implemented',1)
    IF (okvan) CALL infomsg('exx_stress','USPP stress not tested')

    nrxxs = exx_fft%dfftt%nnr
    ngm   = exx_fft%ngmt
    delta = reshape( (/1._dp,0._dp,0._dp, 0._dp,1._dp,0._dp, 0._dp,0._dp,1._dp/), (/3,3/))
    exx_stress_ = 0._dp
    allocate( tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
    allocate( fac_tens(3,3,ngm), fac_stress(ngm) )

    IF ( nks > 1 ) rewind( iunigk_exx )
    !
    nqi=nqs
    !
    ! loop over k-points
    DO ikk = 1, nks
        current_k = ikk
        IF (lsda) current_spin = isk(ikk)
        npw = ngk(ikk)

        IF (nks > 1) THEN
            read(iunigk_exx) igk_exx
            CALL get_buffer(evc_exx, nwordwfc_exx, iunwfc_exx, ikk)
        ENDIF

        ! loop over bands
        DO jbnd = 1, nbnd
            !
            temppsic(:) = ( 0._dp, 0._dp )
!$omp parallel do default(shared), private(ig)
            DO ig = 1, npw
                temppsic(exx_fft%nlt(igk_exx(ig))) = evc(ig,jbnd)
            ENDDO
!$omp end parallel do
            !
            IF(gamma_only) THEN
!$omp parallel do default(shared), private(ig)
                DO ig = 1, npw
                    temppsic(exx_fft%nltm(igk_exx(ig))) = conjg(evc(ig,jbnd))
                ENDDO
!$omp end parallel do
            ENDIF

            !<<<
            !CALL invfft ('CustomWave', temppsic, exx_fft%dfftt)       
            CALL invfft ('CustomWave', temppsic, exx_fft%dfftt, is_exx=.TRUE.)
            !>>>

            DO iqi = 1, nqi
                ! 
                iq=iqi
                !
                ikq  = index_xkq(current_k,iq)
                ik   = index_xk(ikq)
                isym = abs(index_sym(ikq))

                ! FIXME: use cryst_to_cart and company as above..
                xk_cryst(:)=at(1,:)*xk(1,ik)+at(2,:)*xk(2,ik)+at(3,:)*xk(3,ik)
                IF (index_sym(ikq) < 0) xk_cryst = -xk_cryst
                sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                         s(:,2,isym)*xk_cryst(2) + &
                         s(:,3,isym)*xk_cryst(3) 
                xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

                !CALL start_clock ('exxen2_ngmloop')

!$omp parallel do default(shared), private(ig, beta, q, qq, on_double_grid, x)
                DO ig = 1, ngm
                  q(1)= xk(1,current_k) - xkq(1) + g(1,ig)
                  q(2)= xk(2,current_k) - xkq(2) + g(2,ig)
                  q(3)= xk(3,current_k) - xkq(3) + g(3,ig)

                  q = q * tpiba
                  qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) )

                  DO beta = 1, 3
                      fac_tens(1:3,beta,ig) = q(1:3)*q(beta)
                  ENDDO

                  IF (x_gamma_extrapolation) THEN
                      on_double_grid = .true.
                      x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                      x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                      x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                  ELSE
                      on_double_grid = .FALSE.
                  ENDIF

                  IF (use_coulomb_vcut_ws) THEN
                      fac(ig) = vcut_get(vcut, q)
                      fac_stress(ig) = 0._dp   ! not implemented
                      IF (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)

                  ELSE IF ( use_coulomb_vcut_spheric ) THEN
                      fac(ig) = vcut_spheric_get(vcut, q)
                      fac_stress(ig) = 0._dp   ! not implemented
                      IF (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig) 

                  ELSE IF (gau_scrlen > 0) then
                      fac(ig)=e2*((pi/gau_scrlen)**(1.5d0))* &
                            exp(-qq/4.d0/gau_scrlen) * grid_factor
                      fac_stress(ig) =  e2*2.d0/4.d0/gau_scrlen * &
                            exp(-qq/4.d0/gau_scrlen) *((pi/gau_scrlen)**(1.5d0))* &
                                                                     grid_factor
                      IF (gamma_only) fac(ig) = 2.d0 * fac(ig)
                      IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)
                      IF (on_double_grid) fac(ig) = 0._dp
                      IF (on_double_grid) fac_stress(ig) = 0._dp

                  ELSE IF (qq > 1.d-8) THEN
                      IF ( erfc_scrlen > 0 ) THEN
                        fac(ig)=e2*fpi/qq*(1._dp-exp(-qq/4.d0/erfc_scrlen**2)) * grid_factor
                        fac_stress(ig) = -e2*fpi * 2.d0/qq**2 * ( &
                            (1._dp+qq/4.d0/erfc_scrlen**2)*exp(-qq/4.d0/erfc_scrlen**2) - 1._dp) * &
                            grid_factor
                      ELSE
                        fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor
                        fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2 * grid_factor
                      ENDIF

                      IF (gamma_only) fac(ig) = 2.d0 * fac(ig)
                      IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)
                      IF (on_double_grid) fac(ig) = 0._dp
                      IF (on_double_grid) fac_stress(ig) = 0._dp

                  ELSE
                      fac(ig)= -exxdiv ! or rather something else (see f.gygi)
                      fac_stress(ig) = 0._dp  ! or -exxdiv_stress (not yet implemented)
                      IF ( yukawa> 0._dp .and. .not. x_gamma_extrapolation) THEN
                        fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
                        fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2
                      ENDIF
                      IF (erfc_scrlen > 0._dp .and. .not. x_gamma_extrapolation) THEN
                        fac(ig) = e2*fpi / (4.d0*erfc_scrlen**2)
                        fac_stress(ig) = e2*fpi / (8.d0*erfc_scrlen**4)
                      ENDIF
                  ENDIF
                ENDDO
!$omp end parallel do
                !CALL stop_clock ('exxen2_ngmloop')

                IF (gamma_only) THEN
                    !
                    h_ibnd = ibnd_start/2
                    !
                    IF(MOD(ibnd_start,2)==0) THEN
                      h_ibnd=h_ibnd-1
                      ibnd_loop_start=ibnd_start-1
                    ELSE
                      ibnd_loop_start=ibnd_start
                    ENDIF
                    !
                    DO ibnd = ibnd_loop_start, ibnd_end, 2     !for each band of psi
                        !
                        h_ibnd = h_ibnd + 1
                        !
                        IF( ibnd < ibnd_start ) THEN
                            x1 = 0._dp
                        ELSE
                            x1 = x_occupation(ibnd,  ik)
                        ENDIF

                        IF( ibnd == ibnd_end) THEN
                            x2 = 0._dp
                        ELSE
                            x2 = x_occupation(ibnd+1,  ik)
                        ENDIF
                        IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE
                        !
                        ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                        DO ir = 1, nrxxs
                            tempphic(ir) = exxbuff(ir,h_ibnd,ikq)
                            rhoc(ir)     = CONJG(tempphic(ir))*temppsic(ir) / omega
                        ENDDO
!$omp end parallel do
                        ! bring it to G-space
                        !<<<
                        !CALL fwfft ('Custom', rhoc, exx_fft%dfftt)
                        CALL fwfft ('Custom', rhoc, exx_fft%dfftt, is_exx=.TRUE.)
                        !>>>
    
                        vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                        DO ig = 1, ngm
                            !
                            vc(:,:) = vc(:,:) + fac(ig) * x1 * &
                                      abs( rhoc(exx_fft%nlt(ig)) + &
                                      CONJG(rhoc(exx_fft%nltm(ig))))**2 * &
                                      (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                            vc(:,:) = vc(:,:) + fac(ig) * x2 * &
                                      abs( rhoc(exx_fft%nlt(ig)) - &
                                      CONJG(rhoc(exx_fft%nltm(ig))))**2 * &
                                      (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                        enddo
!$omp end parallel do
                        vc = vc / nqs / 4.d0
! Debug Zhenfei Liu 9/26/2015
                        exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
!                        exx_stress_ = exx_stress_ + vc * wg(jbnd,ikk)
! DONE Debug.
                    ENDDO

                ELSE

                    DO ibnd = ibnd_start, ibnd_end    !for each band of psi
                      !
                      IF ( abs(x_occupation(ibnd,ik)) < 1.d-6) cycle
                      !
                      ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                          tempphic(ir) = exxbuff(ir,ibnd,ikq)
                          rhoc(ir)     = CONJG(tempphic(ir))*temppsic(ir) / omega
                      ENDDO
!$omp end parallel do

                      ! bring it to G-space
                      !<<<
                      !CALL fwfft ('Custom', rhoc, exx_fft%dfftt)
                      CALL fwfft ('Custom', rhoc, exx_fft%dfftt, is_exx=.TRUE.)
                      !>>>

                      vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                      DO ig = 1, ngm
                          vc(:,:) = vc(:,:) + rhoc(exx_fft%nlt(ig))  * &
                                        CONJG(rhoc(exx_fft%nlt(ig)))* &
                                    (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                      ENDDO
!$omp end parallel do
                      vc = vc * x_occupation(ibnd,ik) / nqs / 4.d0
! Debug Zhenfei Liu 9/26/2015
                      exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
!                      exx_stress_ = exx_stress_ + vc * wg(jbnd,ikk)
! DONE Debug.

                    ENDDO

                ENDIF ! gamma or k-points

            ENDDO ! iqi
        ENDDO ! jbnd
    ENDDO ! ikk

    DEALLOCATE(tempphic, temppsic, rhoc, fac, fac_tens, fac_stress )
    !
    CALL mp_sum( exx_stress_, intra_egrp_comm )
    CALL mp_sum( exx_stress_, inter_egrp_comm )
    CALL mp_sum( exx_stress_, inter_pool_comm )
    exx_stress = exx_stress_

    CALL stop_clock ('exx_stress')
    !-----------------------------------------------------------------------
  END FUNCTION exx_stress
  !-----------------------------------------------------------------------











  !-----------------------------------------------------------------------
  SUBROUTINE convert_evc(ik_in)
  !-----------------------------------------------------------------------
    ! ... generic, k-point version of vexx
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega, at, bg, tpiba2
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, npw, igk, current_k, ecutwfc, g2kin, nbnd
    USE control_flags,  ONLY : gamma_only
    USE klist,          ONLY : xk, nks, nkstot, ngk
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_exx,       ONLY : inter_egrp_comm, intra_egrp_comm, my_egrp_id, &
                               negrp, nproc_egrp, me_egrp
    USE gvect,              ONLY : ig_l2g, mill_g
    USE wavefunctions_module, ONLY : evc
    USE mp,             ONLY : mp_sum, mp_barrier, mp_bcast, mp_size, mp_rank
    USE uspp,           ONLY : nkb, okvan, vkb
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE recvec_subs,    ONLY : ggen 
    USE fft_base,             ONLY : dfftp, dffts
    USE mp_pools
    USE io_files,             ONLY : nwordwfc, iunwfc, iunigk, &
                                     iunigk_exx, iunwfc_exx
    USE buffers,              ONLY : open_buffer, get_buffer, save_buffer, &
         close_buffer
    USE control_flags,        ONLY : io_level
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik_in
    INTEGER :: lda
    INTEGER :: n, m
    !COMPLEX(DP) :: psi(lda*npol,m) 
    !COMPLEX(DP) :: hpsi(lda*npol,m)
    
    COMPLEX(DP), ALLOCATABLE :: psi_sum(:), psi_temp(:)

    !INTEGER, ALLOCATABLE :: 
    INTEGER :: i, j, k, im, ig, ik
    INTEGER :: total_lda, prev_lda, total_lda_exx, prev_lda_exx
    INTEGER :: lda_local_evc

    LOGICAL :: exst_mem, exst_file

    CALL start_clock ('conv_evc')


    IF (negrp.eq.1) THEN

       !get evc_exx
       IF(.not.allocated(evc_exx))THEN
          ALLOCATE(evc_exx(npwx,nbnd))
       END IF
       evc_exx = evc

       !get vkb_exx
       IF(.not.allocated(vkb_exx)) THEN
          ALLOCATE( vkb_exx( npwx, nkb ) )
          vkb_exx = vkb
       END IF

       !get igk_exx
       IF(.not.allocated(igk_exx)) THEN
          ALLOCATE( igk_exx( npwx ) )
          igk_exx = igk
       END IF

       !this is to ensure that the correct buffer is used
       iunwfc_exx = iunwfc
       iunigk_exx = iunigk
       nwordwfc_exx = nwordwfc

       CALL stop_clock ('conv_evc')

       RETURN

    END IF


    ! NOTE: POSSIBLE SOURCE OF ERRORS HERE
    lda = npwx
    n = npwx 
    npwx_local = npwx
    !npw_local = npw
    IF( .not.allocated(ngk_local) ) allocate(ngk_local(nks))
    ngk_local = ngk

    IF ( .not.allocated(comm_recv) ) THEN
       !initialize all of the conversion maps
       CALL initialize_local_to_exact_map(lda, nbnd)
    ELSE
       !just change the data structure
       CALL change_data_structure(.TRUE.)
    END IF


    ! NOTE: POSSIBLE SOURCE OF ERRORS HERE
    lda = npwx
    n = npwx 
    npwx_exx = npwx
    !npw_exx = npw
    IF( .not.allocated(ngk_exx) ) allocate(ngk_exx(nks))
    ngk_exx = ngk
    
    !get evc_exx
    IF(.not.allocated(evc_exx))THEN
       ALLOCATE(evc_exx(lda,nbnd))
       !<<<
  !
  ! ... open files/buffer for wavefunctions (nwordwfc set in openfil)
  ! ... io_level > 1 : open file, otherwise: open buffer
  !
       nwordwfc_exx  = nbnd*npwx_exx*npol
       CALL open_buffer( iunwfc_exx, 'wfc_exx', nwordwfc_exx, io_level, &
            exst_mem, exst_file )
       !>>>
!    ELSE
!       CALL close_buffer ( iunwfc_exx, 'DELETE' )
!       CALL open_buffer( iunwfc_exx, 'wfc_exx', nwordwfc_exx, io_level, &
!            exst_mem, exst_file )
    END IF
!    CALL reconstruct_for_exact(lda, n, nbnd, ik, evc, evc_exx, 1)

    !get vkb_exx
    IF(.not.allocated(vkb_exx)) THEN
       ALLOCATE( vkb_exx( npwx, nkb ) )
       vkb_exx = 0.0_DP

!       DO ik=1, nks
          CALL reconstruct_for_exact(lda, n, nkb, 1, vkb, vkb_exx, 2)
!       END DO
       
    END IF

!    IF ( nks > 1 ) REWIND( iunigk )
    DO ik=1, nks
!       current_k = ik

       IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       CALL reconstruct_for_exact(lda, n, nbnd, ik, evc, evc_exx, 1)
       IF ( nks > 1 ) CALL save_buffer ( evc_exx, nwordwfc_exx, iunwfc_exx, ik )
    END DO

    !<<<
    !IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, current_k)
    !>>>
    



    CALL stop_clock ('conv_evc')

  END SUBROUTINE convert_evc








  !-----------------------------------------------------------------------
  SUBROUTINE start_exx_parallelization(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    USE becmod,       ONLY : bec_type
    USE mp_exx,       ONLY : negrp
    USE wvfct,        ONLY : current_k, npwx
    !
    !
    IMPLICIT NONE
    !
    Integer :: lda
    INTEGER :: n, m
    COMPLEX(DP) :: psi(lda*npol,m) 
    COMPLEX(DP) :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL   :: becpsi ! or call a calbec(...psi) instead
    
    INTEGER :: ig
    
    CALL start_clock ('start_exxp')

    npwx_local = npwx
    IF ( .not.allocated(comm_recv) ) THEN
       !initialize all of the conversion maps
       CALL initialize_local_to_exact_map(lda, m)
    ELSE
       !just change the data structure
       CALL change_data_structure(.TRUE.)
    END IF
    npwx_exx = npwx
    
    !!!!
    !get igk
    CALL update_igk(.TRUE.)
    !!!!

    !get psi_exx
    CALL reconstruct_for_exact(lda, n, m, current_k, psi, psi_exx, 0)

    !get hpsi_exx
!    CALL reconstruct_for_exact(lda, n, m, current_k, hpsi, hpsi_exx, 0)
    hpsi_exx = 0.d0

    CALL stop_clock ('start_exxp')

  END SUBROUTINE start_exx_parallelization







  !-----------------------------------------------------------------------
  SUBROUTINE end_exx_parallelization(lda, n, m, psi, hpsi)
  !-----------------------------------------------------------------------
    USE mp_exx,       ONLY : negrp
    !
    IMPLICIT NONE
    !
    INTEGER :: lda
    INTEGER :: n, m
    COMPLEX(DP) :: psi(lda_original*npol,m)
    COMPLEX(DP) :: hpsi(lda_original*npol,m)
 
    INTEGER :: ig
    COMPLEX(DP) :: hpsi_temp(lda_original*npol,m)

    CALL start_clock ('end_exxp')

    CALL change_data_structure(.FALSE.)

    !!!!
    !get igk
    CALL update_igk(.FALSE.)
    !!!!
    
    !CALL deconstruct_for_exact(m,hpsi_exx,hpsi)
    !CALL deconstruct_for_exact(m,hpsi_exx,hpsi_temp)
    !hpsi = hpsi + hpsi_temp
    CALL deconstruct_for_exact(m,hpsi_exx,hpsi)

    CALL stop_clock ('end_exxp')

  END SUBROUTINE end_exx_parallelization












  !-----------------------------------------------------------------------
  SUBROUTINE initialize_local_to_exact_map(lda, m)
  !-----------------------------------------------------------------------
    ! ... generic, k-point version of vexx
    !
    USE io_files,       ONLY : iunigk, iunigk_exx
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega, at, bg, tpiba2
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, npw, igk, current_k, ecutwfc, g2kin, nbnd
    USE control_flags,  ONLY : gamma_only
    USE klist,          ONLY : xk, nks, nkstot, ngk
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_exx,       ONLY : inter_egrp_comm, intra_egrp_comm, my_egrp_id, &
                               negrp, nproc_egrp, me_egrp
    USE gvect,              ONLY : ig_l2g, mill_g
    USE mp,             ONLY : mp_sum, mp_barrier, mp_bcast, mp_size, mp_rank
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE recvec_subs,        ONLY : ggen 
    USE fft_base,             ONLY : dfftp, dffts
    USE mp_pools
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_STATUS_SIZE, MPI_DOUBLE_COMPLEX
#endif
    !
    !
    IMPLICIT NONE
    !
    Integer :: lda
    INTEGER :: n, m
    
    INTEGER, ALLOCATABLE :: local_map(:,:), exx_map(:,:)
    INTEGER, ALLOCATABLE :: l2e_map(:,:), e2l_map(:,:)
    INTEGER, ALLOCATABLE :: psi_source(:), psi_source_exx(:)
    INTEGER :: current_index

    INTEGER :: i, j, k, ik, im, ig, count, iproc, prev
    INTEGER :: total_lda(nks), prev_lda(nks)
    INTEGER :: total_lda_exx(nks), prev_lda_exx(nks)
    INTEGER :: n_local
    INTEGER :: lda_max_local, lda_max_exx

    INTEGER :: request_send(nproc_egrp), request_recv(nproc_egrp)
    INTEGER :: ierr
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
    INTEGER :: egrp_base, total_lda_egrp(nks), prev_lda_egrp(nks)
    INTEGER :: igk_loc(npwx)
    
    CALL start_clock ('init_exxp')

    IF ( .not.allocated(comm_recv) ) THEN
       ALLOCATE(comm_recv(nproc_egrp,nks),comm_send(nproc_egrp,nks))
    END IF
    IF ( .not.allocated(lda_local) ) THEN
       ALLOCATE(lda_local(nproc_pool,nks))
       ALLOCATE(lda_exx(nproc_egrp,nks))
    END IF

    ! store the original values of lda and n
    lda_original = lda
    n_original = n

    !construct the local map
    !allocate(all_lda(nproc_pool,nks))
    !all_lda = 0
    lda_local = 0
    IF ( nks > 1 ) REWIND( iunigk )
    IF ( nks.eq.1 ) igk_loc = igk
    DO ik = 1, nks
       IF ( nks > 1 ) READ( iunigk ) igk_loc

       n = 0
       DO i = 1, size(igk_loc)
          IF(igk_loc(i).gt.0)n = n + 1
       END DO
       lda_local(me_pool+1,ik) = n
       CALL mp_sum(lda_local(:,ik),intra_pool_comm)
       total_lda(ik) = sum(lda_local(:,ik))
       prev_lda(ik) = sum(lda_local(1:me_pool,ik))
    END DO

    ALLOCATE(local_map(maxval(total_lda),nks))
    local_map = 0
    IF ( nks > 1 ) REWIND( iunigk )
    DO ik = 1, nks
       IF ( nks > 1 ) READ( iunigk ) igk_loc
    
       local_map(prev_lda(ik)+1:prev_lda(ik)+lda_local(me_pool+1,ik),ik) = &
            ig_l2g(igk_loc(1:lda_local(me_pool+1,ik)))

    END DO
    CALL mp_sum(local_map,intra_pool_comm)

    !-----------------------------------------!
    ! Switch to the exx data structure        !
    !-----------------------------------------!
    CALL change_data_structure(.TRUE.)
    ! NOTE: somehow, lda gets changed by change_data_structure



    !construct the exx map
    lda_exx = 0
    IF ( nks > 1 ) REWIND( iunigk_exx )
    DO ik = 1, nks
       IF ( nks > 1 ) READ( iunigk_exx ) igk_exx

       n = 0
       DO i = 1, size(igk_exx)
          IF(igk_exx(i).gt.0)n = n + 1
       END DO
       lda_exx(me_egrp+1,ik) = n
       CALL mp_sum(lda_exx(:,ik),intra_egrp_comm)
       total_lda_exx(ik) = sum(lda_exx(:,ik))
       prev_lda_exx(ik) = sum(lda_exx(1:me_egrp,ik))

    END DO

    ALLOCATE(exx_map(maxval(total_lda_exx),nks))
    exx_map = 0
    IF ( nks > 1 ) REWIND( iunigk_exx )
    DO ik = 1, nks
       IF ( nks > 1 ) READ( iunigk_exx ) igk_exx

       exx_map(prev_lda_exx(ik)+1:prev_lda_exx(ik)+lda_exx(me_egrp+1,ik),ik) = &
            ig_l2g(igk_exx(1:lda_exx(me_egrp+1,ik)))
    
    END DO
    CALL mp_sum(exx_map,intra_egrp_comm)

    !construct the l2e_map
    !NOTE: For multiple k-points, need a different l2e_map for each k-point
    allocate( l2e_map(maxval(total_lda_exx),nks) )
    l2e_map = 0
    DO ik = 1, nks
       DO ig = 1, lda_exx(me_egrp+1,ik)
          DO j=1, total_lda(ik)
             IF( local_map(j,ik).EQ.exx_map(ig+prev_lda_exx(ik),ik) ) exit
          END DO
          l2e_map(ig+prev_lda_exx(ik),ik) = j
       END DO
    END DO
    CALL mp_sum(l2e_map,intra_egrp_comm)







    lda_max_local = maxval(lda_local)
    lda_max_exx = maxval(lda_exx)
    allocate(psi_source(maxval(total_lda_exx)))

    DO ik = 1, nks !!!!!!

    !determine where each value is coming from
    psi_source = 0
    DO ig = 1, lda_exx(me_egrp+1,ik)
       j = 1
       DO i = 1, nproc_pool
          j = j + lda_local(i,ik)
          IF( j.gt.l2e_map(ig+prev_lda_exx(ik),ik) ) exit
       END DO
       psi_source(ig+prev_lda_exx(ik)) = i-1
    END DO
    CALL mp_sum(psi_source,intra_egrp_comm)

    !allocate communication packets to recieve psi and hpsi
    DO iproc=0, nproc_egrp-1
       
       !determine how many values need to come from the target
       count = 0
       DO ig=1, lda_exx(me_egrp+1,ik)
          IF ( MODULO(psi_source(ig+prev_lda_exx(ik)),nproc_egrp).eq.iproc ) THEN
             count = count + 1
          END IF
       END DO

       !allocate the communication packet
       comm_recv(iproc+1,ik)%size = count
       IF (count.gt.0) THEN
          IF (.not.ALLOCATED(comm_recv(iproc+1,ik)%msg)) THEN
             ALLOCATE(comm_recv(iproc+1,ik)%indices(count))
             ALLOCATE(comm_recv(iproc+1,ik)%msg(count,m))
             ALLOCATE(comm_recv(iproc+1,ik)%msg_evc(count,nbnd))
             ALLOCATE(comm_recv(iproc+1,ik)%msg_vkb(count,nkb))
          END IF
       END IF

       !determine which values need to come from the target
       count = 0
       DO ig=1, lda_exx(me_egrp+1,ik)
          IF ( MODULO(psi_source(ig+prev_lda_exx(ik)),nproc_egrp).eq.iproc ) THEN
             count = count + 1
             comm_recv(iproc+1,ik)%indices(count) = ig
          END IF
       END DO

    END DO

    !allocate communication packets to send psi and hpsi
    prev = 0
    DO iproc=0, nproc_egrp-1
       
       !determine how many values need to be sent to the target
       count = 0
       DO ig=1, lda_exx(iproc+1,ik)
          IF ( MODULO(psi_source(ig+prev),nproc_egrp).eq.me_egrp ) THEN
             count = count + 1
          END IF
       END DO

       !allocate the communication packet
       comm_send(iproc+1,ik)%size = count
       IF (count.gt.0) THEN
          IF (.not.ALLOCATED(comm_send(iproc+1,ik)%msg)) THEN
             ALLOCATE(comm_send(iproc+1,ik)%indices(count))
             ALLOCATE(comm_send(iproc+1,ik)%msg(count,m))
             ALLOCATE(comm_send(iproc+1,ik)%msg_evc(count,nbnd))
             ALLOCATE(comm_send(iproc+1,ik)%msg_vkb(count,nkb))
          END IF
       END IF
          
       !determine which values need to be sent to the target
       count = 0
       DO ig=1, lda_exx(iproc+1,ik)
          IF ( MODULO(psi_source(ig+prev),nproc_egrp).eq.me_egrp ) THEN
             count = count + 1
             comm_send(iproc+1,ik)%indices(count) = l2e_map(ig+prev,ik)
          END IF
       END DO

       prev = prev + lda_exx(iproc+1,ik)
    END DO

    END DO

    !allocate psi_exx and hpsi_exx
    !!!!!!!!!POSSIBLE ERROR: PSI_EXX REALLY NEEDS LDA (NOT N)
    IF(allocated(psi_exx))DEALLOCATE(psi_exx)
    ALLOCATE(psi_exx(npwx*npol,m))
    IF(allocated(hpsi_exx))DEALLOCATE(hpsi_exx)
    ALLOCATE(hpsi_exx(npwx*npol,m))





















    IF ( .not.allocated(comm_recv_reverse) ) THEN
       ALLOCATE(comm_recv_reverse(nproc_egrp,nks))
       ALLOCATE(comm_send_reverse(nproc_egrp,nks))
    END IF

    egrp_base = my_egrp_id*nproc_egrp
    DO ik = 1, nks
       total_lda_egrp(ik) = &
            sum( lda_local(egrp_base+1:(egrp_base+nproc_egrp),ik) )
       prev_lda_egrp(ik) = &
            sum( lda_local(egrp_base+1:(egrp_base+me_egrp),ik) )
    END DO

    allocate( e2l_map(maxval(total_lda_egrp),nks) )
    e2l_map = 0
    DO ik = 1, nks
       DO ig = 1, lda_local(me_pool+1,ik)
          DO j=1, total_lda_exx(ik)
             IF( local_map(ig+prev_lda(ik),ik).EQ.exx_map(j,ik) ) exit
          END DO
          e2l_map(ig+prev_lda_egrp(ik),ik) = j
       END DO
    END DO
    CALL mp_sum(e2l_map,intra_egrp_comm)

    allocate(psi_source_exx( maxval(total_lda_egrp) ))




    DO ik = 1, nks !!!!!!

    !determine where each value is coming from
    psi_source_exx = 0
    DO ig = 1, lda_local(me_pool+1,ik)
       j = 1
       DO i = 1, nproc_egrp
          j = j + lda_exx(i,ik)
          IF( j.gt.e2l_map(ig+prev_lda_egrp(ik),ik) ) exit
       END DO
       psi_source_exx(ig+prev_lda_egrp(ik)) = i-1
    END DO
    CALL mp_sum(psi_source_exx,intra_egrp_comm)

    !allocate communication packets to recieve psi and hpsi (reverse)
    DO iproc=0, nproc_egrp-1
       
       !determine how many values need to come from the target
       count = 0
       DO ig=1, lda_local(me_pool+1,ik)
          IF ( psi_source_exx(ig+prev_lda_egrp(ik)).eq.iproc ) THEN
             count = count + 1
          END IF
       END DO

       !allocate the communication packet
       comm_recv_reverse(iproc+1,ik)%size = count
       IF (count.gt.0) THEN
          IF (.not.ALLOCATED(comm_recv_reverse(iproc+1,ik)%msg)) THEN
             ALLOCATE(comm_recv_reverse(iproc+1,ik)%indices(count))
             ALLOCATE(comm_recv_reverse(iproc+1,ik)%msg(count,m))
          END IF
       END IF

       !determine which values need to come from the target
       count = 0
       DO ig=1, lda_local(me_pool+1,ik)
          IF ( psi_source_exx(ig+prev_lda_egrp(ik)).eq.iproc ) THEN
             count = count + 1
             comm_recv_reverse(iproc+1,ik)%indices(count) = ig
          END IF
       END DO

    END DO

    !allocate communication packets to send psi and hpsi
    prev = 0
    DO iproc=0, nproc_egrp-1

       !determine how many values need to be sent to the target
       count = 0
       DO ig=1, lda_local(iproc+egrp_base+1,ik)
          IF ( psi_source_exx(ig+prev).eq.me_egrp ) THEN
             count = count + 1
          END IF
       END DO

       !allocate the communication packet
       comm_send_reverse(iproc+1,ik)%size = count
       IF (count.gt.0) THEN
          IF (.not.ALLOCATED(comm_send_reverse(iproc+1,ik)%msg)) THEN
             ALLOCATE(comm_send_reverse(iproc+1,ik)%indices(count))
             ALLOCATE(comm_send_reverse(iproc+1,ik)%msg(count,m))
          END IF
       END IF
          
       !determine which values need to be sent to the target
       count = 0
       DO ig=1, lda_local(iproc+egrp_base+1,ik)
          IF ( psi_source_exx(ig+prev).eq.me_egrp ) THEN
             count = count + 1
             comm_send_reverse(iproc+1,ik)%indices(count) = e2l_map(ig+prev,ik)
          END IF
       END DO

       prev = prev + lda_local( egrp_base+iproc+1, ik )
    END DO

    END DO

    !deallocate arrays
    DEALLOCATE( local_map, exx_map )
    DEALLOCATE( l2e_map, e2l_map )
    DEALLOCATE( psi_source, psi_source_exx )

    CALL stop_clock ('init_exxp')

  END SUBROUTINE initialize_local_to_exact_map












  !-----------------------------------------------------------------------
  SUBROUTINE reconstruct_for_exact(lda, n, m, ik, psi, psi_out, type)
  !-----------------------------------------------------------------------
    USE mp,           ONLY : mp_sum
    USE mp_pools,     ONLY : nproc_pool, me_pool
    USE mp_exx,       ONLY : intra_egrp_comm, inter_egrp_comm, &
         nproc_egrp, me_egrp, negrp, my_egrp_id, nibands, ibands
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_STATUS_SIZE, MPI_DOUBLE_COMPLEX
#endif
    USE klist,        ONLY : xk, wk, nkstot, nks, qnorm
    !
    !
    IMPLICIT NONE
    !
    Integer :: lda
    INTEGER :: n, m
    !!!!!!!POSIBLE ERROR HERE: PSI REALLY NEEDS LDA FOR SIZE (NOT N)
    COMPLEX(DP) :: psi(npwx_local*npol,m) 
    COMPLEX(DP) :: psi_out(npwx_exx*npol,m)
!    COMPLEX(DP) :: psi(lda_local(me_pool+1,ik)*npol,m) 
!    COMPLEX(DP) :: psi_out(lda*npol,m)
    INTEGER, INTENT(in) :: type

    COMPLEX(DP), ALLOCATABLE :: psi_work(:,:,:), psi_gather(:,:)
    INTEGER :: i, j, im, iproc, ig, ik, iegrp
    INTEGER :: prev, lda_max_local

    INTEGER :: request_send(nproc_egrp), request_recv(nproc_egrp)
    INTEGER :: ierr
    INTEGER :: current_ik
    INTEGER, EXTERNAL :: find_current_k
    INTEGER :: recvcount(negrp)
    INTEGER :: displs(negrp)
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
#endif

    lda_max_local = maxval(lda_local)

    !current_ik=find_current_k(ik, nkstot, nks)
    current_ik = ik

    !-------------------------------------------------------!
    !Communication Part 1
    !-------------------------------------------------------!
    CALL start_clock ('comm1')
    allocate(psi_work(lda_max_local,m,negrp))
!    allocate(psi_work(npwx_local,m,negrp))
    !<<<
!    psi_work = 0.0_DP
!    DO im=1, m
!       !gather psi and hpsi onto each exact group
!       psi_work(1:lda_local(me_pool+1,current_ik), im, my_egrp_id+1) = &
!            psi(1:lda_local(me_pool+1,current_ik), im)
!    END DO
    !NOTE: Need to change this to an allgather
!    call mp_sum(psi_work,inter_egrp_comm)
    allocate(psi_gather(lda_max_local,m))
    DO im=1, m
       psi_gather(1:lda_local(me_pool+1,current_ik),im) = psi(:,im)
    END DO
    CALL stop_clock ('comm1')
    CALL start_clock ('comm2')
    IF ( type.eq.0 ) THEN

       recvcount = lda_max_local
       DO iegrp=1, negrp
          displs(iegrp) = (iegrp-1)*(lda_max_local*m)
       END DO
       DO iegrp=1, negrp
          !CALL MPI_GATHER( psi_gather, &
          !     lda_max_local*m, MPI_DOUBLE_COMPLEX, &
          !     psi_work, &
          !     lda_max_local*m, MPI_DOUBLE_COMPLEX, &
          !     iegrp-1, &
          !     inter_egrp_comm, ierr )
          DO im=1, nibands(iegrp)
             IF ( my_egrp_id.eq.(iegrp-1) ) THEN
                DO j=1, negrp
                   displs(j) = (j-1)*(lda_max_local*m) + &
                        lda_max_local*(ibands(im,iegrp)-1)
                END DO
             END IF
             CALL MPI_GATHERV( psi_gather(:, ibands(im,iegrp) ), &
                  lda_max_local, MPI_DOUBLE_COMPLEX, &
                  psi_work, &
                  recvcount, displs, MPI_DOUBLE_COMPLEX, &
                  iegrp-1, &
                  inter_egrp_comm, ierr )
          END DO
       END DO

    ELSE
       
       CALL MPI_ALLGATHER( psi_gather, &
            lda_max_local*m, MPI_DOUBLE_COMPLEX, &
            psi_work, &
            lda_max_local*m, MPI_DOUBLE_COMPLEX, &
            inter_egrp_comm, ierr )
       
    END IF
    CALL stop_clock ('comm2')
    CALL start_clock ('comm3')
    !>>>
    
    !-------------------------------------------------------!
    !Communication Part 2
    !-------------------------------------------------------!
    !send communication packets
    DO iproc=0, nproc_egrp-1
       IF ( comm_send(iproc+1,current_ik)%size.gt.0) THEN
          DO i=1, comm_send(iproc+1,current_ik)%size
             ig = comm_send(iproc+1,current_ik)%indices(i)

             !determine which egrp this corresponds to
             prev = 0
             DO j=1, nproc_pool
                IF ((prev+lda_local(j,current_ik)).ge.ig) THEN 
                   ig = ig - prev
                   exit
                END IF
                prev = prev + lda_local(j,current_ik)
             END DO

             IF ( type.eq.0 ) THEN !psi or hpsi
                !DO im=1, m
                !   comm_send(iproc+1,current_ik)%msg(i,im) = &
                !        psi_work(ig,im,1+(j-1)/nproc_egrp)
                !END DO
                DO im=1, nibands(my_egrp_id+1)
                   comm_send(iproc+1,current_ik)%msg(i,im) = &
                        psi_work(ig,ibands(im,my_egrp_id+1),1+(j-1)/nproc_egrp)
                END DO
             ELSE IF (type.eq.1) THEN !evc
                DO im=1, m
                   comm_send(iproc+1,current_ik)%msg_evc(i,im) = &
                        psi_work(ig,im,1+(j-1)/nproc_egrp)
                END DO
             ELSE IF (type.eq.2) THEN !vkb
                DO im=1, m
                   comm_send(iproc+1,current_ik)%msg_vkb(i,im) = &
                        psi_work(ig,im,1+(j-1)/nproc_egrp)
                END DO
             END IF

          END DO


          !send the message
          IF ( type.eq.0 ) THEN !psi or hpsi
             CALL MPI_ISEND( comm_send(iproc+1,current_ik)%msg, &
                  comm_send(iproc+1,current_ik)%size*nibands(my_egrp_id+1), &
                  MPI_DOUBLE_COMPLEX, &
                  iproc, 100+iproc*nproc_egrp+me_egrp, &
                  intra_egrp_comm, request_send(iproc+1), ierr )
          ELSE IF (type.eq.1) THEN !evc
             CALL MPI_ISEND( comm_send(iproc+1,current_ik)%msg_evc, &
                  comm_send(iproc+1,current_ik)%size*m, MPI_DOUBLE_COMPLEX, &
                  iproc, 100+iproc*nproc_egrp+me_egrp, &
                  intra_egrp_comm, request_send(iproc+1), ierr )
          ELSE IF (type.eq.2) THEN ! vkb
             CALL MPI_ISEND( comm_send(iproc+1,current_ik)%msg_vkb, &
                  comm_send(iproc+1,current_ik)%size*m, MPI_DOUBLE_COMPLEX, &
                  iproc, 100+iproc*nproc_egrp+me_egrp, &
                  intra_egrp_comm, request_send(iproc+1), ierr )
          END IF

          
       END IF
    END DO
    
    !begin recieving the communication packets
    DO iproc=0, nproc_egrp-1
       IF ( comm_recv(iproc+1,current_ik)%size.gt.0) THEN

          !recieve the message
          IF (type.eq.0) THEN !psi or hpsi
             CALL MPI_IRECV( comm_recv(iproc+1,current_ik)%msg, &
                  comm_recv(iproc+1,current_ik)%size*nibands(my_egrp_id+1), &
                  MPI_DOUBLE_COMPLEX, &
                  iproc, 100+me_egrp*nproc_egrp+iproc, &
                  intra_egrp_comm, request_recv(iproc+1), ierr )
          ELSE IF (type.eq.1) THEN !evc
             CALL MPI_IRECV( comm_recv(iproc+1,current_ik)%msg_evc, &
                  comm_recv(iproc+1,current_ik)%size*m, MPI_DOUBLE_COMPLEX, &
                  iproc, 100+me_egrp*nproc_egrp+iproc, &
                  intra_egrp_comm, request_recv(iproc+1), ierr )
          ELSE IF (type.eq.2) THEN !vkb
             CALL MPI_IRECV( comm_recv(iproc+1,current_ik)%msg_vkb, &
                  comm_recv(iproc+1,current_ik)%size*m, MPI_DOUBLE_COMPLEX, &
                  iproc, 100+me_egrp*nproc_egrp+iproc, &
                  intra_egrp_comm, request_recv(iproc+1), ierr )
          END IF
          
       END IF
    END DO

    !assign psi
    DO iproc=0, nproc_egrp-1
       IF ( comm_recv(iproc+1,current_ik)%size.gt.0 ) THEN
          
          CALL MPI_WAIT(request_recv(iproc+1), istatus, ierr)

          DO i=1, comm_recv(iproc+1,current_ik)%size
             ig = comm_recv(iproc+1,current_ik)%indices(i)

             !SET PSI_EXX HERE
             IF (type.eq.0) THEN !psi or hpsi
                !DO im=1, m
                !   psi_out(ig,im) = comm_recv(iproc+1,current_ik)%msg(i,im)
                !END DO
                DO im=1, nibands(my_egrp_id+1)
                   psi_out(ig,ibands(im,my_egrp_id+1)) = &
                        comm_recv(iproc+1,current_ik)%msg(i,im)
                        !comm_recv(iproc+1,current_ik)%msg(i,ibands(im,my_egrp_id+1))
                END DO
             ELSE IF (type.eq.1) THEN !evc
                DO im=1, m
                   psi_out(ig,im) = comm_recv(iproc+1,current_ik)%msg_evc(i,im)
                END DO
             ELSE IF (type.eq.2) THEN !vkb
                DO im=1, m
                   psi_out(ig,im) = comm_recv(iproc+1,current_ik)%msg_vkb(i,im)
                END DO
             END IF

          END DO

       END IF
    END DO

    !now wait for everything to finish sending
    DO iproc=0, nproc_egrp-1
       IF ( comm_send(iproc+1,current_ik)%size.gt.0 ) THEN
          CALL MPI_WAIT(request_send(iproc+1), istatus, ierr)
       END IF
    END DO

    !deallocate arrays
    DEALLOCATE( psi_work, psi_gather )

    CALL stop_clock ('comm3')

  END SUBROUTINE reconstruct_for_exact











  !-----------------------------------------------------------------------
  SUBROUTINE deconstruct_for_exact(m, psi, psi_out)
  !-----------------------------------------------------------------------
    USE mp,           ONLY : mp_sum
    USE mp_pools,     ONLY : nproc_pool, me_pool
    USE mp_exx,       ONLY : intra_egrp_comm, inter_egrp_comm, &
         nproc_egrp, me_egrp, negrp, my_egrp_id
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_STATUS_SIZE, MPI_DOUBLE_COMPLEX
#endif
    USE klist,        ONLY : xk, wk, nkstot, nks, qnorm
    USE wvfct,        ONLY : current_k
    !
    !
    IMPLICIT NONE
    !
    INTEGER :: m
    !!!!!!POSSIBLE ERROR HERE: PSI REALLY NEEDS LDA (NOT N)
    COMPLEX(DP) :: psi(npwx_exx*npol,m) 
    COMPLEX(DP) :: psi_out(npwx_local*npol,m)
    !COMPLEX(DP) :: psi(lda_exx(me_egrp+1,ik)*npol,m)
    !COMPLEX(DP) :: psi_out(lda_local(me_pool+1,ik)*npol,m) 

    INTEGER :: i, j, im, iproc, ig, ik, current_ik
    INTEGER :: prev, lda_max_local, prev_lda_exx

    INTEGER :: request_send(nproc_egrp), request_recv(nproc_egrp)
    INTEGER :: ierr
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
    INTEGER, EXTERNAL :: find_current_k

    current_ik = current_k
!    current_ik=find_current_k(ik, nkstot, nks)
    prev_lda_exx = sum( lda_exx(1:me_egrp,current_ik) )

    !send communication packets
    DO iproc=0, nproc_egrp-1
       IF ( comm_send_reverse(iproc+1,current_ik)%size.gt.0) THEN
          DO i=1, comm_send_reverse(iproc+1,current_ik)%size
             ig = comm_send_reverse(iproc+1,current_ik)%indices(i)
             ig = ig - prev_lda_exx
             
             DO im=1, m
                comm_send_reverse(iproc+1,current_ik)%msg(i,im) = psi(ig,im)
             END DO

          END DO

          !send the message
          CALL MPI_ISEND( comm_send_reverse(iproc+1,current_ik)%msg, &
               comm_send_reverse(iproc+1,current_ik)%size*m, MPI_DOUBLE_COMPLEX, &
               iproc, 100+iproc*nproc_egrp+me_egrp, &
               intra_egrp_comm, request_send(iproc+1), ierr )
          
       END IF
    END DO
        
    !begin recieving the communication packets
    DO iproc=0, nproc_egrp-1
       IF ( comm_recv_reverse(iproc+1,current_ik)%size.gt.0) THEN

          !recieve the message
          CALL MPI_IRECV( comm_recv_reverse(iproc+1,current_ik)%msg, &
               comm_recv_reverse(iproc+1,current_ik)%size*m, MPI_DOUBLE_COMPLEX, &
               iproc, 100+me_egrp*nproc_egrp+iproc, &
               intra_egrp_comm, request_recv(iproc+1), ierr )
          
       END IF
    END DO

    !assign psi
    DO iproc=0, nproc_egrp-1
       IF ( comm_recv_reverse(iproc+1,current_ik)%size.gt.0 ) THEN
          
          CALL MPI_WAIT(request_recv(iproc+1), istatus, ierr)

          DO i=1, comm_recv_reverse(iproc+1,current_ik)%size
             ig = comm_recv_reverse(iproc+1,current_ik)%indices(i)

             !set psi_out
             DO im=1, m
                psi_out(ig,im) = psi_out(ig,im) + &
                     comm_recv_reverse(iproc+1,current_ik)%msg(i,im)
             END DO

          END DO

       END IF
    END DO

    !now wait for everything to finish sending
    DO iproc=0, nproc_egrp-1
       IF ( comm_send_reverse(iproc+1,current_ik)%size.gt.0 ) THEN
          CALL MPI_WAIT(request_send(iproc+1), istatus, ierr)
       END IF
    END DO

  END SUBROUTINE deconstruct_for_exact







  !-----------------------------------------------------------------------
  SUBROUTINE change_data_structure(is_exx)
  !-----------------------------------------------------------------------
    ! ... generic, k-point version of vexx
    !
    USE io_files,       ONLY : nwordwfc, iunwfc, iunigk, iunigk_exx, seqopn
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega, at, bg, tpiba2
    USE wvfct,          ONLY : npwx, npw, igk, current_k, ecutwfc, g2kin
    USE control_flags,  ONLY : gamma_only
    USE klist,          ONLY : xk, nks, nkstot, ngk
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_bands,       ONLY : root_bgrp, me_bgrp, nproc_bgrp, ntask_groups, &
                               intra_bgrp_comm
    USE mp_exx,         ONLY : inter_egrp_comm, intra_egrp_comm, my_egrp_id, &
                               negrp, nproc_egrp, me_egrp, exx_mode, root_egrp
    USE gvect,              ONLY : ig_l2g, g, gg, ngm, ngm_g, gcutm, &
                                   mill,  nl, gstart, gkcut, &
                                   eigts1, eigts2, eigts3, ngl, igtongl, &
                                   gvect_init, deallocate_gvect_exx
    USE gvecs,          ONLY : ngms, gcutms, ngms_g, nls, gvecs_init, &
                               deallocate_gvecs
    USE uspp,           ONLY : nkb, okvan, vkb
    USE mp,             ONLY : mp_sum, mp_barrier, mp_bcast, mp_size, mp_rank
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE recvec_subs,        ONLY : ggen 
    USE fft_base,             ONLY : dfftp, dffts
    USE cellmd,             ONLY : lmovecell
    USE io_global,      ONLY : stdout
    USE stick_set,      ONLY : pstickset
    USE mp_pools
    !
    !
    IMPLICIT NONE
    !
    LOGICAL, intent(in) :: is_exx
    COMPLEX(DP), ALLOCATABLE :: work_space(:)
    !INTEGER :: comm
    INTEGER :: ik, i
    INTEGER :: ngm_, ngs_, ngw_
    LOGICAL exst

    !!!!!!!!!!
    IF (negrp.eq.1) RETURN
    !!!!!!!!!!

    CALL start_clock ('cds')

    !-----------------------------------------!
    ! Switch to the exx data structure        !
    !-----------------------------------------!

    CALL start_clock ('cds_fft')

    IF (first_data_structure_change) THEN
       allocate( ig_l2g_loc(ngm), g_loc(3,ngm), gg_loc(ngm) )
       allocate( mill_loc(3,ngm), nl_loc(ngm) )
       allocate( nls_loc(ngms) )
       ig_l2g_loc = ig_l2g
       g_loc = g
       gg_loc = gg
       mill_loc = mill
       nl_loc = nl
       nls_loc = nls
       ngm_loc = ngm
       ngm_g_loc = ngm_g
       gstart_loc = gstart
       ngms_loc = ngms
       ngms_g_loc = ngms_g
    END IF

    !try to generate the correct gvectors for the exact exchange calculation
    IF (is_exx) THEN
       exx_mode = 1
       !call data_structure( gamma_only )
       IF(first_data_structure_change)THEN
          dfftp_loc = dfftp
          dffts_loc = dffts
          !CALL pstickset( gamma_only, bg, gcutm, gkcut, gcutms, &
          !     dfftp, dffts, ngw_ , ngm_ , ngs_ , me_egrp, &
          !     root_egrp, nproc_egrp, intra_egrp_comm, ntask_groups, ionode, &
          !     stdout )
          CALL pstickset( gamma_only, bg, gcutm, gkcut, gcutms, &
               dfftp, dffts, ngw_ , ngm_ , ngs_ , me_egrp, &
               root_egrp, nproc_egrp, intra_egrp_comm, 1, ionode, &
               stdout, exx_mode = 1 )
          dfftp_exx = dfftp
          dffts_exx = dffts
          ngm = ngm_
          ngms = ngs_
       ELSE
          dfftp = dfftp_exx
          dffts = dffts_exx
          ngm = ngm_exx
          ngms = ngms_exx
       END IF
       call deallocate_gvect_exx()
       call deallocate_gvecs()
       call gvect_init( ngm , intra_egrp_comm )
       call gvecs_init( ngms , intra_egrp_comm )
    ELSE
       exx_mode = 2
       dfftp = dfftp_loc
       dffts = dffts_loc
       ngm = ngm_loc
       ngms = ngms_loc
       call deallocate_gvect_exx()
       call deallocate_gvecs()
       call gvect_init( ngm , intra_bgrp_comm )
       call gvecs_init( ngms , intra_bgrp_comm )
       exx_mode = 0
    END IF

    CALL stop_clock ('cds_fft')

    CALL start_clock ('cds_ggn')

    IF (first_data_structure_change) THEN
       CALL ggen( gamma_only, at, bg, intra_egrp_comm, no_global_sort = .FALSE. )
       allocate( ig_l2g_exx(ngm), g_exx(3,ngm), gg_exx(ngm) )
       allocate( mill_exx(3,ngm), nl_exx(ngm) )
       allocate( nls_exx(ngms) )
       ig_l2g_exx = ig_l2g
       g_exx = g
       gg_exx = gg
       mill_exx = mill
       nl_exx = nl
       nls_exx = nls
       ngm_exx = ngm
       ngm_g_exx = ngm_g
       gstart_exx = gstart
       ngms_exx = ngms
       ngms_g_exx = ngms_g
    ELSE IF ( is_exx ) THEN
       ig_l2g = ig_l2g_exx
       g = g_exx
       gg = gg_exx
       mill = mill_exx
       nl = nl_exx
       nls = nls_exx
       ngm = ngm_exx
       ngm_g = ngm_g_exx
       gstart = gstart_exx
       ngms = ngms_exx
       ngms_g = ngms_g_exx
    ELSE ! not is_exx
       ig_l2g = ig_l2g_loc
       g = g_loc
       gg = gg_loc
       mill = mill_loc
       nl = nl_loc
       nls = nls_loc
       ngm = ngm_loc
       ngm_g = ngm_g_loc
       gstart = gstart_loc
       ngms = ngms_loc
       ngms_g = ngms_g_loc
    END IF

    CALL stop_clock ('cds_ggn')



    CALL start_clock ('cds_npw')

    !get npwx
    IF ( is_exx.and.npwx_exx.gt.0 ) THEN
       npwx = npwx_exx
       !npw = npw_exx
       ngk = ngk_exx
    ELSE IF ( .not.is_exx.and.npwx_local.gt.0 ) THEN
       npwx = npwx_local
       !npw = npw_local
       ngk = ngk_local
    ELSE
       call n_plane_waves (ecutwfc, tpiba2, nks, xk, g, ngm, npwx, ngk)
    END IF

    CALL stop_clock ('cds_npw')




    CALL start_clock ('cds_igk')
    
    !get igk
    IF( first_data_structure_change ) THEN
       allocate(igk_exx(npwx),work_space(npwx))
       first_data_structure_change = .FALSE.
       IF ( nks.eq.1 ) THEN
          CALL gk_sort( xk, ngm, g, ecutwfc / tpiba2, npw, igk_exx, work_space )
       END IF
       IF ( nks > 1 ) THEN
          CALL seqopn( iunigk_exx, 'igk_exx', 'UNFORMATTED', exst )
          REWIND( iunigk_exx )
          DO ik = 1, nks
             CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk_exx, &
                  work_space )
             WRITE( iunigk_exx ) igk_exx
          END DO
          
       END IF
       DEALLOCATE( work_space )
    END IF

    CALL stop_clock ('cds_igk')




    !generate ngl and igtongl
    CALL start_clock ('cds_ngl')
    CALL gshells( lmovecell )
    CALL stop_clock ('cds_ngl')

    CALL stop_clock ('cds')

  END SUBROUTINE change_data_structure













  !-----------------------------------------------------------------------
  SUBROUTINE update_igk(is_exx)
  !-----------------------------------------------------------------------
    ! ... generic, k-point version of vexx
    !
    USE io_files,       ONLY : nwordwfc, iunwfc, iunigk, iunigk_exx, seqopn
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega, at, bg, tpiba2
    USE gvect,          ONLY : ngm, g, eigts1, eigts2, eigts3, ngl, igtongl
    USE wvfct,          ONLY : npwx, npw, igk, current_k, ecutwfc, g2kin
    USE control_flags,  ONLY : gamma_only
    USE klist,          ONLY : xk, nks, nkstot, ngk
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_exx,       ONLY : inter_egrp_comm, intra_egrp_comm, my_egrp_id, &
                               negrp, nproc_egrp, me_egrp, exx_mode
    USE gvect,              ONLY : ig_l2g, mill_g
    USE uspp,           ONLY : nkb, okvan, vkb
    USE mp,             ONLY : mp_sum, mp_barrier, mp_bcast, mp_size, mp_rank
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE recvec_subs,        ONLY : ggen 
    USE fft_base,             ONLY : dfftp, dffts
    USE cellmd,             ONLY : lmovecell
    USE mp_pools
    !
    !
    IMPLICIT NONE
    !
    LOGICAL, intent(in) :: is_exx
    COMPLEX(DP), ALLOCATABLE :: work_space(:)
    INTEGER :: comm
    INTEGER :: ik, i
    LOGICAL exst

    !!!!!!!!!!
    IF (negrp.eq.1) RETURN
    !!!!!!!!!!

    !get igk
    allocate(work_space(npwx))
    
    !NOTE: If nks > 1, should probably read igk from file to be faster
    ik = current_k
    IF(is_exx) THEN
       CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk_exx, work_space )
       !npw_exx = npw
    ELSE
       CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, work_space )
    END IF

    DEALLOCATE( work_space )

  END SUBROUTINE update_igk





  !------------------------------------------------------------------------
  SUBROUTINE communicate_exxbuff (ipair, request_send, request_recv)
    USE mp_exx,       ONLY : iexx_start, iexx_end, inter_egrp_comm, &
                               intra_egrp_comm, my_egrp_id, negrp, &
                               max_pairs, egrp_pairs
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_STATUS_SIZE, MPI_DOUBLE_COMPLEX
#endif
    USE io_global,      ONLY : stdout
    INTEGER, intent(in)      :: ipair
    INTEGER                  :: nrxxs
    INTEGER, intent(out)     :: request_send, request_recv
    INTEGER                  :: dest, sender, ierr, jnext, jnext_dest
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
    
    nrxxs= exx_fft%dfftt%nnr

    IF (ipair.lt.max_pairs) THEN
       
       IF (ipair.gt.1) THEN
#if defined(__MPI)
          CALL MPI_WAIT(request_send, istatus, ierr)
          CALL MPI_WAIT(request_recv, istatus, ierr)
#endif
       END IF
       
       sender = my_egrp_id + 1
       IF (sender.ge.negrp) sender = 0
       jnext = egrp_pairs(2,ipair+1,sender+1)

       dest = my_egrp_id - 1
       IF (dest.lt.0) dest = negrp - 1
       jnext_dest = egrp_pairs(2,ipair+1,dest+1)
       
#if defined(__MPI)
       CALL MPI_ISEND( exxbuff(:,:,jnext_dest), nrxxs*npol*nqs, &
            MPI_DOUBLE_COMPLEX, dest, 101, inter_egrp_comm, request_send, ierr )
#endif

#if defined(__MPI)
       CALL MPI_IRECV( exxbuff(:,:,jnext), nrxxs*npol*nqs, &
            MPI_DOUBLE_COMPLEX, sender, 101, inter_egrp_comm, &
            request_recv, ierr )
#endif

    END IF
    
  END SUBROUTINE communicate_exxbuff






  SUBROUTINE result_sum (n, m, data)
    USE mp_exx,       ONLY : iexx_start, iexx_end, inter_egrp_comm, &
                               intra_egrp_comm, my_egrp_id, negrp, &
                               max_pairs, egrp_pairs, max_contributors, &
                               contributed_bands, all_end, &
                               iexx_istart, iexx_iend, band_roots
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_STATUS_SIZE, MPI_DOUBLE_COMPLEX
#endif
    USE io_global,      ONLY : stdout
    USE mp,                   ONLY : mp_sum, mp_bcast

    INTEGER, INTENT(in) :: n, m
    COMPLEX(DP), INTENT(inout) :: data(n,m)
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
    COMPLEX(DP) :: recvbuf(n*max_contributors,(iexx_end-iexx_start+1))
    COMPLEX(DP) :: data_sum(n,m), test(negrp)
    INTEGER :: im, iegrp, ibuf, i, j, nsending(m)
    INTEGER :: ncontributing(m)
    INTEGER :: contrib_this(negrp,m), displs(negrp,m)

    INTEGER sendcount, sendtype, ierr, root, request(m)
    INTEGER sendc(negrp), sendd(negrp)

    IF (negrp.eq.1) RETURN

    !gather data onto the correct nodes
    CALL start_clock ('sum1')
    displs = 0
    ibuf = 0
    nsending = 0
    contrib_this = 0
    DO im=1, m
       
       IF(contributed_bands(im,my_egrp_id+1)) THEN
          sendcount = n
       ELSE
          sendcount = 0
       END IF

       root = band_roots(im)

       IF(my_egrp_id.eq.root) THEN
          !determine the number of sending processors
          ibuf = ibuf + 1
          ncontributing(im) = 0
          DO iegrp=1, negrp
             IF(contributed_bands(im,iegrp)) THEN
                ncontributing(im) = ncontributing(im) + 1
                contrib_this(iegrp,im) = n
                IF(iegrp.lt.negrp) displs(iegrp+1,im) = displs(iegrp,im) + n
                nsending(im) = nsending(im) + 1
             ELSE
                contrib_this(iegrp,im) = 0
                IF(iegrp.lt.negrp) displs(iegrp+1,im) = displs(iegrp,im)
             END IF
          END DO
       END IF
       
       CALL MPI_IGATHERV(data(:,im), sendcount, MPI_DOUBLE_COMPLEX, &
            recvbuf(:,ibuf), contrib_this(:,im), &
            displs(:,im), MPI_DOUBLE_COMPLEX, &
            root, inter_egrp_comm, request(im), ierr)

    END DO
    CALL stop_clock ('sum1')

    CALL start_clock ('sum2')
    DO im=1, m
       CALL MPI_WAIT(request(im), istatus, ierr)
    END DO
    CALL stop_clock ('sum2')

    !perform the sum
    CALL start_clock ('sum3')
    DO im=iexx_istart(my_egrp_id+1), iexx_iend(my_egrp_id+1)
       IF(im.eq.0)exit
       data_sum(:,im) = 0._dp
       ibuf = im - iexx_istart(my_egrp_id+1) + 1
       DO j=1, nsending(im)
          DO i=1, n
             data_sum(i,im) = data_sum(i,im) + recvbuf(i+n*(j-1),ibuf)
          END DO
       END DO
    END DO
    CALL stop_clock ('sum3')

    CALL start_clock ('sum4')
    sendd = 0
    DO iegrp=1, negrp
       IF (iexx_istart(iegrp).eq.0) THEN
          sendc(iegrp) = 0
          IF(iegrp.lt.negrp) sendd(iegrp+1) = sendd(iegrp)
       ELSE
          sendc(iegrp) = n*(iexx_iend(iegrp) - iexx_istart(iegrp) + 1)
          IF(iegrp.lt.negrp) sendd(iegrp+1) = sendd(iegrp) + sendc(iegrp)
       END IF
    END DO
    CALL stop_clock ('sum4')
    CALL start_clock ('sum5')
    CALL MPI_Allgatherv(data_sum(1,max(1,iexx_istart(my_egrp_id+1))), &
         sendc(my_egrp_id+1), &
         MPI_DOUBLE_COMPLEX, &
         data, sendc, sendd, &
         MPI_DOUBLE_COMPLEX, inter_egrp_comm, ierr)
    CALL stop_clock ('sum5')

  END SUBROUTINE result_sum









!>>>
!-----------------------------------------------------------------------
END MODULE exx
!-----------------------------------------------------------------------
