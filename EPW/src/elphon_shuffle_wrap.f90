  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE elphon_shuffle_wrap
  !-----------------------------------------------------------------------
  !
  ! Electron-phonon calculation with Wannier functions: load all phonon q's
  !
  ! 09/2009 This subroutine is the main driver of the electron-phonon 
  ! calculation. It first calculates the electron-phonon matrix elements
  ! on the coarse mesh and then passes the data off to ephwann_shuffle
  ! to perform the interpolation.
  !
  !-----------------------------------------------------------------------
  !
#ifdef __PARA
  USE mp_global,     ONLY : my_pool_id, inter_pool_comm, root_pool, &
                            intra_pool_comm,npool
  USE mp_world,      ONLY : mpime, root
  USE mp,            ONLY : mp_barrier, mp_bcast
  USE io_global,     ONLY : ionode_id
#endif
  USE us,            ONLY : nqxq, dq, qrad
  USE gvect,         ONLY : gcutm
  USE cellmd,        ONLY : cell_factor
  USE uspp_param,    ONLY : lmaxq, nbetam
  USE io_files,      ONLY : prefix, tmp_dir
  USE wavefunctions_module, ONLY: evc
  USE ions_base,     ONLY : nat, nsp, tau, ityp
  USE control_flags, ONLY : iverbosity
  USE io_global,     ONLY : stdout, ionode
  USE io_epw,        ONLY : iuepb 
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : et, xk, ibrav, igk, nks, nbnd, nkstot, ngm
  USE cell_base,     ONLY : at, bg
  USE symm_base,     ONLY : irt, s, nsym, ftau, sname, invs, &
                            s_axis_to_cart, sr, nrot, copy_sym, set_sym_bl, find_sym, inverse_s
  USE start_k,       ONLY : nk1, nk2, nk3
  USE phcom,         ONLY : dpsi, dvpsi, evq, nq1, nq3, nq2 
  USE qpoint,        ONLY : igkq
  USE control_lr,    ONLY : lgamma 
  USE qpoint,        ONLY : xq
  USE modes,         ONLY : nmodes
  USE lr_symm_base,  ONLY : minus_q, rtau, gi, gimq, irotmq, nsymq, invsymq
  USE epwcom,        ONLY : epbread, epbwrite, epwread, phinterp, &
                            phonselfen, elecselfen, nbndsub, elinterp,    &
                            iswitch, kmaps, nest_fn, eig_read, &
                            band_plot, specfun, dvscf_dir, lpolar
  USE elph2,         ONLY : epmatq, dynq, sumr, et_all, xk_all, et_mb, et_ks, &
                            zstar, epsi, cu, cuq, lwin, lwinq, bmat
  USE constants_epw, ONLY : ryd2ev 
  USE fft_base,      ONLY : dfftp
  USE control_ph,    ONLY : u_from_file
  USE noncollin_module, ONLY : m_loc
  USE iotk_module,   ONLY : iotk_open_read, iotk_scan_dat, iotk_free_unit,&
                            iotk_close_read
#ifdef __NAG
  USE f90_unix_io,    ONLY : flush
#endif
  implicit none

  !
  integer :: sym_smallq(48) 
  !
  real(kind=DP), allocatable :: xqc_irr(:,:), wqlist_irr(:), xqc(:,:), wqlist(:)
  ! the qpoints in the irr wedge
  ! the corresponding weigths
  ! the qpoints in the uniform mesh
  ! the corresponding weigths
  integer :: nqc_irr, nqc, max, nqxq_tmp, ibnd, ik, ios, &
             dummy1, dummy2, ik_start, ik_stop
  ! number of qpoints in the irreducible wedge
  ! number of qpoints on the uniform grid
  ! 
  ! symmetry-related variables
  !
  integer :: gmapsym(ngm,48)
  !  correspondence G -> S(G)
  complex(kind=DP) :: eigv (ngm, 48)
  ! e^{ iGv} for 1...nsym (v the fractional translation)
  complex(kind=DP) :: cz1( nmodes, nmodes), cz2(nmodes, nmodes)
  !  the eigenvectors for the first q in the star
  !  the rotated eigenvectors, for the current q in the star
  !
  integer :: nq, isq (48), imq 
  ! degeneracy of the star of q
  ! index of q in the star of a given sym.op.
  ! index of -q in the star of q (0 if not present)
  integer :: sym_sgq(48)
  ! the symmetries giving the q point iq in the star
  real(kind=DP) :: sxq (3, 48), et_tmp(nbnd, nkstot)
  ! list of vectors in the star of q
  integer :: i, j, iq, iq_irr, isym, &
     iq_first, jsym, ism1, nsq, ipol, jpol, ierr, iunpun
  real(kind=DP) xq0(3), aq(3), saq(3), raq(3), ft1, ft2, ft3
!  REAL(DP) :: w2(3*nat) ! dummy here
  logical :: sym(48),  eqvect_strict, nog, symmo, exst
  character (len=256) :: tempfile, dirname,filename
#ifdef __PARA
  character (len=3) :: filelab
#endif
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  CALL start_clock ( 'elphon_wrap' )
  !
  IF ( elinterp .and. (.not.phinterp ) ) CALL errore &
        ('elphon_shuffle_wrap','elinterp requires phinterp' ,1)
  !
  ! READ qpoint list from stdin
  !
#ifdef __PARA
  IF (mpime.eq.ionode_id) &
#endif
  READ(5,*) nqc_irr
#ifdef __PARA
  CALL mp_bcast (nqc_irr, ionode_id, inter_pool_comm)
  CALL mp_bcast (nqc_irr, root_pool, intra_pool_comm)
#endif
  allocate ( xqc_irr(3,nqc_irr), wqlist_irr(nqc_irr) )
  allocate ( xqc(3,nq1*nq2*nq3), wqlist(nq1*nq2*nq3) )
  !  
#ifdef __PARA
  IF (mpime.eq.ionode_id) then
#endif
    DO iq = 1, nqc_irr
      READ (5,*) xqc_irr (:,iq), wqlist_irr (iq)
    ENDDO
#ifdef __PARA
  ENDIF
  CALL mp_bcast (xqc_irr, ionode_id, inter_pool_comm)
  CALL mp_bcast (xqc_irr, root_pool, intra_pool_comm)
#endif
  !
  ! fix for uspp
  max = nqxq
  DO iq = 1, nqc_irr
     nqxq_tmp = INT( ( (sqrt(gcutm) + sqrt(xqc_irr(1,iq)**2 + &
          xqc_irr(2,iq)**2 + xqc_irr(3,iq)**2) ) &
          / dq + 4) * cell_factor )
     IF (nqxq_tmp .gt. max)  max = nqxq_tmp
  ENDDO
  IF (max .gt. nqxq) then
     IF (allocated(qrad)) deallocate(qrad)
     allocate (qrad (max, nbetam*(nbetam+1)/2,lmaxq, nsp))
  ENDIF
  IF (nkstot .ne. nk1*nk2*nk3 ) &
       CALL errore('elphon_shuffle_wrap','nscf run inconsistent with epw input',1)  
  !
  ! READ in external electronic eigenvalues. e.g. GW 
  !
  IF ( .not. ALLOCATED(et_ks) ) ALLOCATE(et_ks(nbnd,nks))
  IF ( .not. ALLOCATED(et_mb) ) ALLOCATE(et_mb(nbnd,nks))
  et_ks(:,:) = 0.d0
  et_mb(:,:) = 0.d0
  IF (eig_read) then
#ifdef __PARA
  IF (mpime.eq.ionode_id) then
#endif
   WRITE (stdout,'(5x,a,i5,a,i5,a)') "Reading external electronic eigenvalues (", &
        nbnd, ",", nkstot,")"
   tempfile=trim(prefix)//'.eig'
   OPEN(1, file=tempfile, form='formatted', action='read', iostat=ios)
   IF (ios /= 0) CALL errore ('elphon_shuffle_wrap','error opening' // tempfile, 1)
    DO ik = 1, nkstot
       DO ibnd = 1, nbnd
          READ (1,*) dummy1, dummy2, et_tmp (ibnd,ik)
          IF (dummy1.ne.ibnd) CALL errore('elphon_shuffle_wrap', "Incorrect eigenvalue file", 1)
          IF (dummy2.ne.ik)   CALL errore('elphon_shuffle_wrap', "Incorrect eigenvalue file", 1)
       ENDDO
    ENDDO
    CLOSE(1)
    ! from eV to Ryd
    et_tmp = et_tmp / ryd2ev
#ifdef __PARA
    ENDIF
    CALL mp_bcast (et_tmp, ionode_id, inter_pool_comm)
    CALL mp_bcast (et_tmp, root_pool, intra_pool_comm)
#endif
    !
    CALL ckbounds(ik_start, ik_stop)
    et_ks(:,:)  = et(:,1:nks)
    et(:,1:nks) = et_tmp(:,ik_start:ik_stop)
    et_mb(:,:)  = et(:,1:nks)
 ENDIF
  !
  ! compute coarse grid dipole matrix elements.  Very fast 
  CALL compute_pmn_para
  !
  !  gather electronic eigenvalues for subsequent shuffle
  !  
  allocate ( xk_all( 3, nkstot) , et_all( nbnd, nkstot) )
  CALL poolgather (    3, nkstot, nks, xk(:,1:nks),      xk_all)
  CALL poolgather ( nbnd, nkstot, nks, et(1:nbnd,1:nks), et_all)
  !
  !
  IF (.not.kmaps) then
     CALL start_clock('kmaps')
     xq0(:)=0.d0
     CALL createkmap_pw2(xk_all,nkstot, xq0)
     CALL stop_clock('kmaps')
     CALL print_clock('kmaps')
  ELSE
     ! 
     ! 26/06/2012 RM
     ! if we do not have epmatq already on file then epbread=.false.
     ! .kgmap is used from disk and .kmap is regenerated for each q-point 
     ! 
     WRITE (stdout, '(/5x,a)')  'Using kmap and kgmap from disk'
  ENDIF
#ifdef __PARA
  CALL mp_barrier(inter_pool_comm)
#endif
  !  if we start with a gamma point calculation, ../PW/set_kplusq.f90
  !  is not active and the gmap has not been produced...
  !
  IF (lgamma) CALL errore &
    ('elphon_shuffle_wrap','tshuffle2 requires q!=0 starting nscf calculation',1)
  !
  !  allocate dynamical matrix and ep matrix for all q's
  !
  allocate ( dynq (nmodes, nmodes, nq1*nq2*nq3), &
       epmatq (nbnd, nbnd, nks, nmodes, nq1*nq2*nq3), &
       epsi(3,3), zstar(3,3,nat), bmat(nbnd, nbnd, nks, nq1*nq2*nq3), &
       cu ( nbnd, nbndsub, nks), cuq ( nbnd, nbndsub, nks), & 
       lwin ( nbnd, nks ), lwinq ( nbnd, nks ) )
  !
  epsi=0.d0
  zstar=0.d0
  !
  ! SP: The symmetries are now consistent with QE 5. This means that the order of the q in the star
  !     should be the same as in the .dyn files produced by QE 5.
  ! 
  !     First we start by setting up the lattice & crystal symm. as done in PH/init_representations.f90
  ! 
  ! ~~~~~~~~ setup bravais lattice symmetry ~~~~~~~~
  CALL set_sym_bl ( ) ! This should define the s matrix
  WRITE(stdout, '(5x,a,i3)') "Symmetries of bravais lattice: ", nrot
  !
  ! ~~~~~~~~ setup crystal symmetry ~~~~~~~~ 
  CALL find_sym ( nat, tau, ityp, dfftp%nr1,dfftp%nr2,dfftp%nr3, .false., m_loc )
  WRITE(stdout, '(5x,a,i3)') "Symmetries of crystal:         ", nsym
  !   
  ! The following loop is required to propertly set up the symmetry matrix s. 
  ! We here copy the calls made in PH/init_representations.f90 to have the same s as in QE 5.
  DO iq_irr = 1, nqc_irr
    xq = xqc_irr(:,iq_irr)
    CALL set_small_group_of_q(nsymq,invsymq,minus_q)
    CALL sgam_ph_new (at, bg, nsym, s, irt, tau, rtau, nat)
    CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
  ENDDO
  ! 
  ! CV: if we read the .fmt files we don't need to read the .epb anymore
  !
  IF (.not. epbread .and. .not. epwread) THEN
    ! 
    ! In the loop over irr q-point, we need to read the pattern that
    ! correspond to the dvscf file computed with QE 5.
    !
    nqc = 0
    iq_first = 1
    DO iq_irr = 1, nqc_irr
      u_from_file = .TRUE.
      !tmp_dir_ph = './_ph0/'
      !
      !  read the displacement patterns
      !
      IF (u_from_file) THEN
         ierr=0
         IF ( ionode ) THEN
         !
         ! ... look for an empty unit (only ionode needs it)
         !
            CALL iotk_free_unit( iunpun, ierr )
         !
         END IF
         dirname = TRIM(dvscf_dir) // TRIM( prefix ) // '.phsave'
         filename= TRIM( dirname ) // '/patterns.' // &
                   TRIM(int_to_char(iq_irr)) // '.xml'
         INQUIRE( FILE=TRIM(filename), EXIST=exst )
         IF (.NOT.exst) CALL errore( 'ph_restart_set_filename ', &
                   'cannot open file for reading or writing', ierr )
         CALL iotk_open_read( iunpun, FILE = TRIM( filename ), &
                                          BINARY = .FALSE., IERR = ierr )
         CALL read_modes(iunpun,iq_irr, ierr )
         IF (ierr /= 0) CALL errore('epw_setup', 'problem with modes file',1)
         IF (ionode) CALL iotk_close_read( iunpun )
      ENDIF
      !  
      WRITE(stdout,'(//5x,a)') repeat('=',67) 
      WRITE(stdout,'(5x,"irreducible q point # ",i4)') iq_irr
      WRITE(stdout,'(5x,a/)') repeat('=',67) 
      CALL flush(6)
      !
      xq = xqc_irr(:,iq_irr)
      !
      ! SP : The following is largely inspiered by PH/q2qstar.f90
      ! 
      ! ~~~~~~~~ setup small group of q symmetry ~~~~~~~~ 
      !
      minus_q = .true.
      sym = .false.
      sym(1:nsym) = .true.
      CALL smallg_q(xq, 0, at, bg, nsym, s, ftau, sym, minus_q) ! s is intent(in)
      !
      ! SP: Notice that the function copy_sym reshuffle the s matrix for each irr_q.  
      !     This is why we then need to call gmap_sym for each irr_q [see below]. 
      nsymq = copy_sym(nsym, sym)    
      !
      ! Recompute the inverses as the order of sym.ops. has changed
      CALL inverse_s ( )       
      call s_axis_to_cart ( )
      !
      ! This computes gi, gimq
      call set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
      WRITE(stdout, '(5x,a,i3)') "Symmetries of small group of q:", nsymq
      IF(minus_q) WRITE(stdout, '(10x,a)') "in addition sym. q -> -q+G:"
      ! 
      ! Finally this does some of the above again and also computes rtau...
      CALL sgam_ph_new(at, bg, nsym, s, irt, tau, rtau, nat)
      !
      ! ######################### star of q #########################
      ! 
      sym_smallq(:) = 0
      CALL star_q2(xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, .true., sym_smallq )
      !
      ! The reason for xq instead of xq0 in the above is because xq is pass to QE through module  
      xq0 = xq
      !
      !  determine the G vector map S(G) -> G 
      !  SP: The mapping need to be done for each irr_q because the QE 5 symmetry routine
      !      reshuffle the s matrix for each irr_q [putting the sym of the small group of q first].
      !
      !  [I checked that gmapsym(gmapsym(ig,isym),invs(isym)) = ig]
      CALL gmap_sym ( nsym, s, ftau, gmapsym, eigv, invs)
      !
      !  Re-set the variables needed for the pattern representation
      !  and the symmetries of the small group of irr-q
      !  (from phq_setup.f90)
      !
      DO isym = 1, nsym
        sym (isym) = .true.
      ENDDO
      !
      CALL sgam_ph_new (at, bg, nsym, s, irt, tau, rtau, nat)
      !
#ifdef __PARA
     IF ( .not. allocated(sumr) ) allocate ( sumr(2,3,nat,3) )
     IF (mpime.eq.ionode_id) then
#endif
        CALL readmat_shuffle2 ( iq_irr, nqc_irr, nq, iq_first, sxq, imq,isq,&
                              invs, s, irt, rtau)
#ifdef __PARA
     ENDIF
     CALL mp_barrier(inter_pool_comm)
     CALL mp_barrier(intra_pool_comm)
     CALL mp_bcast (zstar, ionode_id, inter_pool_comm)
     CALL mp_bcast (zstar, root_pool, intra_pool_comm)
     CALL mp_bcast (epsi, ionode_id, inter_pool_comm)
     CALL mp_bcast (epsi, root_pool, intra_pool_comm)
     CALL mp_bcast (dynq, ionode_id, inter_pool_comm)
     CALL mp_bcast (dynq, root_pool, intra_pool_comm)
     CALL mp_bcast (sumr, ionode_id, inter_pool_comm)
     CALL mp_bcast (sumr, root_pool, intra_pool_comm)
#endif
      !
      ! now dynq is the cartesian dyn mat (NOT divided by the masses)
      !
      minus_q = (iswitch .gt. -3)  
      !
      !  loop over the q points of the star 
      !
      DO iq = 1, nq
        ! SP: First the vlocq need to be initialized propertly with the first
        !     q in the star
        xq = xq0         
        CALL epw_init(.false.)
        !
        ! retrieve the q in the star
        xq = sxq(:,iq)                               
        !
        ! and populate the uniform grid
        nqc = nqc + 1
        xqc(:,nqc) = xq
        !
        IF (iq.eq.1) write(6,*)
        WRITE( stdout, 5) nqc, xq
        !
        !  prepare the gmap for the refolding
        !
        CALL createkmap ( xq )                      
        !
        IF (iverbosity.eq.1) then
          !
          !   description of symmetries
          !
          WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
          CALL s_axis_to_cart() ! give sr (:,:, isym)
          DO isym = 1, nsym
            WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isym, sname(isym)
          !  CALL s_axis_to_cart (s(1,1,isym), sr, at, bg)
            IF (ftau(1,isym).ne.0.or.ftau(2,isym).ne.0.or.ftau(3,isym).ne.0) then
                ft1 = at(1,1)*ftau(1,isym)/dfftp%nr1 + at(1,2)*ftau(2,isym)/dfftp%nr2 + &
                      at(1,3)*ftau(3,isym)/dfftp%nr3
                ft2 = at(2,1)*ftau(1,isym)/dfftp%nr1 + at(2,2)*ftau(2,isym)/dfftp%nr2 + &
                      at(2,3)*ftau(3,isym)/dfftp%nr3
                ft3 = at(3,1)*ftau(1,isym)/dfftp%nr1 + at(3,2)*ftau(2,isym)/dfftp%nr2 + &
                      at(3,3)*ftau(3,isym)/dfftp%nr3
                WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                      &        " )    f =( ",f10.7," )")') &
                      isym, (s(1,ipol,isym),ipol=1,3), dble(ftau(1,isym))/dble(dfftp%nr1)
                WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )")') &
                            (s(2,ipol,isym),ipol=1,3), dble(ftau(2,isym))/dble(dfftp%nr2)
                WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )"/)') &
                            (s(3,ipol,isym),ipol=1,3), dble(ftau(3,isym))/dble(dfftp%nr3)
                WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7, &
                      &        " )    f =( ",f10.7," )")') &
                      isym, (sr(1,ipol,isym),ipol=1,3), ft1
                WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                            (sr(2,ipol,isym),ipol=1,3), ft2
                WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)') &
                            (sr(3,ipol,isym),ipol=1,3), ft3
            ELSE
                WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), " )")') &
                                       isym,  (s (1, ipol, isym) , ipol = 1,3)
                WRITE( stdout, '(17x," (",3(i6,5x)," )")')  (s(2,ipol,isym), ipol=1,3)
                WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') (s(3,ipol,isym), ipol=1,3)
                WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7," )")') &
                                              isym,  (sr (1, ipol,isym) , ipol = 1, 3)
                WRITE( stdout, '(17x," (",3f11.7," )")')  (sr (2, ipol,isym) , ipol = 1, 3)
                WRITE( stdout, '(17x," (",3f11.7," )"/)') (sr (3, ipol,isym) , ipol = 1, 3)
            ENDIF
          ENDDO
          !
        ENDIF
        !
        ! isq(isym)=iq means: when we apply symmetry isym to the originating q 
        ! of the star, we get the iq-th member of the star. There are as many 
        ! matches as the degeneracy of the star.
        !
        ! We now need to pick up the q in the small group of q* so that Sxq0+G=iq with G=0.
        ! If we choose another element in the small group
        ! the actual q-point may be Sq+G and we screw up the q-vector below to generate
        ! k+q from k and for the KB projectors
        !
        nsq = 0 ! nsq is the degeneracy of the small group for this iq in the star
        !
        DO jsym = 1, nsym
          IF ( isq(jsym) .eq. iq ) then
             nsq = nsq + 1
             sym_sgq(nsq) = jsym
          ENDIF
        ENDDO
        IF ( nsq*nq .ne. nsym ) CALL errore ('elphon_shuffle_wrap', 'wrong degeneracy', iq)
        ! 
        IF (iverbosity.eq.1) then
          !
          WRITE(stdout,*) 'iq, i, isym, nog, symmo'
          DO i = 1, nsq
            !
            isym = sym_sgq(i)
            ism1 = invs (isym)
            !
            !  check for G such that Sq = q* + G 
            ! 
            aq  = xq0
            saq = xq
            CALL cryst_to_cart (1, aq, at, -1)
            DO j = 1, 3
              raq (j) = s (j, 1, ism1) * aq (1) &
                      + s (j, 2, ism1) * aq (2) &
                      + s (j, 3, ism1) * aq (3)
            ENDDO
            CALL cryst_to_cart (1, saq, at, -1)
            nog = eqvect_strict( raq, saq) 
            !
            !  check whether the symmetry belongs to a symmorphic group
            !
            symmo = (ftau(1,isym).eq.0 .and. ftau(2,isym).eq.0 .and. ftau(3,isym).eq.0)
            !
            WRITE(stdout,'(3i5,2a)') iq, i, isym, nog, symmo
            !
          ENDDO  
          !
        ENDIF
        ! 
        ! SP: We now need to select one symmetry among the small group of q (i.e. that respect 
        !     Sq0+G=q ) that has G=0. There should always be such symmetry. 
        !     We enforce this for later easiness. 
        ! 
        aq = xq0
        saq = xq
        call cryst_to_cart (1, aq, at, - 1)
        CALL cryst_to_cart (1, saq, at, -1)
        !write(*,*)'xq0 ',aq
        !write(*,*)'xq ',saq
        DO jsym=1, nsq
          ism1 = invs (sym_sgq(jsym))
          raq = 0.d0
          DO ipol = 1, 3
             DO jpol = 1, 3
                raq (ipol) = raq (ipol) + s (ipol, jpol, ism1) * aq (jpol)
             ENDDO
          ENDDO
          nog = eqvect_strict( raq, saq)
          IF (nog) THEN ! This is the symmetry such that Sq=q
            isym = sym_sgq(jsym)
            EXIT
          ENDIF
          ! If we enter into that loop it means that we have not found 
          ! such symmetry within the small group of Q. 
          IF (jsym == nsq) THEN
            CALL errore( 'elphon_shuffle_wrap ', 'No sym. such that Sxq0=iq was found in the sgq !', 1 )
          ENDIF
        ENDDO
        !
        !DBSP
        !write(*,*)'isym ',isym 
        !END
        !
        !
        CALL loadumat ( nbnd, nbndsub, nks, nkstot, xq, cu, cuq, lwin, lwinq )
        !
        ! Calculate overlap U_k+q U_k^\dagger
        IF (lpolar) CALL compute_bmn_para3 ( nbnd, nbndsub, nks, cu, cuq,bmat(:,:,:,nqc) )
        !
        !   calculate the sandwiches
        !
        ! a more accurate way of doing is to symmetrize the matrix element w.r.t.
        ! the small group of the given q in the star. I'm not doint this here.
        ! (but I checked that even without symm the result of full zone and irr zone
        ! are equal to 5+ digits)
        ! For any volunteers, please write to giustino@civet.berkeley.edu
        !
        CALL elphon_shuffle ( iq_irr, nqc_irr, nqc, gmapsym, eigv, isym, invs, xq0, .false. )
        !
        !  bring epmatq in the mode representation of iq_first, 
        !  and then in the cartesian representation of iq
        !
        CALL rotate_eigenm ( iq_first, iq, nqc, isym, nsym, s, invs, irt, &
           rtau, xq, isq, cz1, cz2 )
        !
        CALL rotate_epmat ( cz1, cz2, xq, nqc, lwin, lwinq )
  !DBSP
  !      write(*,*)'epmatq(:,:,2,:,nqc)',SUM(epmatq(:,:,2,:,nqc))
  !      write(*,*)'epmatq(:,:,2,:,nqc)**2',SUM((REAL(REAL(epmatq(:,:,2,:,nqc))))**2)+&
  !        SUM((REAL(AIMAG(epmatq(:,:,2,:,nqc))))**2)
  !END
        ! SP: Now we treat separately the case imq == 0
        IF (imq .eq. 0) then
          !
          ! SP: First the vlocq need to be initialized propertly with the first
          !     q in the star
          xq = -xq0
          CALL epw_init(.false.)
          !
          ! retrieve the q in the star
          xq = -sxq(:,iq)
          !
          ! and populate the uniform grid
          nqc = nqc + 1
          xqc(:,nqc) = xq
          !
          IF (iq.eq.1) write(stdout,*)
          WRITE( stdout, 5) nqc, xq
          !
          !  prepare the gmap for the refolding
          !
          CALL createkmap ( xq )
          !
          xq0 = -xq0
          !
          CALL elphon_shuffle ( iq_irr, nqc_irr, nqc, gmapsym, eigv, isym, invs, xq0, .true. )
          !
          !  bring epmatq in the mode representation of iq_first, 
          !  and then in the cartesian representation of iq
          !
          CALL rotate_eigenm ( iq_first, iq, nqc, isym, nsym, s, invs, irt, &
             rtau, xq, isq, cz1, cz2 )
          !
          CALL rotate_epmat ( cz1, cz2, xq, nqc, lwin, lwinq )
          !
    !DBSP
    !      write(*,*)'epmatq(:,:,2,:,nqc)',SUM(epmatq(:,:,2,:,nqc))
    !      write(*,*)'epmatq(:,:,2,:,nqc)**2',SUM((REAL(REAL(epmatq(:,:,2,:,nqc))))**2)+&
    !        SUM((REAL(AIMAG(epmatq(:,:,2,:,nqc))))**2)
    !END
          xq0 = -xq0
        ENDIF ! end imq == 0  
        !
      ENDDO
      !
      iq_first = iq_first + nq
      if (imq .eq. 0) iq_first = iq_first + nq
      !
    ENDDO ! irr-q loop
    ! 
    !
    IF (nqc.ne.nq1*nq2*nq3) &
       CALL errore ('elphon_shuffle_wrap','nqc .ne. nq1*nq2*nq3',nqc)
    wqlist = dble(1)/dble(nqc)
    !
  ENDIF
  !
  IF ( epbread .or. epbwrite ) THEN
    !
    ! write the e-ph matrix elements and other info in the Bloch representation
    ! (coarse mesh)
    ! in .epb files (one for each pool)
    !
    tempfile = trim(tmp_dir) // trim(prefix) // '.epb' 
#ifdef __PARA
     CALL set_ndnmbr (0,my_pool_id+1,1,npool,filelab)
     tempfile = trim(tmp_dir) // trim(prefix) // '.epb' // filelab
#endif
     !
     IF (epbread)  THEN
        inquire(file = tempfile, exist=exst)
        IF (.not. exst ) CALL errore( 'elphon_shuffle_wrap', 'epb files not found ', 1)
        OPEN  (iuepb, file = tempfile, form = 'unformatted')
        WRITE(stdout,'(/5x,"Reading epmatq from .epb files"/)') 
        READ  (iuepb) nqc, xqc, et, dynq, epmatq, zstar, epsi
        CLOSE (iuepb)
     ENDIF
     !
     IF (epbwrite) THEN
        OPEN  (iuepb, file = tempfile, form = 'unformatted')
        WRITE(stdout,'(/5x,"Writing epmatq on .epb files"/)') 
        WRITE (iuepb) nqc, xqc, et, dynq, epmatq, zstar, epsi
        CLOSE (iuepb)
     ENDIF
  ENDIF
  !
  IF ( .not.epbread .and. epwread ) THEN
  !  CV: need dummy nqc, xqc for the ephwann_shuffle call
     nqc=1
     xqc=0.d0
     WRITE(stdout,'(/5x,"Do not need to read .epb files; read .fmt files"/)')
  !
  ENDIF
  !
  !
#ifdef __PARA
  CALL mp_barrier(inter_pool_comm)
#endif
  !
  !   now dynq is the cartesian dyn mat ( NOT divided by the masses)
  !   and epmatq is the epmat in cartesian representation (rotation in elphon_shuffle)
  !
  !
  ! free up some memory
  !
  IF ( ASSOCIATED (evq)  ) NULLIFY    (evq)
  IF ( ALLOCATED  (evc)  ) DEALLOCATE (evc)
  IF ( ASSOCIATED (igkq) ) NULLIFY    (igkq)
  IF ( ALLOCATED  (igk)  ) DEALLOCATE (igk)
  IF ( ALLOCATED  (dvpsi)) DEALLOCATE (dvpsi)
  IF ( ALLOCATED  (dpsi) ) DEALLOCATE (dpsi)
  IF ( ALLOCATED  (sumr) ) DEALLOCATE (sumr)
  IF ( ALLOCATED (cu) )    DEALLOCATE (cu)
  IF ( ALLOCATED (cuq) )   DEALLOCATE (cuq)
  IF ( ALLOCATED (lwin) )  DEALLOCATE (lwin)
  IF ( ALLOCATED (lwinq) ) DEALLOCATE (lwinq)
  IF ( ALLOCATED (bmat) )  DEALLOCATE (bmat)
  ! 
  CALL stop_clock ( 'elphon_wrap' )
!DBSP
!  DO iq = 1, nqc
!    write(*,*) iq, xqc(:,iq)
!    write(*,*)'epmatq(:,:,2,:,iq)',SUM(epmatq(:,:,2,:,iq))
!    write(*,*)'epmatq(:,:,2,:,iq)**2',SUM((REAL(REAL(epmatq(:,:,2,:,iq))))**2)+&
!               SUM((REAL(AIMAG(epmatq(:,:,2,:,iq))))**2)
!  ENDDO
!END
  !
  ! the electron-phonon wannier interpolation
  !
  CALL ephwann_shuffle ( nqc, xqc )
  !
5 format (8x,"q(",i5," ) = (",3f12.7," )") 
  !
  RETURN
  END SUBROUTINE elphon_shuffle_wrap
  !
  !---------------------------------------------------------------------------
  SUBROUTINE irotate ( x, s, sx)
  !---------------------------------------------------------------------------
  !
  ! a simple symmetry operation in crystal coordinates ( s is integer!)
  !
  USE kinds, ONLY : DP
  implicit none
  real(kind=DP) x(3), sx(3)
  integer :: s(3,3),i
  !
  DO i = 1, 3
     sx (i) = dble(s (i, 1)) * x (1) &
            + dble(s (i, 2)) * x (2) &
            + dble(s (i, 3)) * x (3)
  ENDDO
  !
  RETURN
  end
  !---------------------------------------------------------------------------
  logical function eqvect_strict (x, y)
  !-----------------------------------------------------------------------
  !
  !   This function test if two tridimensional vectors are equal
  !
  USE kinds
  implicit none
  real(kind=DP) :: x (3), y (3)
  ! input: input vector
  ! input: second input vector
  real(kind=DP) :: accep
  ! acceptance PARAMETER
  PARAMETER (accep = 1.0d-5)
  !
  eqvect_strict = abs( x(1)-y(1) ).lt.accep .and. &
                  abs( x(2)-y(2) ).lt.accep .and. &
                  abs( x(3)-y(3) ).lt.accep
  END FUNCTION eqvect_strict
  !---------------------------------------------------------------------------
  SUBROUTINE read_modes(iunpun,current_iq, ierr )
  !---------------------------------------------------------------------------
  !
  ! This routine reads the displacement patterns.
  !
  USE modes,        ONLY : nirr, npert, u, name_rap_mode, num_rap_mode
  USE lr_symm_base, ONLY : minus_q, nsymq  
  USE iotk_module,  ONLY : iotk_index, iotk_scan_dat, iotk_scan_begin, &
                           iotk_scan_end
  USE io_global,    ONLY : ionode, ionode_id
  USE mp_images,    ONLY : intra_image_comm
  USE mp,           ONLY : mp_bcast

  IMPLICIT NONE

  INTEGER,          INTENT(IN) :: current_iq,iunpun
  INTEGER,          INTENT(OUT) :: ierr
  INTEGER :: imode0, imode, irr, ipert, iq 
  !
  ierr=0
  IF (ionode) THEN
     CALL iotk_scan_begin( iunpun, "IRREPS_INFO" )
     !
     CALL iotk_scan_dat(iunpun,"QPOINT_NUMBER",iq)
  ENDIF
  CALL mp_bcast( iq,  ionode_id, intra_image_comm )
  IF (iq /= current_iq) CALL errore('read_modes', &
            'problems with current_iq', 1 )

  IF (ionode) THEN

     CALL iotk_scan_dat(iunpun, "QPOINT_GROUP_RANK", nsymq)
     CALL iotk_scan_dat(iunpun, "MINUS_Q_SYM", minus_q)
     CALL iotk_scan_dat(iunpun, "NUMBER_IRR_REP", nirr)
     imode0=0
     DO irr=1,nirr
        CALL iotk_scan_begin( iunpun, "REPRESENTION"// &
                                   TRIM( iotk_index( irr ) ) )
        CALL iotk_scan_dat(iunpun,"NUMBER_OF_PERTURBATIONS", npert(irr))
        DO ipert=1,npert(irr)
           imode=imode0+ipert
           CALL iotk_scan_begin( iunpun, "PERTURBATION"// &
                                   TRIM( iotk_index( ipert ) ) )
           CALL iotk_scan_dat(iunpun,"SYMMETRY_TYPE_CODE", &
                                                      num_rap_mode(imode))
           CALL iotk_scan_dat(iunpun,"SYMMETRY_TYPE", name_rap_mode(imode))
           CALL iotk_scan_dat(iunpun,"DISPLACEMENT_PATTERN",u(:,imode))
           CALL iotk_scan_end( iunpun, "PERTURBATION"// &
                                   TRIM( iotk_index( ipert ) ) )
        ENDDO
        imode0=imode0+npert(irr)
        CALL iotk_scan_end( iunpun, "REPRESENTION"// &
                                   TRIM( iotk_index( irr ) ) )
     ENDDO
     !
     CALL iotk_scan_end( iunpun, "IRREPS_INFO" )
     !
  ENDIF

  CALL mp_bcast( nirr,  ionode_id, intra_image_comm )
  CALL mp_bcast( npert,  ionode_id, intra_image_comm )
  CALL mp_bcast( nsymq,  ionode_id, intra_image_comm )
  CALL mp_bcast( minus_q,  ionode_id, intra_image_comm )
  CALL mp_bcast( u,  ionode_id, intra_image_comm )
  CALL mp_bcast( name_rap_mode,  ionode_id, intra_image_comm )
  CALL mp_bcast( num_rap_mode,  ionode_id, intra_image_comm )

  RETURN
  END SUBROUTINE read_modes
