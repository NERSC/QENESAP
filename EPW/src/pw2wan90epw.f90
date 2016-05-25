  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PP/pw2wannier - Quantum-ESPRESSO group               
  !------------------------------------------------------------------------
  SUBROUTINE pw2wan90epw 
  !------------------------------------------------------------------------
  ! This is the interface to the Wannier90 code: see http://www.wannier.org
  !
  !
  ! 10/2008  Parellel computation of Amn and Mmn 
  ! 12/2008  Added phase setting of overlap matrix elements
  ! 02/2009  works with standard nk1*nk2*nk3 grids
  ! 12/2009  works with USPP 
  !
  ! RM - Nov/Dec 2014
  ! Imported the noncolinear case implemented by xlzhang
  !
  !------------------------------------------------------------------------
  USE io_global,  ONLY : stdout
  USE klist,      ONLY : nkstot
  USE io_files,   ONLY : prefix
  USE epwcom,     ONLY : write_wfn
  USE noncollin_module, ONLY : noncolin
  USE wannier,    ONLY : seedname2, wvfn_formatted, reduce_unk, ispinw, &
                         ikstart, ikstop, iknum
!  use w90_wannierise,    ONLY : seedname
  !
  IMPLICIT NONE
  CHARACTER(LEN=4)   :: spin_component
  CHARACTER(len=256) :: outdir
  !
  !
  outdir = './'
  seedname2 = prefix
  spin_component = 'none'
  wvfn_formatted = .false.
  reduce_unk= .false.
  !
  !
  !
  WRITE(stdout,*)
  SELECT CASE ( TRIM( spin_component ) )
  CASE ( 'up' )
     WRITE(stdout,*) '    Spin CASE ( up )'
     ispinw  = 1
     ikstart = 1
     ikstop  = nkstot/2
     iknum   = nkstot/2
  CASE ( 'down' )
     WRITE(stdout,*) '    Spin CASE ( down )'
     ispinw = 2
     ikstart = nkstot/2 + 1
     ikstop  = nkstot
     iknum   = nkstot/2
  CASE DEFAULT
     IF (noncolin) THEN
        WRITE(stdout,*) '    Spin CASE ( non-collinear )'
     ELSE
        WRITE(stdout,*) '    Spin CASE ( default = unpolarized )'
     ENDIF
     ispinw = 0
     ikstart = 1
     ikstop  = nkstot
     iknum   = nkstot
  END SELECT
  !
  WRITE(stdout,*) 
  WRITE(stdout,*) '    Initializing Wannier90'
  WRITE(stdout,*) 
  CALL setup_nnkp
  CALL ylm_expansion
  CALL compute_amn_para
  CALL compute_mmn_para
  !
  CALL phases_a_m
  !
  CALL write_band
  !
  IF(write_wfn) CALL write_plot
  !
  WRITE(stdout,*)
  WRITE(stdout,*) '    Running Wannier90'
  CALL run_wannier
  !
  CALL lib_dealloc
  !
  END SUBROUTINE pw2wan90epw
!
!-----------------------------------------------------------------------
SUBROUTINE lib_dealloc
  !-----------------------------------------------------------------------
  !
  USE wannier
  !
  implicit none
  IF (ALLOCATED(m_mat) )     DEALLOCATE(m_mat)
  IF (ALLOCATED(u_mat) )     DEALLOCATE(u_mat)
  IF (ALLOCATED(u_mat_opt) ) DEALLOCATE(u_mat_opt)
  IF (ALLOCATED(a_mat) )     DEALLOCATE(a_mat)
  IF (ALLOCATED(eigval) )    DEALLOCATE(eigval)
  IF (ALLOCATED(lwindow) )   DEALLOCATE(lwindow)
  !

END SUBROUTINE lib_dealloc
!
!-----------------------------------------------------------------------
SUBROUTINE setup_nnkp (  )
  !-----------------------------------------------------------------------
  !
#ifdef __PARA
  USE io_global, ONLY : ionode
#endif
  USE io_global, ONLY : stdout, ionode_id
  USE mp_world,  ONLY : world_comm
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps6, tpi
  USE cell_base, ONLY : at, bg, alat
  USE gvect,     ONLY : g, gg
  USE ions_base, ONLY : nat, tau, ityp, atm
  USE mp,        ONLY : mp_bcast
  USE wvfct,     ONLY : nbnd, npwx
  USE wannier,   ONLY : num_nnmax, mp_grid, atcart, atsym, kpb, g_kpb, &
                        center_w, alpha_w, l_w, mr_w, r_w, zaxis,  &
                        xaxis, excluded_band, rlatt, glatt, gf,    &
                        csph, ig_, iknum, seedname2, kpt_latt, nnb, &
                        num_bands, n_wannier, nexband, nnbx, n_proj
  USE noncollin_module, ONLY : noncolin
  USE constants_epw,    ONLY : bohr
  implicit none
  real(DP) :: g_(3), gg_
  integer  :: ik, ib, ig, iw, ia, indexb, type
  real(DP) :: xnorm, znorm, coseno
  integer  :: exclude_bands(nbnd)

! SP: An interface is required because the Wannier routine has optional
!     arguments
  Interface
    Subroutine wannier_setup(seed__name,mp_grid_loc,num_kpts_loc,&
     real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_tot, &
     num_atoms_loc,atom_symbols_loc,atoms_cart_loc, gamma_only_loc,spinors_loc,&
     nntot_loc,nnlist_loc,nncell_loc,num_bands_loc,num_wann_loc, &
     proj_site_loc,proj_l_loc,proj_m_loc,proj_radial_loc,proj_z_loc, &
     proj_x_loc,proj_zona_loc,exclude_bands_loc,proj_s_loc,proj_s_qaxis_loc)

     USE kinds,    ONLY : dp
     USE wannier,  ONLY : num_nnmax

     implicit none

     character(len=*), intent(in) :: seed__name
     integer, dimension(3), intent(in) :: mp_grid_loc
     integer, intent(in) :: num_kpts_loc
     real(kind=dp), dimension(3,3), intent(in) :: real_lattice_loc
     real(kind=dp), dimension(3,3), intent(in) :: recip_lattice_loc
     real(kind=dp), dimension(3,num_kpts_loc), intent(in) :: kpt_latt_loc
     integer, intent(in) :: num_bands_tot
     integer, intent(in) :: num_atoms_loc
     character(len=*), dimension(num_atoms_loc), intent(in) :: atom_symbols_loc
     real(kind=dp), dimension(3,num_atoms_loc), intent(in) :: atoms_cart_loc
     logical, intent(in) :: gamma_only_loc
     logical, intent(in) :: spinors_loc
     integer, intent(out) :: nntot_loc
     integer, dimension(num_kpts_loc,num_nnmax), intent(out) :: nnlist_loc
     integer,dimension(3,num_kpts_loc,num_nnmax), intent(out) :: nncell_loc
     integer, intent(out) :: num_bands_loc
     integer, intent(out) :: num_wann_loc
     real(kind=dp), dimension(3,num_bands_tot), intent(out) :: proj_site_loc
     integer, dimension(num_bands_tot), intent(out) :: proj_l_loc
     integer, dimension(num_bands_tot), intent(out) :: proj_m_loc
     integer, dimension(num_bands_tot), intent(out) :: proj_radial_loc
     real(kind=dp), dimension(3,num_bands_tot), intent(out) :: proj_z_loc
     real(kind=dp), dimension(3,num_bands_tot), intent(out) :: proj_x_loc
     real(kind=dp), dimension(num_bands_tot), intent(out) :: proj_zona_loc
     integer, dimension(num_bands_tot), intent(out) :: exclude_bands_loc
     integer, dimension(num_bands_tot), optional, intent(out) :: proj_s_loc
     real(kind=dp), dimension(3,num_bands_tot), optional, intent(out) :: proj_s_qaxis_loc

    End Subroutine wannier_setup
  End Interface

  num_nnmax = 32

  ! aam: translations between PW2Wannier90 and Wannier90
  ! pw2wannier90   <==>   Wannier90
  !    nbnd                num_bands_tot
  !    n_wannier           num_wann
  !    num_bands           num_bands
  !    nat                 num_atoms
  !    iknum               num_kpts
  !    rlatt               transpose(real_lattice)
  !    glatt               transpose(recip_lattice)
  !    kpt_latt            kpt_latt
  !    nnb                 nntot
  !    kpb                 nnlist
  !    g_kpb               nncell
  !    mp_grid             mp_grid
  !    center_w            proj_site
  !    l_w,mr_w,r_w        proj_l,proj_m,proj_radial
  !    xaxis,zaxis         proj_x,proj_z
  !    alpha_w             proj_zona
  !    exclude_bands       exclude_bands
  !    atcart              atoms_cart
  !    atsym               atom_symbols

  ALLOCATE( atcart(3,nat), atsym(nat) )
  ALLOCATE( kpb(iknum,num_nnmax), g_kpb(3,iknum,num_nnmax) )
  ALLOCATE( center_w(3,nbnd), alpha_w(nbnd), l_w(nbnd), &
       mr_w(nbnd), r_w(nbnd), zaxis(3,nbnd), xaxis(3,nbnd) )
  ALLOCATE( excluded_band(nbnd) )

  ! real lattice (Cartesians, Angstrom)
  rlatt(:,:) = transpose(at(:,:))*alat*bohr
  ! reciprocal lattice (Cartesians, Angstrom)
  glatt(:,:) = transpose(bg(:,:))*tpi/(alat*bohr)
  ! atom co-ordinates in Cartesian co-ords and Angstrom units
  atcart(:,:) = tau(:,:)*bohr*alat
  ! atom symbols
  DO ia=1,nat
     type=ityp(ia)
     atsym(ia)=atm(type)
  ENDDO

#ifdef __PARA
   IF (ionode) THEN
#endif
      CALL wannier_setup(seedname2, mp_grid, iknum,       &  ! input
           rlatt, glatt, kpt_latt, nbnd,                 &  ! input
           nat, atsym, atcart, .false., noncolin,        &  ! input
           nnb, kpb, g_kpb, num_bands, n_wannier,        &  ! output
           center_w, l_w, mr_w, r_w, zaxis,              &  ! output
           xaxis, alpha_w, exclude_bands)                   ! output
#ifdef __PARA
   ENDIF
#endif
   
   CALL mp_bcast(nnb,          ionode_id, world_comm )
   CALL mp_bcast(kpb,          ionode_id, world_comm )
   CALL mp_bcast(g_kpb,        ionode_id, world_comm )
   CALL mp_bcast(num_bands,    ionode_id, world_comm )
   CALL mp_bcast(n_wannier,    ionode_id, world_comm )
   CALL mp_bcast(center_w,     ionode_id, world_comm )
   CALL mp_bcast(l_w,          ionode_id, world_comm )
   CALL mp_bcast(mr_w,         ionode_id, world_comm )
   CALL mp_bcast(r_w,          ionode_id, world_comm )
   CALL mp_bcast(zaxis,        ionode_id, world_comm )
   CALL mp_bcast(xaxis,        ionode_id, world_comm )
   CALL mp_bcast(alpha_w,      ionode_id, world_comm )
   CALL mp_bcast(exclude_bands,ionode_id, world_comm )
   CALL mp_bcast(noncolin,     ionode_id, world_comm )
   !
   !
   ! n_proj = nr. of projections (=#WF unless spinors then =#WF/2) 
   IF (noncolin) THEN
      n_proj=n_wannier/2
   ELSE
      n_proj=n_wannier
   ENDIF
   !
   WRITE (stdout,*)
   WRITE (stdout,*) '    Initial Wannier projections'
   WRITE (stdout,*)
   ! RM changed according to QE4.0.3/PP/pw2wannier90 
   DO iw=1,n_proj
   !DO iw=1,n_wannier
      WRITE (stdout, '(5x,"(",3f10.5,") :  l = ",i3, " mr = ", i3)') center_w(:,iw), l_w(iw), mr_w(iw)
   ENDDO
   !
   WRITE(stdout,'(/,"      - Number of bands is (",i3,")")') num_bands 
   WRITE(stdout,'("      - Number of wannier functions is (",i3,")")') n_wannier 
   !
   ! RM changed according to QE4.0.3/PP/pw2wannier90
   ALLOCATE( gf(npwx,n_proj), csph(16,n_proj) )
   !ALLOCATE( gf(npwx,n_wannier), csph(16,n_wannier) ) 
   !
   excluded_band(1:nbnd)=.false.
   nexband=0
   band_loop: DO ib=1,nbnd
      indexb=exclude_bands(ib)
      IF (indexb>nbnd .or. indexb<0) THEN
         CALL errore('setup_nnkp',' wrong excluded band index ', 1)
      ELSEIF (indexb.eq.0) THEN
         exit band_loop
      ELSE
         nexband=nexband+1
         excluded_band(indexb)=.true.
      ENDIF
   ENDDO band_loop
 
   IF ( (nbnd-nexband).ne.num_bands ) &
       CALL errore('setup_nnkp',' something wrong with num_bands',1)
 
   ! RM changed according to QE4.0.3/PP/pw2wannier90 
   DO iw=1,n_proj
   !DO iw=1,n_wannier
      xnorm = sqrt(xaxis(1,iw)*xaxis(1,iw) + xaxis(2,iw)*xaxis(2,iw) + &
           xaxis(3,iw)*xaxis(3,iw))
      IF (xnorm < eps6) CALL errore ('setup_nnkp',' |xaxis| < eps ',1)
      znorm = sqrt(zaxis(1,iw)*zaxis(1,iw) + zaxis(2,iw)*zaxis(2,iw) + &
           zaxis(3,iw)*zaxis(3,iw))
      IF (znorm < eps6) CALL errore ('setup_nnkp',' |zaxis| < eps ',1)
      coseno = (xaxis(1,iw)*zaxis(1,iw) + xaxis(2,iw)*zaxis(2,iw) + &
           xaxis(3,iw)*zaxis(3,iw))/xnorm/znorm
      IF (abs(coseno) > eps6) &
           CALL errore('setup_nnkp',' xaxis and zaxis are not orthogonal !',1)
      IF (alpha_w(iw) < eps6) &
           CALL errore('setup_nnkp',' zona value must be positive', 1)
      ! convert wannier center in cartesian coordinates (in unit of alat)
      CALL cryst_to_cart( 1, center_w(:,iw), at, 1 )
   ENDDO
   WRITE(stdout,*) '     - All guiding functions are given '

   nnbx=0
   nnb=max(nnbx,nnb)
 
   ALLOCATE( ig_(iknum,nnb) )
 
   DO ik=1, iknum
      DO ib = 1, nnb
         g_(:) = REAL( g_kpb(:,ik,ib) )
         CALL cryst_to_cart (1, g_, bg, 1)
         gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
         ig_(ik,ib) = 0
         ig = 1
         DO WHILE  (gg(ig) <= gg_ + eps6)
            IF ( (abs(g(1,ig)-g_(1)) < eps6) .and.  &
                 (abs(g(2,ig)-g_(2)) < eps6) .and.  &
                 (abs(g(3,ig)-g_(3)) < eps6)  ) ig_(ik,ib) = ig
            ig= ig +1
         ENDDO
      ENDDO
   ENDDO
   WRITE(stdout,*) '     - All neighbours are found '
   WRITE(stdout,*)
   !
   RETURN
   !
END SUBROUTINE setup_nnkp
 !
 !-----------------------------------------------------------------------
SUBROUTINE run_wannier
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, ionode_id
#ifdef __PARA
  USE io_global, ONLY : ionode
#endif
!  USE io_epw,    ONLY : iummn
  USE ions_base, ONLY : nat
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  USE cell_base, ONLY : celldm
  USE io_files,  ONLY : prefix
  USE wannier
  USE epwcom,    ONLY : eig_read
  USE wvfct,     ONLY : nbnd
  USE constants_epw, ONLY : czero, bohr

  implicit none

  integer             :: i, ik, ibnd, dummy1, dummy2, ios
  character (len=256) :: tempfile
  !
  ALLOCATE(u_mat(n_wannier,n_wannier,iknum))
  ALLOCATE(u_mat_opt(num_bands,n_wannier,iknum))
  ALLOCATE(lwindow(num_bands,iknum))
  ALLOCATE(wann_centers(3,n_wannier))
  ALLOCATE(wann_spreads(n_wannier))
  !
  u_mat_opt = czero
  !
#ifdef __PARA
   IF (ionode) THEN
#endif
      ! read in external eigenvalues, e.g.  GW
      IF (eig_read) then
         WRITE (stdout,'(5x,a,i5,a,i5,a)') "Reading electronic eigenvalues (", &
              nbnd, ",", iknum,")"
         tempfile=trim(prefix)//'.eig'
         OPEN(1, file=tempfile, form='formatted', action='read', iostat=ios)
         IF (ios /= 0) CALL errore ('run_wannier','error opening' // tempfile, 1)
         !
         ! the form should be band, kpt, eigenvalue
         !
         DO ik = 1, iknum
            DO ibnd = 1, nbnd
               READ (1,*) dummy1, dummy2, eigval (ibnd,ik)
               IF (dummy1.ne.ibnd) CALL errore('run_wannier', "Incorrect eigenvalue file", 1)
               IF (dummy2.ne.ik) CALL errore('run_wannier', "Incorrect eigenvalue file", 1)
            ENDDO
         ENDDO
         CLOSE(1)
      ENDIF
      !
! SP : This file is not used for now. Only required to build the UNK file
!      tempfile=trim(prefix)//'.mmn'
!      OPEN(iummn, file=tempfile, iostat=ios, form='unformatted')
!      WRITE(iummn) m_mat
!      CLOSE(iummn)

      CALL wannier_run(seedname2, mp_grid, iknum,    &                 ! input
           rlatt, glatt, kpt_latt, num_bands,       &                 ! input
           n_wannier, nnb, nat, atsym,              &                 ! input
           atcart, .false., m_mat, a_mat, eigval,   &                 ! input
           u_mat, u_mat_opt, lwindow, wann_centers, &                 ! output
           wann_spreads, spreads)                                     ! output

#ifdef __PARA
   ENDIF
#endif
  !
  CALL mp_bcast(u_mat,       ionode_id, world_comm )
  CALL mp_bcast(u_mat_opt,   ionode_id, world_comm )
  CALL mp_bcast(lwindow,     ionode_id, world_comm )
  CALL mp_bcast(wann_centers,ionode_id, world_comm )
  CALL mp_bcast(wann_spreads,ionode_id, world_comm )
  CALL mp_bcast(spreads,     ionode_id, world_comm )
  !
  !
  ! output the results of the wannierization
  !
  WRITE (stdout,*)
  WRITE (stdout,*) '    Wannier Function centers (cartesian, alat) and spreads (ang):'
  WRITE (stdout,*)
  ! RM - loop is up to n_wannier according to W90/wannierise.F90
  DO i=1,n_wannier
     WRITE (stdout, '(5x,"(",3f10.5,") :  ",f8.5)') wann_centers(:,i)/celldm(1)/bohr, wann_spreads(i)
  ENDDO
  WRITE (stdout,*)
  !
  ! store the final minimisation matrix on disk for later use
  !
  CALL write_filukk
  !
  RETURN
  !
END SUBROUTINE run_wannier
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE compute_amn_para
!-----------------------------------------------------------------------
!  adapted from compute_amn in pw2wannier90.f90
!  parallelization on k-points has been added
!  10/2008 Jesse Noffsinger UC Berkeley
!
   USE io_global,       ONLY : stdout 
   USE kinds,           ONLY : DP
   USE klist,           ONLY : xk, nks
   USE wvfct,           ONLY : nbnd, npw, npwx, igk, g2kin
   USE gvecw,           ONLY : ecutwfc
   USE wavefunctions_module,  ONLY : evc
   USE units_ph,        ONLY : lrwfc, iuwfc
   USE gvect,           ONLY : g, ngm
   USE cell_base,       ONLY : tpiba2
   USE uspp,            ONLY : nkb, vkb
   USE becmod,          ONLY : becp, calbec, deallocate_bec_type, allocate_bec_type
   USE wannier,         ONLY : csph, excluded_band, gf, num_bands, &
                               n_wannier, iknum, n_proj, a_mat 
   USE uspp_param,      ONLY : upf
   USE noncollin_module,ONLY : noncolin, npol
   USE constants_epw,   ONLY : czero
#ifdef __NAG
   USE f90_unix_io,    ONLY : flush
#endif
#ifdef __PARA
   USE mp_global,       ONLY : npool, intra_pool_comm, inter_pool_comm
   USE mp,              ONLY : mp_sum
#endif

   implicit none
  
   complex(DP) :: amn, ZDOTC
   complex(DP), allocatable :: sgf(:,:)
   integer :: amn_tot, ik, ibnd, ibnd1, iw, nkq, nkq_abs, ipool, ik_g, ipol, &
              istart
   logical            :: any_uspp
   real(kind=DP)      :: zero_vect(3)
   !
   !nocolin: we have half as many projections g(r) defined as wannier
   !         functions. We project onto (1,0) (ie up spin) and then onto
   !         (0,1) to obtain num_wann projections. jry
   !
   any_uspp = ANY( upf(:)%tvanp )
   !
   IF (any_uspp .and. noncolin) CALL errore('pw2wan90epw',&
       'noncolin calculation not implimented with USP',1)
   !
   ALLOCATE( a_mat(num_bands,n_wannier,iknum))
   ! RM - changed according to QE4.0.3/PP/pw2wannie90.f90
   ALLOCATE( sgf(npwx,n_proj) )
   !ALLOCATE( sgf(npwx,n_wannier) )
   !
   !
   ! initialize
   a_mat = czero
   zero_vect = 0.d0
   !
   amn_tot = iknum * nbnd * n_wannier
   !
   WRITE (stdout,'(5x,a)') 'AMN'
   !
   IF (any_uspp) then
      CALL deallocate_bec_type ( becp )
      CALL allocate_bec_type ( nkb, n_wannier, becp )
      CALL init_us_1
   ENDIF
   !
#ifdef __PARA
   WRITE(stdout,'(6x,a,i5,a,i4,a)') 'k points = ',iknum, ' in ', npool, ' pools'
#endif
   DO ik=1,nks
      CALL ktokpmq ( xk(:,ik), zero_vect, +1, ipool, nkq, nkq_abs )
      ik_g = nkq_abs
      !
      WRITE (stdout,'(5x,i8, " of ", i4,a)') ik , nks, ' on ionode'
      CALL flush(6)
      CALL davcio( evc, lrwfc, iuwfc, ik, -1 )
      !
      CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      CALL generate_guiding_functions(ik)   ! they are called gf(npw,n_proj)
      !
      !  USPP
      !
      IF (any_uspp) THEN
         CALL init_us_2 (npw, igk, xk (1, ik), vkb)
         ! below we compute the product of beta functions with trial func.
         ! RM - changed according to QE4.0.3/PP/pw2wannie90.f90
         CALL calbec ( npw, vkb, gf, becp, n_proj )
         !CALL calbec(npw, vkb, gf, becp) 
         ! and we use it for the product S|trial_func>
         ! RM - changed according to QE4.0.3/PP/pw2wannie90.f90
         CALL s_psi (npwx, npw, n_proj, gf, sgf)  
         !CALL s_psi (npwx, npw, n_wannier, gf, sgf)
      ELSE
         sgf(:,:) = gf(:,:)
      ENDIF
      !
      ! RM changed according to QE4.0.3/PP/pw2wannier90 
      IF (noncolin) THEN
         ! we do the projection as g(r)*a(r) and g(r)*b(r)
         DO ipol=1,npol
            istart = (ipol-1)*npwx + 1
            DO iw = 1,n_proj
               ibnd1 = 0
               DO ibnd = 1,nbnd
                  IF (excluded_band(ibnd)) CYCLE
                  amn = (0.0d0,0.0d0)
                  amn = ZDOTC(npw,evc(istart,ibnd),1,sgf(1,iw),1)
#ifdef __PARA
                  CALL mp_sum(amn, intra_pool_comm)
#endif
                  ibnd1=ibnd1+1
                  a_mat(ibnd1,iw+n_proj*(ipol-1),ik_g) = amn
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO iw = 1,n_proj
            ibnd1 = 0
            DO ibnd = 1,nbnd
               IF (excluded_band(ibnd)) CYCLE
               amn=(0.0_dp,0.0_dp)
               amn = ZDOTC(npw,evc(1,ibnd),1,sgf(1,iw),1)
#ifdef __PARA
               CALL mp_sum(amn, intra_pool_comm)
#endif
               ibnd1=ibnd1+1
               a_mat(ibnd1,iw,ik_g) = amn
            ENDDO !bands
         ENDDO !wannier fns
      ENDIF
   ENDDO  ! k-points
   !
   DEALLOCATE (sgf,csph)
   IF (any_uspp) CALL deallocate_bec_type ( becp )
   !
#ifdef __PARA
   CALL mp_sum(a_mat, inter_pool_comm)
#endif
   !
   !
   !
   WRITE(stdout,*)
   WRITE(stdout,'(5x,a)') 'AMN calculated'
   !
   RETURN
!-----------------------------------------------------------------------
END SUBROUTINE compute_amn_para
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
SUBROUTINE compute_mmn_para
!-----------------------------------------------------------------------
!
!  adapted from compute_mmn in pw2wannier90.f90
!  parallelization on k-points has been added
!  10/2008 Jesse Noffsinger UC Berkeley
!
   USE io_global,       ONLY : stdout, ionode
   USE io_files,        ONLY : diropn
   USE mp_global,       ONLY : my_pool_id
#ifdef __PARA
   USE mp_global,       ONLY : npool, intra_pool_comm
#endif
   USE kinds,           ONLY : DP
   USE wvfct,           ONLY : nbnd, npw, npwx, igk, g2kin
   USE gvecw,           ONLY : ecutwfc
   USE control_flags,   ONLY : gamma_only
   USE wavefunctions_module, ONLY : evc, psic, psic_nc
   USE units_ph,        ONLY : lrwfc, iuwfc
   USE gvecs,           ONLY : nls, nlsm
   USE fft_base,        ONLY : dffts
   USE fft_interfaces,  ONLY : fwfft, invfft
   USE klist,           ONLY : nkstot, xk, nks
   USE gvect,           ONLY : g, ngm, gstart
   USE cell_base,       ONLY : tpiba2, omega, tpiba, bg
   USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
   USE uspp,            ONLY : nkb, vkb
   USE uspp_param,      ONLY : upf, lmaxq, nh
   USE becmod,          ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type
   USE noncollin_module,ONLY : noncolin, npol
   USE wannier,         ONLY : m_mat, num_bands, nnb, iknum, g_kpb, kpb, ig_, &
                               excluded_band
   USE constants_epw,   ONLY : czero, twopi
#ifdef __NAG
   USE f90_unix_io,     ONLY : flush
#endif
#ifdef __PARA
   USE mp,              ONLY : mp_sum
   USE mp_global,       ONLY : inter_pool_comm
#endif

   implicit none

   integer :: mmn_tot ,ik, ikp, ib, npwq, i, m, n
   integer :: ikb, jkb, ih, jh, na, nt, ijkb0, ind, nbt
   complex(DP), allocatable :: phase(:), aux(:), aux2(:),  evcq(:,:), & 
                               Mkb(:,:), aux_nc(:,:),becp2(:,:)
   complex(DP), allocatable :: qb(:,:,:,:), qgm(:)
   real(DP), allocatable    :: qg(:), ylm(:,:), dxk(:,:)
   real(DP), ALLOCATABLE    :: rbecp2(:,:)
   integer, allocatable     :: igkq(:)
   real(DP)                 :: xktot(3,nkstot)
   real(DP), dimension(3)   :: zero_vect
   complex(DP)              :: mmn, ZDOTC, phase1
   real(DP)                 :: arg, g_(3)
   logical                  :: any_uspp, exst
   integer                  :: nkq, nkq_abs, ipool, ipol, istart, iend
   integer                  :: ik_g, ikp_g, ind0 !, iummn
   ! 
   any_uspp = ANY( upf(:)%tvanp )
   !
   IF (any_uspp .and. noncolin) CALL errore('pw2wan90epw',&
       'noncolin calculation not implimented with USP',1)
   !
   ALLOCATE( phase(dffts%nnr) ) 
   ALLOCATE( igkq(npwx) )
   ! RM changed according to QE4.0.3/PP/pw2wannie90.f90
   ALLOCATE( evcq(npol*npwx,nbnd) )
   ! ALLOCATE( evcq(npwx,nbnd) )
   !
   IF (noncolin) THEN
      ALLOCATE( aux_nc(npwx,npol) )
   ELSE
      ALLOCATE( aux(npwx) )
   ENDIF

   IF (gamma_only) ALLOCATE(aux2(npwx))
   !
   ALLOCATE(m_mat(num_bands,num_bands,nnb,iknum))
   !
   !
   ! close all the wfc files to allow access for each pool to all wfs
   CLOSE (unit = iuwfc,  status = 'keep')
   !
   WRITE (stdout,*)
   WRITE (stdout,'(5x,a)') 'MMN'
   !
   ! Get all the k-vector coords to each pool via xktot
   xktot = 0.d0
   IF (ionode) then
      DO i = 1,nkstot
         xktot(:,i) = xk(:,i)
      ENDDO
   ENDIF
#ifdef __PARA
   CALL mp_sum(xktot, inter_pool_comm)
#endif
   !
   zero_vect = 0.0d0
   m_mat = czero
   !
   mmn_tot = iknum * nnb * nbnd * nbnd
   !
   !   USPP
   !
   !
   IF(any_uspp) THEN
      CALL init_us_1
      CALL allocate_bec_type ( nkb, nbnd, becp )
      IF (gamma_only) THEN
         ALLOCATE ( rbecp2(nkb,nbnd))
      ELSE
         ALLOCATE ( becp2(nkb,nbnd) )
      ENDIF
   ENDIF
   !
   nbt = nnb * iknum
   !
   ALLOCATE( qg(nbt) )
   ALLOCATE (dxk(3,nbt))
   !
   ind = 0
   DO ik=1,iknum
      DO ib=1,nnb
         ind = ind + 1
         ikp = kpb(ik,ib) 
         !
         g_(:) = REAL( g_kpb(:,ik,ib) )
         CALL cryst_to_cart (1, g_, bg, 1)
         dxk(:,ind) = xktot(:,ikp) +g_(:) - xktot(:,ik) 
         qg(ind) = dxk(1,ind)*dxk(1,ind)+dxk(2,ind)*dxk(2,ind)+dxk(3,ind)*dxk(3,ind)
      ENDDO
  ENDDO
  !
  !  USPP
  !
  IF (any_uspp) THEN
     !
     ALLOCATE( ylm(nbt,lmaxq*lmaxq), qgm(nbt) )
     ALLOCATE( qb (nkb, nkb, ntyp, nbt) )
     !
     CALL ylmr2 (lmaxq*lmaxq, nbt, dxk, qg, ylm)
     qg(:) = sqrt(qg(:)) * tpiba
     !
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           DO ih = 1, nh (nt)
              DO jh = 1, nh (nt)
                 CALL qvan2 (nbt, ih, jh, nt, qg, qgm, ylm)
                 qb (ih, jh, nt, 1:nbt) = omega * qgm(1:nbt)
              ENDDO
           ENDDO
        ENDIF
     ENDDO
     !
     DEALLOCATE (qg, qgm, ylm )
     !
  ENDIF
  !
  !
  
  ALLOCATE( Mkb(nbnd,nbnd) )
  !
#ifdef __PARA
  WRITE(stdout,'(6x,a,i5,a,i4,a)') 'k points = ',iknum, ' in ', npool, ' pools'
#endif
  ! get the first k-point in this pool 
  CALL ktokpmq( xk(:, 1), zero_vect, +1, ipool, nkq, nkq_abs )
  ind0 = (nkq_abs-1) * nnb
  !
  ind = ind0
  DO ik=1,nks 
     CALL ktokpmq( xk(:, ik), zero_vect, +1, ipool, nkq, nkq_abs)
     ik_g = nkq_abs
     !
     WRITE (stdout,'(5x,i8, " of ", i4,a)') ik , nks, ' on ionode'
     CALL flush(6)
     !
     !
     !
     CALL readwfc(my_pool_id+1, ik, evc)
     !
     CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     !
     !  USPP
     !
     IF (any_uspp) THEN
        CALL init_us_2 (npw, igk, xk(1,ik), vkb)
        CALL calbec(npw, vkb, evc, becp) 
     ENDIF
     !
     !
     DO ib=1,nnb  ! loop on finite diff vectors
        ind = ind + 1
        !
        ikp = kpb(ik_g,ib)
        ikp_g = ikp 
        !
        CALL ktokpmq( xk(:, ik), xktot(:,ikp_g)-xk(:,ik), +1, ipool, nkq, nkq_abs )
        !
        ! read wfc at k+b
        !
        CALL readwfc( ipool, nkq, evcq )
        !
        CALL gk_sort (xktot(1,ikp_g), ngm, g, ecutwfc / tpiba2, npwq, igkq, g2kin)
        !
        ! compute the phase
        !
        phase(:) = (0.d0,0.d0)
        IF ( ig_(ik_g,ib)>0) phase( nls(ig_(ik_g,ib)) ) = (1.d0,0.d0)
        CALL invfft ('Wave', phase, dffts)
        !
        !  USPP
        !
        IF(any_uspp) THEN
           CALL init_us_2 (npwq, igkq, xk(1,ikp), vkb)
           ! below we compute the product of beta functions with |psi>
           IF (gamma_only) THEN
              CALL calbec ( npwq, vkb, evcq, rbecp2 )
           ELSE
              CALL calbec ( npwq, vkb, evcq, becp2 )
           ENDIF
        ENDIF
        !
        !
        Mkb(:,:) = (0.0d0,0.0d0) 
        !
        IF (any_uspp) THEN
           ijkb0 = 0
           DO nt = 1, ntyp
              IF ( upf(nt)%tvanp ) THEN
                 DO na = 1, nat
                    !
                    arg = DOT_PRODUCT( dxk(:,ind), tau(:,na) ) * twopi
                    phase1 = CMPLX ( COS(arg), -SIN(arg) )
                    !
                    IF ( ityp(na) == nt ) THEN
                       DO jh = 1, nh(nt)
                          jkb = ijkb0 + jh
                          DO ih = 1, nh(nt)
                             ikb = ijkb0 + ih
                             !
                             DO m = 1,nbnd
                               IF (excluded_band(m)) CYCLE
                               IF (gamma_only) THEN
                                 DO n=1,m ! Mkb(m,n) is symmetric in m and n for gamma_only case
                                    IF (excluded_band(n)) CYCLE
                                    Mkb(m,n) = Mkb(m,n) + &
                                         phase1 * qb(ih,jh,nt,ind) * &
                                         becp%r(ikb,m) * rbecp2(jkb,n)
                                 ENDDO
                               ELSE
                                 DO n=1,nbnd
                                    IF (excluded_band(n)) CYCLE
                                    Mkb(m,n) = Mkb(m,n) + &
                                         phase1 * qb(ih,jh,nt,ind) * &
                                         conjg( becp%k(ikb,m) ) *becp2(jkb,n)

                                 ENDDO
                               ENDIF
                             ENDDO ! m
                          ENDDO !ih
                       ENDDO !jh
                       ijkb0 = ijkb0 + nh(nt)
                    ENDIF  !ityp
                 ENDDO  !nat 
              ELSE  !tvanp
                 DO na = 1, nat
                    IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                 ENDDO
              ENDIF !tvanp
           ENDDO !ntyp
        ENDIF ! any_uspp
        !
        !
        DO m=1,nbnd  ! loop on band m
           IF (excluded_band(m)) CYCLE
           !
           IF (noncolin) THEN
              psic_nc(:,:) = (0.d0, 0.d0)
              DO ipol=1,2!npol
                 istart=(ipol-1)*npwx+1
                 iend=istart+npw-1
                 psic_nc(nls (igk (1:npw) ),ipol ) = evc(istart:iend, m)
                 CALL invfft ('Wave', psic_nc(:,ipol), dffts)
                 psic_nc(1:dffts%nnr,ipol) = psic_nc(1:dffts%nnr,ipol) * &
                                               phase(1:dffts%nnr)
                 CALL fwfft ('Wave', psic_nc(:,ipol), dffts)
                 aux_nc(1:npwq,ipol) = psic_nc(nls (igkq(1:npwq) ),ipol )
              ENDDO
           ELSE
              psic(:) = (0.d0, 0.d0)
              psic(nls (igk (1:npw) ) ) = evc (1:npw, m)
              IF(gamma_only) psic(nlsm(igk (1:npw) ) ) = conjg(evc (1:npw, m))
              CALL invfft ('Wave', psic, dffts)
              psic(1:dffts%nnr) = psic(1:dffts%nnr) * phase(1:dffts%nnr)
              CALL fwfft ('Wave', psic, dffts)
              aux(1:npwq)  = psic(nls (igkq(1:npwq) ) )
           ENDIF
           IF(gamma_only) THEN
              IF (gstart==2) psic(nlsm(1)) = (0.d0,0.d0)
              aux2(1:npwq) = conjg(psic(nlsm(igkq(1:npwq) ) ) )
           ENDIF
         !  aa = 0.d0
           !
           !  Mkb(m,n) = Mkb(m,n) + \sum_{ijI} qb_{ij}^I * e^-i(b*tau_I)
           !             <psi_m,k1| beta_i,k1 > < beta_j,k2 | psi_n,k2 > 
           !
           IF (gamma_only) THEN
             DO n=1,m ! Mkb(m,n) is symmetric in m and n for gamma_only case
                IF (excluded_band(n)) CYCLE
                mmn = zdotc (npwq, aux,1,evcq(1,n),1) &
                     + conjg(zdotc(npwq,aux2,1,evcq(1,n),1))
#ifdef __PARA
                CALL mp_sum(mmn, intra_pool_comm)
#endif
                Mkb(m,n) = mmn + Mkb(m,n)
                IF (m/=n) Mkb(n,m) = Mkb(m,n) ! fill other half of matrix by symmetry
             ENDDO
           ELSEIF (noncolin) THEN
              DO n=1,nbnd
                 IF (excluded_band(n)) CYCLE
                 mmn=(0.d0, 0.d0)
!                 DO ipol=1,2
!                    mmn = mmn+ZDOTC (npwq, aux_nc(1,ipol),1,evcq_nc(1,ipol,n),1)
                 mmn = mmn + ZDOTC (npwq, aux_nc(1,1),1,evcq(1,n),1) &
                       + ZDOTC (npwq, aux_nc(1,2),1,evcq(npwx+1,n),1)
!                 ENDDO
#ifdef __PARA
                 CALL mp_sum(mmn, intra_pool_comm)
#endif
                 Mkb(m,n) = mmn + Mkb(m,n)
               !  aa = aa + abs(mmn)**2
              ENDDO
           ELSE
              DO n=1,nbnd
                 IF (excluded_band(n)) CYCLE
                 mmn = ZDOTC (npwq, aux,1,evcq(1,n),1)
#ifdef __PARA
                 CALL mp_sum(mmn, intra_pool_comm)
#endif
                 Mkb(m,n) = mmn + Mkb(m,n)
!                 aa = aa + abs(mmn)**2
              ENDDO
           ENDIF
        ENDDO   ! m
        !
        DO n=1,nbnd
           IF (excluded_band(n)) CYCLE
           DO m=1,nbnd
              IF (excluded_band(m)) CYCLE
              m_mat(m,n,ib,ik_g)=Mkb(m,n)
           ENDDO
        ENDDO
        !
        !
     ENDDO !ib
  ENDDO !ik
  !
#ifdef __PARA
  CALL mp_sum(m_mat, inter_pool_comm)
#endif
  !
  !
  IF (gamma_only) DEALLOCATE(aux2)
  DEALLOCATE (Mkb, dxk, phase, evcq, igkq)
  IF (noncolin) THEN
     DEALLOCATE(aux_nc)
  ELSE
     DEALLOCATE(aux)
  ENDIF
  IF(any_uspp) THEN
     DEALLOCATE (  qb)
     CALL deallocate_bec_type (becp)
     IF (gamma_only) THEN
         DEALLOCATE (rbecp2)
      ELSE
         DEALLOCATE (becp2)
      ENDIF
   ENDIF

  !
  !
  WRITE(stdout,'(5x,a)') 'MMN calculated'
  !
  ! reopen wfc here, leaving unit=20 in the same state
  iuwfc = 20
  CALL diropn(iuwfc,'wfc',lrwfc,exst)  
  !
  RETURN
END SUBROUTINE compute_mmn_para
!
!------------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE compute_pmn_para
!-----------------------------------------------------------------------
!  adapted from compute_amn_para
!  06/2010  Jesse Noffsinger
!  
!
   USE io_global,       ONLY : stdout
#ifdef __PARA
   USE mp_global,       ONLY : intra_pool_comm
   USE mp,              ONLY : mp_sum
#endif
   USE kinds,           ONLY : DP
   USE klist,           ONLY : xk, nks
   USE wvfct,           ONLY : nbnd, npw, npwx, igk, g2kin
   USE gvecw,           ONLY : ecutwfc
   USE wavefunctions_module,  ONLY : evc
   USE units_ph,        ONLY : lrwfc, iuwfc
   USE gvect,           ONLY : g, ngm
   USE cell_base,       ONLY : tpiba2, tpiba
   USE noncollin_module,ONLY : noncolin
   USE elph2,           ONLY : dmec
   USE constants_epw,   ONLY : czero
   implicit none
  
   integer :: ik, ibnd, ig, jbnd
   COMPLEX(DP) :: dipole_aux(3,nbnd,nbnd), caux 
   !
   !
   ALLOCATE( dmec(3,nbnd,nbnd,nks) )
   !
   ! initialize
   dmec = czero
   dipole_aux(:,:,:) = (0.0d0,0.0d0)
   !
   !
   DO ik=1,nks
      !
      ! read wfc for the given kpt
      CALL davcio( evc, lrwfc, iuwfc, ik, -1 )
      !
      ! setup k+G grids for each kpt
      CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      !
      dipole_aux = (0.d0, 0.d0)
      DO jbnd = 1,nbnd
         DO ibnd = 1,nbnd
            !
            IF ( ibnd .eq. jbnd ) CYCLE
            !
            ! taken from PP/epsilon.f90 SUBROUTINE dipole_calc
            DO ig = 1, npw
               IF (igk(ig) .gt. SIZE(g,2) .or. igk(ig).lt.1) CYCLE
               !
               caux = conjg(evc(ig,ibnd))*evc(ig,jbnd) 
               !
               dipole_aux(:,ibnd,jbnd) = dipole_aux(:,ibnd,jbnd) + &
                                       ( g(:,igk(ig)) ) * caux
               !
               ! RM - this should cover the noncolin case
               IF (noncolin) THEN
                  !
                  caux = conjg(evc(ig+npwx,ibnd))*evc(ig+npwx,jbnd)
                  !
                  dipole_aux(:,ibnd,jbnd) = dipole_aux(:,ibnd,jbnd) + &
                                          ( g(:,igk(ig)) ) * caux
                  !
               ENDIF
               !
            ENDDO
            !
         ENDDO !bands i
      ENDDO ! bands j
      ! metal diagonal part
      DO ibnd = 1, nbnd
         DO ig = 1, npw
            IF (igk(ig) .gt. SIZE(g,2) .or. igk(ig).lt.1) CYCLE
            !
            caux = conjg(evc(ig,ibnd))*evc(ig,ibnd) 
            !
            dipole_aux(:,ibnd,ibnd) = dipole_aux(:,ibnd,ibnd) + &
                                    ( g(:,igk(ig)) + xk(:,ik) ) * caux
            !
            ! RM - this should cover the noncolin case
            IF (noncolin) THEN
               !
               caux = conjg(evc(ig+npwx,ibnd))*evc(ig+npwx,ibnd)
               !
               dipole_aux(:,ibnd,ibnd) = dipole_aux(:,ibnd,ibnd) + &
                                       ( g(:,igk(ig)) + xk(:,ik) ) * caux
               !
            ENDIF
            !
         ENDDO
      ENDDO
      ! need to divide by 2pi/a to fix the units
      dmec(:,:,:,ik) = dipole_aux(:,:,:) * tpiba
     !
   ENDDO  ! k-points
   !
   !
   WRITE(stdout,'(/5x,a)') 'Dipole matrix elements calculated'
   WRITE(stdout,*)
   !
   RETURN
END SUBROUTINE compute_pmn_para
!-----------------------------------------------------------------------
!!
!-----------------------------------------------------------------------
SUBROUTINE write_filukk
!-----------------------------------------------------------------------
!
!  Here we compute and write out the final ukk matrix which is used by
!  epw.x to localize the electron wavefuctions (and therefore the ep-matrix 
!  elements)
!  10/2008 Jesse Noffsinger UC Berkeley
!  07/2010 Fixed the rotation for ndimwin when lower bands are not included
!
   USE kinds,        ONLY : DP
   USE io_epw,       ONLY : iuukk
   USE wvfct,        ONLY : nbnd
   USE wannier,      ONLY : n_wannier, iknum, u_mat, u_mat_opt, lwindow
   USE epwcom,       ONLY : filukk
   USE constants_epw,ONLY : czero
#ifdef __PARA
   USE io_global,    ONLY : ionode
#endif
   !
   implicit none
   !
   complex(kind=DP), allocatable :: u_kc(:,:,:)
   integer :: jbnd, k_wan, ik, ndimwin(iknum)
   !
   !
   !
#ifdef __PARA
   IF (ionode) THEN
#endif
      !
      ndimwin(:) = 0
      DO ik = 1, iknum
         DO jbnd = 1, nbnd
            IF (lwindow(jbnd,ik)) ndimwin(ik) = ndimwin(ik) + 1
         ENDDO
      ENDDO
      ALLOCATE( u_kc(nbnd, n_wannier, iknum) )
      u_kc = czero
      !
      ! get the final rotation matrix, which is the product of the optimal
      ! subspace and the rotation among the n_wannier wavefunctions
      DO ik = 1, iknum
         !
         u_kc(1:ndimwin(ik),1:n_wannier,ik) = &
              matmul (u_mat_opt(1:ndimwin(ik),:,ik), u_mat(:,1:n_wannier,ik))
         !
      ENDDO
      !
      open (unit = iuukk, file = filukk, form = 'formatted')
      DO ik = 1,iknum
         DO jbnd = 1, nbnd
            DO k_wan = 1, n_wannier
               WRITE (iuukk,*) u_kc(jbnd,k_wan,ik)
            ENDDO
         ENDDO
      ENDDO
      ! needs also lwindow when disentanglement is used
      DO ik = 1,iknum
         DO jbnd = 1,nbnd
            WRITE (iuukk,*) lwindow(jbnd,ik)
         ENDDO
      ENDDO
      close (iuukk)
      IF ( ALLOCATED(u_kc) ) DEALLOCATE(u_kc)
#ifdef __PARA
   ENDIF
#endif
   !
!-----------------------------------------------------------------------
END SUBROUTINE write_filukk
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
SUBROUTINE phases_a_m
!-----------------------------------------------------------------------
! will set phases here on the matrices.  Should not affect the spreads and
! centers found in w90, but will leave the u_mat_opt and u_mat to reflect the
! known phases
!
!
   USE mp_global,       ONLY : inter_pool_comm
   USE mp,              ONLY : mp_sum
   USE kinds,           ONLY : DP
   USE io_global,       ONLY : ionode
   USE klist,           ONLY : nkstot, xk, nks
   USE wvfct,           ONLY : nbnd
   USE wannier,         ONLY : a_mat, m_mat, n_wannier, nnb, kpb, iknum
   USE elph2,           ONLY : umat, umat_all
   USE constants_epw,   ONLY : czero, cone

   implicit none

   integer                      ::  ik, ipool, ib, ikb, i,nkq, &
                                    ik_g
   real(kind=DP)                ::  xktot(3,nkstot)
   complex(kind=DP), allocatable ::  a_mat_tmp(:,:,:), m_mn_tmp1(:,:), &
       m_mn_tmp2(:,:), m_mn_tmp3(:,:,:,:)
   real(DP), dimension(3)   :: zero_vect

   xktot = 0.d0
   IF (ionode) then
      DO i = 1,nkstot
         xktot(:,i) = xk(:,i)
      ENDDO
   ENDIF
   CALL mp_sum(xktot, inter_pool_comm)
   !
   ALLOCATE(m_mn_tmp3(nbnd,nbnd,nnb,iknum))
   ALLOCATE(a_mat_tmp(nbnd,n_wannier,iknum)) 
   ALLOCATE(m_mn_tmp1(nbnd,nbnd))
   ALLOCATE(m_mn_tmp2(nbnd,nbnd))
   !
   ! zero all temporary/work quantities
   !
   zero_vect = 0.0d0
   a_mat_tmp = czero
   m_mn_tmp1 = czero
   m_mn_tmp2 = czero
   m_mn_tmp3 = czero
   !
   DO ik=1,nks
      CALL ktokpmq ( xk(:,ik), zero_vect, +1, ipool, nkq, ik_g )
      !
      !  GF_n are the guiding functions which are our initial guesses 
      !  Amn(k) = <psi_k,m|GF_n>.  
      !  We want U(k)^dagger<psi_k,m|GF_m>
      !
      CALL zgemm ('c', 'n', nbnd, n_wannier, nbnd, cone, umat(:,:,ik), & 
           nbnd, a_mat(:,:,ik_g), nbnd, czero, a_mat_tmp(:,:,ik_g), nbnd)
      !
   ENDDO
   CALL mp_sum(a_mat_tmp, inter_pool_comm)
   !
   a_mat(:,:,:) = a_mat_tmp(:,:,:)
   !
   !
   DO ik=1,nks
      CALL ktokpmq ( xk(:,ik), zero_vect, +1, ipool, nkq, ik_g )
      !
      DO ib = 1,nnb
         ikb = kpb(ik_g,ib)
         !
         ! Mmn(k,k+b)  = <psi_k_m| psi_(k+b)_n> so we need
         !  (U(k)^dagger <psi_k_m| ) * (|psi_k+b_n> U(k+b)
         ! = U(k)^dagger (M_mn) = m_mat_tmp, 
         ! Mmn(k,k+b)' = m_mat_tmp*U(k+b) 
         !
         CALL zgemm ('c', 'n', nbnd, nbnd, nbnd, cone, umat(:,:,ik), & 
              nbnd, m_mat(:,:,ib,ik_g), nbnd, czero, m_mn_tmp1(:,:), nbnd)
         CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, m_mn_tmp1(:,:), & 
              nbnd, umat_all(:,:,ikb), nbnd, czero, m_mn_tmp2(:,:), nbnd)
         !         !
         !         m_mn_tmp1 = matmul ( conjg( transpose (umat(:,:,ik) )), m_mat(:,:,ib,ik_g ) )
         !         m_mn_tmp2 = matmul ( m_mn_tmp1, umat_g(:,:,ikb) )
         !
         m_mn_tmp3(:,:,ib,ik_g) = m_mn_tmp2
      ENDDO
   ENDDO
   CALL mp_sum(m_mn_tmp3, inter_pool_comm)

   m_mat(:,:,:,:) = m_mn_tmp3(:,:,:,:)
   !
   DEALLOCATE(m_mn_tmp3, m_mn_tmp2, m_mn_tmp1, a_mat_tmp)
   !
   !
   RETURN
!-----------------------------------------------------------------------
END SUBROUTINE phases_a_m
!-----------------------------------------------------------------------
!
!---------------------------------------
!
! SUBROUTINEs below are largely unchanged
! from pw2wannier90.x of esp-4.0.1
!
!---------------------------------------
!
SUBROUTINE generate_guiding_functions(ik)
   !
   USE kinds,          ONLY : DP
   USE mp_global,      ONLY : intra_pool_comm
   USE mp,             ONLY : mp_sum
   USE constants,      ONLY : tpi, eps8
   USE wvfct,          ONLY : npw, igk
   USE gvect,          ONLY : g
   USE cell_base,      ONLY : tpiba
   USE wannier,        ONLY : n_proj, gf, center_w, csph, alpha_w, &
                              r_w
   USE klist,          ONLY : xk 

   implicit none

   integer, parameter :: lmax=3, lmax2=(lmax+1)**2
   integer :: iw, ig, ik, l
   integer :: lm, iig
   real(DP) :: arg, anorm
   complex(DP) :: ZDOTC, lphase
   real(DP), allocatable :: gk(:,:), qg(:), ylm(:,:), radial(:,:)
   complex(DP), allocatable :: sk(:) 
   !
   ALLOCATE( gk(3,npw), qg(npw), ylm(npw,lmax2), sk(npw), radial(npw,0:lmax) )
   !
   DO ig = 1, npw
      gk (1,ig) = xk(1, ik) + g(1, igk(ig) )
      gk (2,ig) = xk(2, ik) + g(2, igk(ig) )
      gk (3,ig) = xk(3, ik) + g(3, igk(ig) )
      qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
   ENDDO

   CALL ylmr2 (lmax2, npw, gk, qg, ylm)
   ! define qg as the norm of (k+g) in a.u.
   qg(:) = sqrt(qg(:)) * tpiba
   ! RM changed according to QE4.0.3/PP/pw2wannier90 
   Do iw = 1, n_proj
   !DO iw = 1, n_wannier
      !
      gf(:,iw) = (0.d0,0.d0)

      CALL radialpart(npw, qg, alpha_w(iw), r_w(iw), lmax, radial) 

      DO lm = 1, lmax2
         IF ( abs(csph(lm,iw)) < eps8 ) CYCLE
         l = int (sqrt( lm-1.d0))
         lphase = (0.d0,-1.d0)**l
         !
         DO ig=1,npw
            gf(ig,iw) = gf(ig,iw) + csph(lm,iw) * ylm(ig,lm) * radial(ig,l) * lphase
         ENDDO !ig
      ENDDO ! lm
      DO ig=1,npw
         iig = igk(ig)
         arg = ( gk(1,ig)*center_w(1,iw) + gk(2,ig)*center_w(2,iw) + &
                                           gk(3,ig)*center_w(3,iw) ) * tpi
         ! center_w are cartesian coordinates in units of alat 
         sk(ig) = CMPLX(cos(arg), -sin(arg) )
         gf(ig,iw) = gf(ig,iw) * sk(ig) 
      ENDDO
      anorm = REAL(ZDOTC(npw,gf(1,iw),1,gf(1,iw),1))
      CALL mp_sum(anorm, intra_pool_comm)
      gf(:,iw) = gf(:,iw) / dsqrt(anorm)
   ENDDO
   !
   DEALLOCATE ( gk, qg, ylm, sk, radial)
   !
   RETURN
END SUBROUTINE generate_guiding_functions

SUBROUTINE write_band
   USE wvfct,         ONLY : nbnd, et
   USE constants_epw, ONLY: ryd2ev
   USE wannier

   implicit none

   integer ik, ibnd, ibnd1, ikevc

   ALLOCATE(eigval(num_bands,iknum))

   DO ik=ikstart,ikstop
      ikevc = ik - ikstart + 1
      ibnd1=0
      DO ibnd=1,nbnd
         IF (excluded_band(ibnd)) CYCLE
         ibnd1=ibnd1 + 1
         eigval(ibnd1,ikevc) = et(ibnd,ik)*ryd2ev
      ENDDO
   ENDDO
   !
   RETURN
END SUBROUTINE write_band

!---------------------------------------
SUBROUTINE write_plot
!---------------------------------------
!
! JN 06/2009: 
! added a couple of calls -- now works with multiple
! pools/procs (but one proc per pool)
!
   USE kinds,           ONLY : DP
   USE io_global,       ONLY : stdout
   USE io_epw,          ONLY : iun_plot
   USE wvfct,           ONLY : nbnd, npw, igk, g2kin
   USE gvecw,           ONLY : ecutwfc
   USE control_flags,   ONLY : gamma_only
   USE wavefunctions_module, ONLY : evc, psic
   USE units_ph,        ONLY : lrwfc, iuwfc
   USE wannier,         ONLY : reduce_unk, wvfn_formatted, ispinw, nexband, &
                               excluded_band 
   USE gvecs,           ONLY : nls, nlsm
   USE klist,           ONLY : xk, nks
   USE gvect,           ONLY : g, ngm 
   USE cell_base,       ONLY : tpiba2
   USE fft_base,        ONLY : dffts
   USE fft_interfaces,  ONLY : invfft
   USE noncollin_module,ONLY : noncolin
#ifdef __PARA
   USE scatter_mod,        ONLY : gather_grid
#endif
   !
   implicit none
   integer ik, ibnd, ibnd1, j, spin, nkq, nkq_abs, ipool
   character(len=20) wfnname

   ! aam: 1/5/06: for writing smaller unk files 
   integer :: n1by2,n2by2,n3by2,i,k,idx,pos
   COMPLEX(DP),allocatable :: psic_small(:)   
   real(kind=DP)      :: zero_vect(3)
   !-------------------------------------------!
#ifdef __PARA
   integer nxxs
   COMPLEX(DP),allocatable :: psic_all(:)
   nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x
   ALLOCATE(psic_all(nxxs) )
#endif
   !
   IF (noncolin) CALL errore('pw2wan90epw',&
       'write_unk not implemented with noncolin',1)
   !
   zero_vect = 0.d0
   WRITE (stdout,*) '    Writing out UNK plot files'
   !
   IF (reduce_unk) then
      WRITE(stdout,'(3(a,i5))') 'nr1s =',dffts%nr1,'nr2s=',dffts%nr2,'nr3s=',dffts%nr3
      n1by2=(dffts%nr1+1)/2;n2by2=(dffts%nr2+1)/2;n3by2=(dffts%nr3+1)/2
      WRITE(stdout,'(3(a,i5))') 'n1by2=',n1by2,'n2by2=',n2by2,'n3by2=',n3by2
      ALLOCATE(psic_small(n1by2*n2by2*n3by2))   
   ENDIF
   !
   DO ik=1,nks
      CALL ktokpmq ( xk(:,ik), zero_vect, +1, ipool, nkq, nkq_abs )
      !
      spin=ispinw
      IF(ispinw.eq.0) spin=1
      WRITE(wfnname,200) nkq_abs, spin
200   FORMAT ('UNK',i5.5,'.',i1)
      !
      IF (wvfn_formatted) THEN 
         OPEN (unit=iun_plot, file=wfnname,form='formatted')
         IF (reduce_unk) THEN
            WRITE(iun_plot,*)  n1by2,n2by2,n3by2, nkq_abs, nbnd-nexband
         ELSE
            WRITE(iun_plot,*)  dffts%nr1,dffts%nr2,dffts%nr3, nkq_abs, nbnd-nexband
         ENDIF
      ELSE
         OPEN (unit=iun_plot, file=wfnname,form='unformatted')
         IF (reduce_unk) THEN
            WRITE(iun_plot)  n1by2,n2by2,n3by2, nkq_abs, nbnd-nexband
         ELSE
            WRITE(iun_plot)  dffts%nr1,dffts%nr2,dffts%nr3, nkq_abs, nbnd-nexband
         ENDIF
      ENDIF
      !
      CALL davcio (evc, lrwfc, iuwfc, ik, -1 )
      CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      !
      ibnd1 = 0
      DO ibnd=1,nbnd
         IF (excluded_band(ibnd)) CYCLE
         ibnd1=ibnd1 + 1
         psic(:) = (0.d0, 0.d0)
         psic(nls (igk (1:npw) ) ) = evc (1:npw, ibnd)
         IF (gamma_only)  psic(nlsm(igk (1:npw) ) ) = conjg(evc (1:npw, ibnd))
         CALL invfft ('Wave', psic, dffts)
         IF (reduce_unk) pos=0
#ifdef __PARA
         CALL gather_grid(dffts,psic,psic_all)
         !
         IF (reduce_unk) THEN
            DO k=1,dffts%nr3,2
               DO j=1,dffts%nr2,2
                  DO i=1,dffts%nr1,2
                     idx = (k-1)*dffts%nr3*dffts%nr2 + (j-1)*dffts%nr2 + i
                     pos=pos+1
                     psic_small(pos) = psic_all(idx) 
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         !
         IF (wvfn_formatted) THEN
            IF (reduce_unk) THEN
               WRITE (iun_plot,'(2ES20.10)') (psic_small(j),j=1,n1by2*n2by2*n3by2)
            ELSE
               WRITE (iun_plot,*) (psic_all(j),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
            ENDIF
         ELSE
            IF (reduce_unk) THEN
               WRITE (iun_plot) (psic_small(j),j=1,n1by2*n2by2*n3by2)
            ELSE
               WRITE (iun_plot) (psic_all(j),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
            ENDIF
         ENDIF
#else
         IF (reduce_unk) THEN
            DO k=1,dffts%nr3,2
               DO j=1,dffts%nr2,2
                  DO i=1,dffts%nr1,2
                     idx = (k-1)*dffts%nr3*dffts%nr2 + (j-1)*dffts%nr2 + i
                     pos=pos+1
                     psic_small(pos) = psic(idx) 
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         IF (wvfn_formatted) THEN
            IF (reduce_unk) THEN
               WRITE (iun_plot,'(2ES20.10)') (psic_small(j),j=1,n1by2*n2by2*n3by2)
            ELSE
               WRITE (iun_plot,*) (psic(j),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
            ENDIF
         ELSE
            IF (reduce_unk) THEN
               WRITE (iun_plot) (psic_small(j),j=1,n1by2*n2by2*n3by2)
            ELSE
               WRITE (iun_plot) (psic(j),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
            ENDIF
         ENDIF
#endif
      ENDDO !ibnd
      !
      CLOSE (iun_plot)
      !
   ENDDO  !ik
   !
   IF (reduce_unk) DEALLOCATE(psic_small)   
   !
#ifdef __PARA
   DEALLOCATE( psic_all )
#endif
   !
   RETURN
END SUBROUTINE write_plot
!---------------------------------------
SUBROUTINE ylm_expansion 
   USE kinds,         ONLY : DP
   USE random_numbers,ONLY : randy
   USE wannier,       ONLY : n_proj , xaxis, zaxis, csph, l_w, mr_w
   USE matrix_inversion, ONLY: invmat
   implicit none
   ! local variables
   integer, parameter :: lmax2=16
   integer ::  i, ir, iw
   real(DP) :: capel
   real(DP), allocatable :: r(:,:), rr(:), rp(:,:), ylm_w(:), ylm(:,:), mly(:,:)
   real(DP) :: u(3,3)

   ALLOCATE (r(3,lmax2), rp(3,lmax2), rr(lmax2), ylm_w(lmax2))
   ALLOCATE (ylm(lmax2,lmax2), mly(lmax2,lmax2) )

   ! generate a set of nr=lmax2 random vectors
   DO ir=1,lmax2
      DO i=1,3
         r(i,ir) = randy() -0.5d0
      ENDDO
   ENDDO
   rr(:) = r(1,:)*r(1,:) + r(2,:)*r(2,:) + r(3,:)*r(3,:)
   !- compute ylm(ir,lm)
   CALL ylmr2(lmax2, lmax2, r, rr, ylm)
   !- store the inverse of ylm(ir,lm) in mly(lm,ir)
   CALL invmat(lmax2, ylm, mly, capel)
   !- check that r points are independent
   CALL check_inverse(lmax2, ylm, mly)
   !
   ! RM changed according to QE4.0.3/PP/pw2wannier90 
   DO iw=1, n_proj 
   !DO iw=1, n_wannier

      !- define the u matrix that rotate the reference frame
      CALL set_u_matrix (xaxis(:,iw),zaxis(:,iw),u)
      !- find rotated r-vectors 
      rp(:,:) = matmul ( u(:,:) , r(:,:) )
      !- set ylm funtion according to wannier90 (l,mr) indexing in the rotaterd points
      CALL ylm_wannier(ylm_w,l_w(iw),mr_w(iw),rp,lmax2) 

      csph(:,iw) = matmul (mly(:,:), ylm_w(:))

!      write (stdout,*) 
!      write (stdout,'(2i4,2(2x,3f6.3))') l_w(iw), mr_w(iw), xaxis(:,iw), zaxis(:,iw)
!      write (stdout,'(16i6)')   (lm, lm=1,lmax2)
!      write (stdout,'(16f6.3)') (csph(lm,iw), lm=1,lmax2)

   ENDDO
   DEALLOCATE (r, rp, rr, ylm_w, ylm, mly )
   !
   RETURN
END SUBROUTINE ylm_expansion
!-----------------------------------------------
SUBROUTINE check_inverse(lmax2, ylm, mly)
   USE kinds,     ONLY : DP
   USE constants, ONLY : eps8
   implicit none
   ! I/O variables
   integer :: lmax2
   real(DP) :: ylm(lmax2,lmax2), mly(lmax2,lmax2)
   ! local variables
   real(DP), allocatable :: uno(:,:)
   real(DP) :: capel
   integer :: lm
   !
   ALLOCATE (uno(lmax2,lmax2) )
   uno = matmul(mly, ylm)
   capel = 0.d0
   DO lm = 1, lmax2
      uno(lm,lm) = uno(lm,lm) - 1.d0
   ENDDO
   capel = capel + SUM ( abs(uno(1:lmax2,1:lmax2) ) )
!   write (stdout,*) "capel = ", capel
   IF (capel > eps8) CALL errore('ylm_expansion', &
                    ' inversion failed: r(*,1:nr) are not all independent !!',1)
   DEALLOCATE (uno)
   !
   RETURN
END SUBROUTINE check_inverse
!--------------------------------------------   
SUBROUTINE set_u_matrix(x,z,u)
   USE kinds,     ONLY : DP
   USE constants, ONLY : eps8
   implicit none
   ! I/O variables
   real(DP) :: x(3),z(3),u(3,3)
   ! local variables
   real(DP) :: xx, zz, y(3), coseno

   xx = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
   IF (xx < eps8) CALL errore ('set_u_matrix',' |xaxis| < eps ',1)
!   x(:) = x(:)/xx
   zz = sqrt(z(1)*z(1) + z(2)*z(2) + z(3)*z(3))
   IF (zz < eps8) CALL errore ('set_u_matrix',' |zaxis| < eps ',1)
!   z(:) = z(:)/zz

   coseno = (x(1)*z(1) + x(2)*z(2) + x(3)*z(3))/xx/zz
   IF (abs(coseno) > eps8) CALL errore('set_u_matrix',' xaxis and zaxis are not orthogonal !',1)

   y(1) = (z(2)*x(3) - x(2)*z(3))/xx/zz
   y(2) = (z(3)*x(1) - x(3)*z(1))/xx/zz
   y(3) = (z(1)*x(2) - x(1)*z(2))/xx/zz

   u(1,:) = x(:)/xx
   u(2,:) = y(:)
   u(3,:) = z(:)/zz

!   write (stdout,'(3f10.7)') u(:,:)


END SUBROUTINE set_u_matrix

SUBROUTINE ylm_wannier(ylm,l,mr,r,nr) 
!
! this routine returns in ylm(r) the values at the nr points r(1:3,1:nr) 
! of the spherical harmonic identified  by indices (l,mr) 
! in table 3.1 of the wannierf90 specification.
! 
! No reference to the particular ylm ordering internal to quantum-espresso
! is assumed. 
!
! If ordering in wannier90 code is changed or extended this should be the 
! ONLY place to be modified accordingly
!
   USE kinds,     ONLY : DP
   USE constants, ONLY : pi, eps8
   implicit none
! I/O variables
!
   integer :: l, mr, nr
   real(DP) :: ylm(nr), r(3,nr)
!
! local variables
!
   real(DP), external :: s, p_z,px,py, dz2, dxz, dyz, dx2my2, dxy, &
                        fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
   real(DP) :: rr, cost, phi
   integer :: ir
   real(DP) :: bs2, bs3, bs6, bs12
   bs2 = 1.d0/sqrt(2.d0)
   bs3=1.d0/sqrt(3.d0)
   bs6 = 1.d0/sqrt(6.d0)
   bs12 = 1.d0/sqrt(12.d0)
!
   IF (l > 3 .OR. l < -5 ) CALL errore('ylm_wannier',' l out of range ', 1)
   IF (l>=0) then
      IF (mr < 1 .OR. mr > 2*l+1) CALL errore('ylm_wannier','mr out of range' ,1)
   ELSE
      IF (mr < 1 .OR. mr > abs(l)+1 ) CALL errore('ylm_wannier','mr out of range',1)
   end if

   DO ir=1, nr
      rr = sqrt( r(1,ir)*r(1,ir) +  r(2,ir)*r(2,ir) + r(3,ir)*r(3,ir) )
      IF (rr < eps8) CALL errore('ylm_wannier',' rr too small ',1)

      cost =  r(3,ir) / rr
      !
      !  beware the arc tan, it is defined modulo pi
      !
      IF (r(1,ir) > eps8) THEN
         phi = atan( r(2,ir)/r(1,ir) )
      ELSEIF (r(1,ir) < -eps8 ) THEN
         phi = atan( r(2,ir)/r(1,ir) ) + pi
      ELSE
         phi = sign( pi/2.d0,r(2,ir) )
      ENDIF

    
      IF (l==0) THEN   ! s orbital
                    ylm(ir) = s()  
      ENDIF
      IF (l==1) THEN   ! p orbitals
         IF (mr==1) ylm(ir) = p_z(cost) 
         IF (mr==2) ylm(ir) = px(cost,phi)
         IF (mr==3) ylm(ir) = py(cost,phi)
      ENDIF
      IF (l==2) THEN   ! d orbitals
         IF (mr==1) ylm(ir) = dz2(cost)
         IF (mr==2) ylm(ir) = dxz(cost,phi)
         IF (mr==3) ylm(ir) = dyz(cost,phi)
         IF (mr==4) ylm(ir) = dx2my2(cost,phi)
         IF (mr==5) ylm(ir) = dxy(cost,phi)
      ENDIF
      IF (l==3) THEN   ! f orbitals
         IF (mr==1) ylm(ir) = fz3(cost)
         IF (mr==2) ylm(ir) = fxz2(cost,phi)
         IF (mr==3) ylm(ir) = fyz2(cost,phi)
         IF (mr==4) ylm(ir) = fzx2my2(cost,phi)
         IF (mr==5) ylm(ir) = fxyz(cost,phi)
         IF (mr==6) ylm(ir) = fxx2m3y2(cost,phi)
         IF (mr==7) ylm(ir) = fy3x2my2(cost,phi)
      ENDIF
      IF (l==-1) THEN  !  sp hybrids
         IF (mr==1) ylm(ir) = bs2 * ( s() + px(cost,phi) ) 
         IF (mr==2) ylm(ir) = bs2 * ( s() - px(cost,phi) ) 
      ENDIF
      IF (l==-2) THEN  !  sp2 hybrids 
         IF (mr==1) ylm(ir) = bs3*s()-bs6*px(cost,phi)+bs2*py(cost,phi)
         IF (mr==2) ylm(ir) = bs3*s()-bs6*px(cost,phi)-bs2*py(cost,phi)
         IF (mr==3) ylm(ir) = bs3*s() +2.d0*bs6*px(cost,phi) 
      ENDIF
      IF (l==-3) THEN  !  sp3 hybrids
         IF (mr==1) ylm(ir) = 0.5d0*(s()+px(cost,phi)+py(cost,phi)+p_z(cost))
         IF (mr==2) ylm(ir) = 0.5d0*(s()+px(cost,phi)-py(cost,phi)-p_z(cost))
         IF (mr==3) ylm(ir) = 0.5d0*(s()-px(cost,phi)+py(cost,phi)-p_z(cost))
         IF (mr==4) ylm(ir) = 0.5d0*(s()-px(cost,phi)-py(cost,phi)+p_z(cost))
      ENDIF
      IF (l==-4) THEN  !  sp3d hybrids
         IF (mr==1) ylm(ir) = bs3*s()-bs6*px(cost,phi)+bs2*py(cost,phi)
         IF (mr==2) ylm(ir) = bs3*s()-bs6*px(cost,phi)-bs2*py(cost,phi)
         IF (mr==3) ylm(ir) = bs3*s() +2.d0*bs6*px(cost,phi) 
         IF (mr==4) ylm(ir) = bs2*p_z(cost)+bs2*dz2(cost)
         IF (mr==5) ylm(ir) =-bs2*p_z(cost)+bs2*dz2(cost)
      ENDIF
      IF (l==-5) THEN  ! sp3d2 hybrids
         IF (mr==1) ylm(ir) = bs6*s()-bs2*px(cost,phi)-bs12*dz2(cost)+.5*dx2my2(cost,phi)
         IF (mr==2) ylm(ir) = bs6*s()+bs2*px(cost,phi)-bs12*dz2(cost)+.5*dx2my2(cost,phi)
         IF (mr==3) ylm(ir) = bs6*s()-bs2*py(cost,phi)-bs12*dz2(cost)-.5*dx2my2(cost,phi)
         IF (mr==4) ylm(ir) = bs6*s()+bs2*py(cost,phi)-bs12*dz2(cost)-.5*dx2my2(cost,phi)
         IF (mr==5) ylm(ir) = bs6*s()-bs2*p_z(cost)+bs3*dz2(cost)
         IF (mr==6) ylm(ir) = bs6*s()+bs2*p_z(cost)+bs3*dz2(cost)
      ENDIF

   ENDDO
   !
   RETURN
END SUBROUTINE ylm_wannier

!======== l = 0 =====================================================================
FUNCTION s()
   USE kinds,         ONLY : DP
   USE constants_epw, ONLY : fpi
   implicit none
   real(DP) :: s
   s = 1.d0/ sqrt(fpi)
   RETURN
END FUNCTION s
!======== l = 1 =====================================================================
FUNCTION p_z(cost)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : fpi
   implicit none
   real(DP) ::p_z, cost
   p_z =  sqrt(3.d0/fpi) * cost
   RETURN
END FUNCTION p_z
FUNCTION px(cost,phi)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : fpi
   implicit none
   real(DP) ::px, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   px =  sqrt(3.d0/fpi) * sint * cos(phi)
   RETURN
END FUNCTION px
FUNCTION py(cost,phi)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : fpi
   implicit none
   real(DP) ::py, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   py =  sqrt(3.d0/fpi) * sint * sin(phi)
   RETURN
END FUNCTION py
!======== l = 2 =====================================================================
FUNCTION dz2(cost)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : fpi
   implicit none
   real(DP) ::dz2, cost
   dz2 =  sqrt(1.25d0/fpi) * (3.d0* cost*cost-1.d0)
   RETURN
END FUNCTION dz2
FUNCTION dxz(cost,phi)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : fpi
   implicit none
   real(DP) ::dxz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxz =  sqrt(15.d0/fpi) * sint*cost * cos(phi)
   RETURN
END FUNCTION dxz
FUNCTION dyz(cost,phi)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : fpi
   implicit none
   real(DP) ::dyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dyz =  sqrt(15.d0/fpi) * sint*cost * sin(phi)
   RETURN
END FUNCTION dyz
FUNCTION dx2my2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : fpi
   implicit none
   real(DP) ::dx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dx2my2 =  sqrt(3.75d0/fpi) * sint*sint * cos(2.d0*phi)
   RETURN
END FUNCTION dx2my2
FUNCTION dxy(cost,phi)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : fpi
   implicit none
   real(DP) ::dxy, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxy =  sqrt(3.75d0/fpi) * sint*sint * sin(2.d0*phi)
   RETURN
END FUNCTION dxy
!======== l = 3 =====================================================================
FUNCTION fz3(cost)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : pi
   implicit none
   real(DP) ::fz3, cost
   fz3 =  0.25d0*sqrt(7.d0/pi) * ( 5.d0 * cost * cost - 3.d0 ) * cost
   RETURN
END FUNCTION fz3
FUNCTION fxz2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : pi
   implicit none
   real(DP) ::fxz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * cos(phi)
   RETURN
END FUNCTION fxz2
FUNCTION fyz2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : pi
   implicit none
   real(DP) ::fyz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fyz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * sin(phi)
   RETURN
END FUNCTION fyz2
FUNCTION fzx2my2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : pi
   implicit none
   real(DP) ::fzx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fzx2my2 =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * cos(2.d0*phi)
   RETURN
END FUNCTION fzx2my2
FUNCTION fxyz(cost,phi)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : pi
   implicit none
   real(DP) ::fxyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxyz =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * sin(2.d0*phi)
   RETURN
END FUNCTION fxyz
FUNCTION fxx2m3y2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : pi
   implicit none
   real(DP) ::fxx2m3y2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * cos(3.d0*phi)
   RETURN
END FUNCTION fxx2m3y2
FUNCTION fy3x2my2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants_epw, ONLY : pi
   implicit none
   real(DP) ::fy3x2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * sin(3.d0*phi)
   RETURN
END FUNCTION fy3x2my2
!
!
!-----------------------------------------------------------------------
SUBROUTINE radialpart(ng, q, alfa, rvalue, lmax, radial)
  !-----------------------------------------------------------------------
  !
  ! This routine computes a table with the radial Fourier transform 
  ! of the radial functions.
  !
  USE kinds,      ONLY : dp
  USE constants,  ONLY : fpi
  USE cell_base,  ONLY : omega
  !
  implicit none
  ! I/O
  integer :: ng, rvalue, lmax
  real(DP) :: q(ng), alfa, radial(ng,0:lmax)
  ! local variables
  real(DP), parameter :: xmin=-6.d0, dx=0.025d0, rmax=10.d0

  real(DP) :: rad_int, pref, x
  integer :: l, ir, ig, mesh_r
  real(DP), allocatable :: bes(:), func_r(:), r(:), rij(:), aux(:)

  mesh_r = nint ( ( log ( rmax ) - xmin ) / dx + 1 )
  ALLOCATE ( bes(mesh_r), func_r(mesh_r), r(mesh_r), rij(mesh_r) )
  ALLOCATE ( aux(mesh_r))
  !
  !    compute the radial mesh
  !
  DO ir = 1, mesh_r
     x = xmin  + DBLE (ir - 1) * dx 
     r (ir) = exp (x) / alfa
     rij (ir) = dx  * r (ir)
  ENDDO
  !
  IF (rvalue==1) func_r(:) = 2.d0 * alfa**(3.d0/2.d0) * exp(-alfa*r(:))
  IF (rvalue==2) func_r(:) = 1.d0/sqrt(8.d0) * alfa**(3.d0/2.d0) * & 
                     (2.0d0 - alfa*r(:)) * exp(-alfa*r(:)*0.5d0)
  IF (rvalue==3) func_r(:) = sqrt(4.d0/27.d0) * alfa**(2.0d0/3.0d0) * &
                     (1.d0 - 1.5d0*alfa*r(:) + 2.d0*(alfa*r(:))**2/27.d0) * &
                                           exp(-alfa*r(:)/3.0d0)
  pref = fpi/sqrt(omega)
  !
  DO l = 0, lmax
     DO ig=1,ng
       CALL sph_bes (mesh_r, r(1), q(ig), l, bes)
       aux(:) = bes(:) * func_r(:) * r(:)
       CALL simpson (mesh_r, aux, rij, rad_int)
       radial(ig,l) = rad_int * pref
     ENDDO
  ENDDO

  DEALLOCATE (bes, func_r, r, rij, aux )
  
  RETURN
END SUBROUTINE radialpart

