!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE kpoint_grid ( nrot, time_reversal, skip_equivalence, s, t_rev, &
                         bg, npk, k1,k2,k3, nk1,nk2,nk3, nks, xk, wk)
!-----------------------------------------------------------------------
!
!  Automatic generation of a uniform grid of k-points
!
  USE kinds, ONLY: DP
  IMPLICIT NONE
  !
  INTEGER, INTENT(in):: nrot, npk, k1, k2, k3, nk1, nk2, nk3, &
                        t_rev(48), s(3,3,48)
  LOGICAL, INTENT(in):: time_reversal, skip_equivalence
  real(DP), INTENT(in):: bg(3,3)
  !
  INTEGER, INTENT(out) :: nks
  real(DP), INTENT(out):: xk(3,npk)
  real(DP), INTENT(out):: wk(npk)
  ! LOCAL:
  real(DP), PARAMETER :: eps=1.0d-5
  real(DP) :: xkr(3), fact, xx, yy, zz
  real(DP), ALLOCATABLE:: xkg(:,:), wkk(:)
  INTEGER :: nkr, i,j,k, ns, n, nk
  INTEGER, ALLOCATABLE :: equiv(:)
  LOGICAL :: in_the_list
  !
  nkr=nk1*nk2*nk3
  ALLOCATE (xkg( 3,nkr),wkk(nkr))
  ALLOCATE (equiv( nkr))
  !
  DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
           xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
           xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
        ENDDO
     ENDDO
  ENDDO

  !  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
  !  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)

  DO nk=1,nkr
     equiv(nk)=nk
  ENDDO

  IF ( skip_equivalence ) THEN
    CALL infomsg('kpoint_grid', 'ATTENTION: skip check of k-points equivalence')
    wkk = 1.d0
  ELSE
    DO nk=1,nkr
    !  check if this k-point has already been found equivalent to another
      IF (equiv(nk) == nk) THEN
        wkk(nk)   = 1.0d0
        !  check if there are equivalent k-point to this in the list
        !  (excepted those previously found to be equivalent to another)
        !  check both k and -k
        DO ns=1,nrot
           DO i=1,3
              xkr(i) = s(i,1,ns) * xkg(1,nk) &
                     + s(i,2,ns) * xkg(2,nk) &
                     + s(i,3,ns) * xkg(3,nk)
              xkr(i) = xkr(i) - nint( xkr(i) )
           ENDDO
           IF(t_rev(ns)==1) xkr = -xkr
           xx = xkr(1)*nk1 - 0.5d0*k1
           yy = xkr(2)*nk2 - 0.5d0*k2
           zz = xkr(3)*nk3 - 0.5d0*k3
           in_the_list = abs(xx-nint(xx))<=eps .and. &
                         abs(yy-nint(yy))<=eps .and. &
                         abs(zz-nint(zz))<=eps
           IF (in_the_list) THEN
              i = mod ( nint ( xkr(1)*nk1 - 0.5d0*k1 + 2*nk1), nk1 ) + 1
              j = mod ( nint ( xkr(2)*nk2 - 0.5d0*k2 + 2*nk2), nk2 ) + 1
              k = mod ( nint ( xkr(3)*nk3 - 0.5d0*k3 + 2*nk3), nk3 ) + 1
              n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
              IF (n>nk .and. equiv(n)==n) THEN
                 equiv(n) = nk
                 wkk(nk)=wkk(nk)+1.0d0
              ELSE
                 IF (equiv(n)/=nk .or. n<nk ) CALL errore('kpoint_grid', &
                    'something wrong in the checking algorithm',1)
              ENDIF
           ENDIF
           IF ( time_reversal ) THEN
              xx =-xkr(1)*nk1 - 0.5d0*k1
              yy =-xkr(2)*nk2 - 0.5d0*k2
              zz =-xkr(3)*nk3 - 0.5d0*k3
              in_the_list=abs(xx-nint(xx))<=eps.and.abs(yy-nint(yy))<=eps &
                                                 .and. abs(zz-nint(zz))<=eps
              IF (in_the_list) THEN
                 i = mod ( nint (-xkr(1)*nk1 - 0.5d0 * k1 + 2*nk1), nk1 ) + 1
                 j = mod ( nint (-xkr(2)*nk2 - 0.5d0 * k2 + 2*nk2), nk2 ) + 1
                 k = mod ( nint (-xkr(3)*nk3 - 0.5d0 * k3 + 2*nk3), nk3 ) + 1
                 n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                 IF (n>nk .and. equiv(n)==n) THEN
                    equiv(n) = nk
                    wkk(nk)=wkk(nk)+1.0d0
                 ELSE
                    IF (equiv(n)/=nk.or.n<nk) CALL errore('kpoint_grid', &
                    'something wrong in the checking algorithm',2)
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  !  count irreducible points and order them

  nks=0
  fact=0.0d0
  DO nk=1,nkr
     IF (equiv(nk)==nk) THEN
        nks=nks+1
        IF (nks>npk) CALL errore('kpoint_grid','too many k-points',1)
        wk(nks) = wkk(nk)
        fact    = fact+wk(nks)
        !  bring back into to the first BZ
        DO i=1,3
           xk(i,nks) = xkg(i,nk)-nint(xkg(i,nk))
        ENDDO
     ENDIF
  ENDDO
  !  go to cartesian axis (in units 2pi/a0)
  CALL cryst_to_cart(nks,xk,bg,1)
  !  normalize weights to one
  DO nk=1,nks
     wk(nk) = wk(nk)/fact
  ENDDO

  DEALLOCATE(equiv)
  DEALLOCATE(xkg,wkk)

  RETURN
END SUBROUTINE kpoint_grid
!
!-----------------------------------------------------------------------
SUBROUTINE tetrahedra ( nsym, s, time_reversal, t_rev, at, bg, npk, &
     k1,k2,k3, nk1,nk2,nk3, nks, xk, wk, ntetra, tetra )
  !-----------------------------------------------------------------------
  !
  ! Tetrahedron method according to P. E. Bloechl et al, PRB49, 16223 (1994)
  !
  USE kinds, ONLY: DP
  IMPLICIT NONE
  ! 
  INTEGER, INTENT(IN):: nks, nsym, t_rev(48), s(3,3,48), npk, &
                        k1, k2, k3, nk1, nk2, nk3, ntetra
  LOGICAL, INTENT (IN) :: time_reversal
  real(DP), INTENT(IN) :: at(3,3), bg(3,3), xk(3,npk), wk(npk)
  !
  INTEGER, INTENT(OUT) :: tetra(4,ntetra)
  ! 
  real(DP) :: xkr(3), deltap(3), deltam(3)
  real(DP), PARAMETER:: eps=1.0d-5
  real(DP), ALLOCATABLE :: xkg(:,:)
  INTEGER :: nkr, i,j,k, ns, n, nk, ip1,jp1,kp1, &
       n1,n2,n3,n4,n5,n6,n7,n8
  INTEGER, ALLOCATABLE:: equiv(:)
  !
  ! Re-generate a uniform grid of k-points xkg
  !
  nkr=nk1*nk2*nk3
  !      ntetra=6*nkr
  ALLOCATE (xkg( 3,nkr))
  ALLOCATE (equiv( nkr))
!
  DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
           xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
           xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
        ENDDO
     ENDDO
  ENDDO

  !  locate k-points of the uniform grid in the list of irreducible k-points
  !  that was previously calculated

  !  bring irreducible k-points to crystal axis
  CALL cryst_to_cart (nks,xk,at,-1)
  !
  DO nk=1,nkr
     DO n=1,nks
        DO ns=1,nsym
           DO i=1,3
              xkr(i) = s(i,1,ns) * xk(1,n) + &
                       s(i,2,ns) * xk(2,n) + &
                       s(i,3,ns) * xk(3,n)
           ENDDO
           IF(t_rev(ns)==1) xkr = -xkr
           !  xkr is the n-th irreducible k-point rotated wrt the ns-th symmetry
           DO i=1,3
              deltap(i) = xkr(i)-xkg(i,nk) - nint (xkr(i)-xkg(i,nk) )
              deltam(i) = xkr(i)+xkg(i,nk) - nint (xkr(i)+xkg(i,nk) )
           ENDDO
           !  deltap is the difference vector, brought back in the first BZ
           !  deltam is the same but with k => -k (for time reversal)
           IF ( sqrt ( deltap(1)**2 + &
                       deltap(2)**2 + &
                       deltap(3)**2 ) < eps .or. ( time_reversal .and. &
                sqrt ( deltam(1)**2 +  &
                       deltam(2)**2 +  &
                       deltam(3)**2 ) < eps ) ) THEN
              !  equivalent irreducible k-point found
              equiv(nk) = n
              GOTO 15
           ENDIF
        ENDDO
     ENDDO
     !  equivalent irreducible k-point found - something wrong
     CALL errore('tetrahedra','cannot locate  k point',nk)
15   CONTINUE
  ENDDO

  DO n=1,nks
     DO nk=1,nkr
        IF (equiv(nk)==n) GOTO 20
     ENDDO
     !  this failure of the algorithm may indicate that the displaced grid
     !  (with k1,k2,k3.ne.0) does not have the full symmetry of the lattice
     CALL errore('tetrahedra','cannot remap grid on k-point list',n)
20   CONTINUE
  ENDDO

  !  bring irreducible k-points back to cartesian axis
  CALL cryst_to_cart (nks,xk,bg, 1)

  !  construct tetrahedra

  DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           !  n1-n8 are the indices of k-point 1-8 forming a cube
           ip1 = mod(i,nk1)+1
           jp1 = mod(j,nk2)+1
           kp1 = mod(k,nk3)+1
           n1 = (  k-1) + (  j-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n2 = (  k-1) + (  j-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n3 = (  k-1) + (jp1-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n4 = (  k-1) + (jp1-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n5 = (kp1-1) + (  j-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n6 = (kp1-1) + (  j-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n7 = (kp1-1) + (jp1-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n8 = (kp1-1) + (jp1-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           !  there are 6 tetrahedra per cube (and nk1*nk2*nk3 cubes)
           n  = 6 * ( (k-1) + (j-1)*nk3 + (i-1)*nk3*nk2 )

           tetra (1,n+1) = equiv(n1)
           tetra (2,n+1) = equiv(n2)
           tetra (3,n+1) = equiv(n3)
           tetra (4,n+1) = equiv(n6)

           tetra (1,n+2) = equiv(n2)
           tetra (2,n+2) = equiv(n3)
           tetra (3,n+2) = equiv(n4)
           tetra (4,n+2) = equiv(n6)

           tetra (1,n+3) = equiv(n1)
           tetra (2,n+3) = equiv(n3)
           tetra (3,n+3) = equiv(n5)
           tetra (4,n+3) = equiv(n6)

           tetra (1,n+4) = equiv(n3)
           tetra (2,n+4) = equiv(n4)
           tetra (3,n+4) = equiv(n6)
           tetra (4,n+4) = equiv(n8)

           tetra (1,n+5) = equiv(n3)
           tetra (2,n+5) = equiv(n6)
           tetra (3,n+5) = equiv(n7)
           tetra (4,n+5) = equiv(n8)

           tetra (1,n+6) = equiv(n3)
           tetra (2,n+6) = equiv(n5)
           tetra (3,n+6) = equiv(n6)
           tetra (4,n+6) = equiv(n7)
        ENDDO
     ENDDO
  ENDDO

  !  check

  DO n=1,ntetra
     DO i=1,4
        IF ( tetra(i,n)<1 .or. tetra(i,n)>nks ) &
             CALL errore ('tetrahedra','something wrong',n)
     ENDDO
  ENDDO

  DEALLOCATE(equiv)
  DEALLOCATE(xkg)

  RETURN
END SUBROUTINE tetrahedra

!-----------------------------------------------------------------------
SUBROUTINE kpoint_grid_efield (at, bg, npk, &
                         k1,k2,k3, nk1,nk2,nk3, nks, xk, wk, nspin)
!-----------------------------------------------------------------------
!
!  Automatic generation of a uniform grid of k-points, for Berry's phase electric field
!
  USE kinds, ONLY : DP
  USE bp,    ONLY : nppstr_3d, nx_el, l3dstring, efield_cart, efield_cry,&
                      transform_el
  USE io_global,  ONLY : stdout
  USE noncollin_module,   ONLY : noncolin

  IMPLICIT NONE
  !
  INTEGER, INTENT(in):: npk, k1, k2, k3, nk1, nk2, nk3,nspin
  real(DP), INTENT(in):: bg(3,3), at(3,3)
  !
  INTEGER, INTENT(out) :: nks
  real(DP), INTENT(out):: xk(3,npk)
  real(DP), INTENT(out):: wk(npk)

  INTEGER :: i,j,k,n,nk,m
  INTEGER :: nppstr_max
  real(DP) :: fact, sca
  real(DP) :: cry_to_cart(3,3)
  real(DP) :: bg_n(3,3)
  !

  !
  DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xk(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
           xk(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
           xk(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
        ENDDO
     ENDDO
  ENDDO

  nks=nk1*nk2*nk3
  !  go to cartesian axis (in units 2pi/a0)
  CALL cryst_to_cart(nks,xk,bg,1)
  fact=1.d0/dble(nks)
  !  normalize weights to one
  DO nk=1,nks
     wk(nk) = fact
  ENDDO

!setup nppstr_3d
  nppstr_3d(1)=nk1
  nppstr_3d(2)=nk2
  nppstr_3d(3)=nk3

!allocate and set up correspondence
  nppstr_max=nk1*nk2*nk3

  IF(noncolin) THEN
     ALLOCATE(nx_el(nppstr_max,3))
  ELSE
     ALLOCATE(nx_el(nppstr_max*nspin,3))
  END IF

 !establih correspondence

   DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           nx_el(n,3)=n
           m = (i-1) + (k-1)*nk1 + (j-1)*nk3*nk1 + 1
           nx_el(m,1)=n
           m = (j-1) + (i-1)*nk2 + (k-1)*nk1*nk2 + 1
           nx_el(m,2)=n
      ENDDO
     ENDDO
  ENDDO

  IF(nspin==2) THEN
     DO i=1,nks
        nx_el(i+nks,:)=nx_el(i,:)+nks
     ENDDO
  ENDIF
  l3dstring=.true.
!setup transfromation matrix
!  do i=1,3
!     cry_to_cart(:,i)=bg(:,i)
!     sca=sqrt(cry_to_cart(1,i)**2.d0+cry_to_cart(2,i)**2.d0+cry_to_cart(3,i)**2.d0)
!     cry_to_cart(:,i)=cry_to_cart(:,i)/sca
! enddo
! call  invmat (3, cry_to_cart, transform_el, sca)


  DO i=1,3
     sca=at(1,i)**2.d0+at(2,i)**2.d0+at(3,i)**2.d0
     sca=dsqrt(sca)
     bg_n(1:3,i)=(1.d0/sca)*at(1:3,i)
  ENDDO


  DO i=1,3
     DO j=1,3
        cry_to_cart(j,i)=bg_n(1,j)*bg_n(1,i)+bg_n(2,j)*bg_n(2,i)+bg_n(3,j)*bg_n(3,i)
     ENDDO
 ENDDO
 CALL  invmat (3, cry_to_cart, transform_el, sca)

!set up electric field

!calculate EFFECTIVE electric field on crystal axis
  efield_cry(:)=0.d0

!  do i=1,3
!     do j=1,3
!        efield_cry(i)=efield_cry(i)+transform_el(i,j)*efield_cart(j)
!     enddo
!  enddo

  DO i=1,3
    ! do j=1,3
        !efield_cry(i)=efield_cry(i)+transform_el(i,j)*(efield_cart(1)*bg_n(1,j)+efield_cart(2)*bg_n(2,j)+efield_cart(3)*bg_n(3,j))
        efield_cry(i)=efield_cry(i)+efield_cart(1)*bg_n(1,i)+efield_cart(2)*bg_n(2,i)+efield_cart(3)*bg_n(3,i)
  !enddo
  ENDDO
  !efield_cry(:)=0.001d0
  !efield_cry(3)=0.001d0
  WRITE(*,*) 'EFIELD CART', efield_cart(1),efield_cart(2), efield_cart(3)
  WRITE(*,*) 'EFIELD CRY', efield_cry(1),efield_cry(2), efield_cry(3)

  WRITE(*,*) 'BG1', bg(1,1),bg(2,1),bg(3,1)
  WRITE(*,*) 'BG1', at(1,1),at(2,1),at(3,1)
!
  WRITE(*,*) 'nx_el1', nx_el(1:nks,1)
!  write(*,*) 'nx_el2', nx_el(1:nks,2)
!  write(*,*) 'nx_el3', nx_el(1:nks,3)
   RETURN
 END SUBROUTINE kpoint_grid_efield
