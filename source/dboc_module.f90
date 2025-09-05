  module dboc_module
!**********************************************************************
!     SHARP Pack module for diagonal born-oppenheimer correction 
!     ref: JCP 144, 154103 (2016)
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
!  use global_module
  use global_module, only : nstates,lpsi2,nb,eps_dboc
  use modelvar_module,only : psi_dboc,d_12,dij_bead
  use models_module

  implicit none

  private

  public :: compute_nac, compute_grad_dboc

  contains

  subroutine compute_nac(nat,ib,r,d_ij,lpsi)
!**********************************************************************
!     SHARP Pack routine to compute NACVs, d_I_ij
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
  implicit none

  integer            :: i,j,k,l,ip
  integer,intent(in) :: nat,ib
  real*8,intent(in)  :: r(nat)
  real*8,intent(out) :: d_ij(nstates,nstates,nat)
  real*8             :: hmat(nstates,nstates)
  real*8             :: enr(nstates)
  real*8             :: delH(nstates,nstates,nat)
  logical            :: lpsi    !if true, save psi for phase correction

  call gethel(r,hmat,delH)
  call Diag(enr,psi_dboc(:,:,2,ib),hmat)

! phase correction of wavefunction  
  do i = 1, nstates
    psi_dboc(:,i,2,ib) = psi_dboc(:,i,2,ib) * &
             dot_product(psi_dboc(:,i,1,ib),psi_dboc(:,i,2,ib))/ &
             abs(dot_product(psi_dboc(:,i,1,ib),psi_dboc(:,i,2,ib)))
  end do

  d_ij=0.d0

  do i=1,nstates-1
    do j=i+1,nstates
      do ip=1,nat
        !if (i.ne.j) then
        !  do k=1,nstates
        !    do l=1,nstates
        !      d_ij(ip) = d_ij(ip)+psi_dboc(k,i)*delH(k,l,ip)*psi_dboc(l,j)
        !    enddo
        !  end do
        !  d_ij(ip) = d_ij(ip)/enr(j)-enr(i)
        !endif
        d_ij(i,j,ip) = sum(psi_dboc(:,i,2,ib)*matmul(delH(:,:,ip),psi_dboc(:,j,2,ib)))
      end do
      d_ij(i,j,:) = d_ij(i,j,:)/(enr(j)-enr(i))
      d_ij(j,i,:) = -d_ij(i,j,:)
    end do
  enddo

  if(lpsi)psi_dboc(:,:,1,ib) = psi_dboc(:,:,2,ib)

  end subroutine compute_nac

  subroutine compute_grad_dboc(nat,ib,r,mass,istate,grad)
!**********************************************************************
!     SHARP Pack routine to compute dboc corrected force by finite
!     difference method
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
  implicit none

  integer            :: j,ip
  integer,intent(in) :: nat,ib,istate
  real*8             :: r(nat)
  real*8,intent(in)  :: mass(nat)
  real*8,intent(out) :: grad(nat)
  real*8             :: rplus(nat),rminus(nat)
  real*8             :: dij(nstates,nstates,nat)
  real*8             :: dijplus(nstates,nstates,nat)
  real*8             :: dijminus(nstates,nstates,nat)
  real*8             :: eps
  logical            :: lpsi

!  eps =0.1d0  !
  eps = eps_dboc
  grad = 0.d0
  rplus = r + eps
  rminus = r - eps

  lpsi = .true.
  call compute_nac(nat,ib,r,dij,lpsi)
! save NAC of current step  
  d_12 = d_12 + dij

! save NAC of current beads
  dij_bead(:,:,:,ib) = dij(:,:,:)

  lpsi=.false.
  call compute_nac(nat,ib,rplus,dijplus,lpsi)

  lpsi=.false.
  call compute_nac(nat,ib,rminus,dijminus,lpsi)

  do j = 1, nstates

    if(j .ne. istate)then
      do ip = 1, nat

        grad(ip) = grad(ip) - 1.d0/mass(ip) * dij(istate,j,ip) * &
                (dijplus(istate,j,ip)-dijminus(istate,j,ip))/(2.d0*eps)

      enddo
    endif
  enddo

  end subroutine compute_grad_dboc

  end module dboc_module
