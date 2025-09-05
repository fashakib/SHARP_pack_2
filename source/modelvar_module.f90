      module modelvar_module
!**********************************************************************
!     SHARP PACK module that contains array variables     
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
      use global_module, only: np, nb, nstates, nsteps, nprint,iskip

      implicit none

      real*8, allocatable    :: redmat(:,:)
      complex*16,allocatable :: redmat_ec(:,:,:)
      complex*16,allocatable :: diabat2(:,:)
      real*8, allocatable    :: diabat1(:,:),diabat3(:,:)
      real*8,allocatable     :: rc(:),vc(:)
      real*8,allocatable     :: rp_samp(:,:),vp_samp(:,:)
      real*8,allocatable     :: rp(:,:),vp(:,:),fp(:,:)
      real*8,allocatable     :: dhel(:,:,:)
      real*8,allocatable     :: dhel_rc(:,:,:)
      real*8,allocatable     :: dev_rc_old(:,:)
      real*8,allocatable     :: vdotd_old(:,:), vdotd_new(:,:)
      real*8,allocatable     :: vdotd(:,:)
      real*8,allocatable     :: devc_rc(:,:,:)
      real*8,allocatable     :: deva(:)
      real*8,allocatable     :: psi(:,:,:) ! old,current or new electronic step
      real*8,allocatable     :: hel(:,:)   
      complex*16,allocatable :: adia(:)
      real*8,allocatable     :: eva(:,:)

      real*8, allocatable    :: redmatR(:,:)
      real*8,allocatable     :: diabat1R(:,:),diabat3R(:,:)
      complex*16,allocatable :: diabat2R(:,:)
      complex*16,allocatable :: redmat_ecR(:,:,:)

      integer,allocatable    :: nJump(:,:), nJumpFail(:,:) 
      integer,allocatable    :: nFrust(:,:), nFrustR(:,:) 
      integer,allocatable    :: nIniStat(:)

      real*8,allocatable     :: cmat(:,:)
      real*8,allocatable     :: mp(:,:)
      real*8,allocatable     :: nfreq(:)

      ! sytem-bath variable for spin-boson model
      real*8,allocatable     :: Wb(:)
      real*8,allocatable     :: cb(:)

      real*8,allocatable     :: pileC1(:),pileC2(:)
      real*8,allocatable     :: rp0(:)

      real*8,allocatable     :: psi_bead(:,:,:,:) ! old,current or new electronic step
      real*8,allocatable     :: psi_dboc(:,:,:,:) ! old,current or new electronic step
      real*8,allocatable     :: d_12(:,:,:)
      real*8,allocatable     :: dij_bead(:,:,:,:)
      contains

      subroutine modelallocat()
!**********************************************************************
!     SHARP PACK subroutine to allocate array variables     
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none
     
      !nprint = nsteps 
      !if(nsteps .gt. 50000)
      nprint = int(nsteps/iskip)

      allocate(redmat(nstates,0:nprint))
      allocate(redmat_ec(nstates,nstates,0:nprint))
      allocate(diabat2(nstates,0:nprint))
      allocate(diabat1(nstates,0:nprint),diabat3(nstates,0:nprint))
      allocate(rc(np),vc(np))
      allocate(rp_samp(np,nb),vp_samp(np,nb))
      allocate(rp(np,nb),vp(np,nb),fp(np,nb))
      allocate(dhel(nstates,nstates,np))
      allocate(dhel_rc(nstates,nstates,np))
      allocate(dev_rc_old(nstates,nstates))
      allocate(vdotd_old(nstates,nstates),vdotd_new(nstates,nstates))
      allocate(vdotd(nstates,nstates))
      allocate(devc_rc(nstates,nstates,np))
      allocate(deva(nstates))
      allocate(psi(nstates,nstates,2)) ! old,current or new electronic step
      allocate(hel(nstates,nstates))   
      allocate(adia(nstates))
      allocate(eva(nstates,2))

      allocate(redmatR(nstates,0:nprint))
      allocate(diabat2R(nstates,0:nprint))
      allocate(diabat1R(nstates,0:nprint),diabat3R(nstates,0:nprint))
      allocate(redmat_ecR(nstates,nstates,0:nprint))

      allocate(nJump(nstates,nstates),nJumpFail(nstates,nstates)) 
      allocate(nFrust(nstates,nstates),nFrustR(nstates,nstates)) 
      allocate(nIniStat(nstates))

      allocate(cmat(nb,nb))
      allocate(mp(np,nb))
      allocate(nfreq(nb))

      allocate(Wb(np))
      allocate(cb(np))

      allocate(pileC1(nb),pileC2(nb))
      
      allocate(rp0(np))
      
      allocate(psi_bead(nstates,nstates,2,nb)) ! old,current or new electronic step
      allocate(psi_dboc(nstates,nstates,2,nb)) ! old,current or new electronic step
      allocate(d_12(nstates,nstates,np))
      allocate(dij_bead(nstates,nstates,np,nb))

      diabat1=0.
      diabat2=cmplx(0.,0.)
      diabat3=0.

      diabat1R=0.
      diabat2R=cmplx(0.,0.)
      diabat3R=0.

      redmat=0.
      redmat_ec=cmplx(0.,0.)

      redmatR=0.
      redmat_ecR=cmplx(0.,0.)

      nJump = 0
      nJumpFail = 0
      nFrust = 0
      nFrustR = 0
      nIniStat = 0

      cmat = 0.d0
      mp = 1.d0

      psi_bead = 0.d0
      psi_dboc = 0.d0
      d_12 = 0.d0
      dij_bead = 0.d0

      return 

      end subroutine modelallocat


      subroutine modeldeallocat()
!**********************************************************************
!     SHARP PACK subroutine to deallocate array variables     
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none

      deallocate(redmat)
      deallocate(redmat_ec)
      deallocate(diabat2)
      deallocate(diabat1,diabat3)
      deallocate(rc,vc)
      deallocate(rp_samp,vp_samp)
      deallocate(rp,vp,fp)
      deallocate(dhel)
      deallocate(dhel_rc)
      deallocate(dev_rc_old)
      deallocate(vdotd_old, vdotd_new)
      deallocate(vdotd)
      deallocate(devc_rc)
      deallocate(deva)
      deallocate(psi) 
      deallocate(hel)   
      deallocate(adia)
      deallocate(eva)

      deallocate(redmatR)
      deallocate(diabat2R)
      deallocate(diabat1R,diabat3R)
      deallocate(redmat_ecR)

      deallocate(nJump,nJumpFail)
      deallocate(nFrust,nFrustR)
      deallocate(nIniStat)

      deallocate(cmat)
      deallocate(mp)
      deallocate(nfreq)

      deallocate(Wb)
      deallocate(cb)

      deallocate(pileC1,pileC2)
      
      deallocate(rp0)

      deallocate(psi_bead)
      deallocate(psi_dboc)
      deallocate(d_12)
      deallocate(dij_bead)

      return 

      end subroutine modeldeallocat

!**********************************************************************

      end module modelvar_module
