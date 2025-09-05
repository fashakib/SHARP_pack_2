      module global_module
!**********************************************************************
!     SHARP Pack module for defining global variables
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
      real*8,parameter    :: hbar=1.d0
      complex*16,parameter:: eye=(0.,1.)

      integer             :: ndim=1
      integer             :: np     
      real*8              :: dt    
      real*8              :: dtq

      integer             :: keymodel
      integer             :: nstates
      integer             :: nsteps
      integer             :: nsample
      integer             :: ntraj
      integer             :: ntrajR
      integer             :: nb
      integer             :: iskip 
      integer             :: istrj 
      integer             :: v0key
      integer             :: R0key
      real*8              :: P0
      real*8              :: R0 
      real*8              :: omega
      character(len=36)   :: modelname
      character(len=8)    :: method
      character(len=8)    :: dyn

      logical             :: ldtl
      logical             :: lfft
  
      integer,parameter   :: nrite_hopp=110
      integer,parameter   :: nrite_dcoup=120
      integer,parameter   :: nrite_bath=140
!      integer,parameter   :: nrite_rgr=201
      integer,parameter   :: nrite_pimd=202
      integer,parameter   :: nrite_samp=203
      integer,parameter   :: nrite_psi=204
      integer,parameter   :: nrite_dboc=205
 
      integer             :: vrkey
      integer             :: nfrust_hop 
      integer             :: nfrust_hop2
      integer             :: ncpu 
      integer             :: vckey
      integer             :: iseed  !random number generator seed

! Conversion factors to the atomic units
! Energy from eV, Time from ps, Distance from angstrom
! Temperature from K, Frequency from cm^-1

      real*8,parameter    :: energy=27.21d0,tim=0.00002419d0,dist=0.5291d0
      real*8,parameter    :: temp=315774d0,freq=219463.343d0

! Classical environment parameters
      real*8,parameter    :: pi=3.14159265d0  !DACOS(-1.d0)
      real*8              :: beta
!      real*8              :: mp 
      real*8              :: kT

!! Spin-Boson parameters
      integer           :: keybath
      real*8            :: eps, delta, enu
      real*8            :: E_r,wc,wmax
      real*8, parameter :: wmb = 3.0d0   ! wm = wm/wc ~ 3 

!!! LinearChainModel !!!
      integer :: dlevel
      integer :: nprint
      real*8  :: gamaLC(64)
      real*8  :: sigmaLC(64)
      real*8  :: d_ab(3,3)     
      real*8  :: v11
      real*8  :: v22
      real*8,parameter :: kJ_mol2au=0.00038125d0
      real*8,parameter :: amu2au=1822.4d0

!!! with 3-states LinearChainModel !!!
      real*8  :: v33

!!! added in v2 for SB model 
      logical :: lpimd
      logical :: ltrpmd
      integer :: mdkey
      real(8) :: taut


!!! added for bead approximation      
      integer :: apkey

!!! added more keywords for bead sampling      
      logical :: lnorm
      logical :: ldsurf

!!! added for decoherence       
      integer :: dckey
      integer :: nequil
      integer :: nskip 

!!! added for diagonal born-oppenheimer correction (2nd derivatice nac)       
      logical :: ldboc
      logical :: lpsi2
      integer :: dvkey
      real*8,parameter :: eps_dboc=0.1d0

      end module global_module
