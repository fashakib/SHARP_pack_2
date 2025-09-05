PROGRAM ring_polymer_surface_hopping

  use global_module
  use sysdef_module
  use modelvar_module
  use models_module
  use initial_module
  use nonadiabatic_module
  use propagation_module
  use runtraj_module
  use print_module
  use rpmd_module

  implicit none

  integer  :: i, icpu
  integer  :: seed_dimension,time

  real*8   :: RANF
  real*8   :: tt0,tt1

  logical  :: lkval

  integer, allocatable :: seed(:)

  call CPU_TIME(tt0)

! read input parameters from "input.in" file
  call sysdef(lkval)

! adjust nTraj based on nCPUs
  ntraj = int(real(ntraj)/real(ncpu))

  if(nequil .eq. 0)nequil=int(0.2*nsample)
  nskip = int((nsample-nequil)/ntraj)

! Initialize RNG seed
  call RANDOM_SEED(size=seed_dimension)
  allocate(seed(seed_dimension))
  do i=1,seed_dimension
    if(iseed .ne. 0)then
       seed(i) = iseed+3*(i-1)
    elseif(iseed .eq. 0)then
       seed(i) = time()+3*(i-1)
    endif
  enddo

  iseed=seed(1)
!  seed = time()
  if(ncpu >1)then
    open(1,file='icpu.in',status='old')
    read(1,*) icpu
    close(1)
    seed = seed + icpu*12345
  endif

!  seed = 123
  call RANDOM_SEED(PUT=seed)

! allocate variables and initialize values
  call modelallocat()

! assign some model specific parameters like mass, etc
  call modelParam(keymodel)

  call cmat_init()

  ! open file for writing purpose into the files
  call openfile()

!=======================initialization==================================

! print input parameters as "param.out"
  call printin()

! bath initialization for spin-boson
  if(keymodel==7)call bath_init(keybath)

! print analytical potential energy surface
  call printHel()

! main trajectory loop
  call runTraj()

! writing hopping statistics
  if(ldtl)call hopping_stat()
  call hopping_stat2()

  call CPU_TIME(tt1)
  if(ldtl)then
    write(nrite_dcoup,*) '#RUN TIME: ',(tt1-tt0), 'seconds'
    write(nrite_dcoup,*) '#RUN TIME: ',(tt1-tt0)/ntraj, 'seconds/traj'
  endif

! close all open files
  call close_file()

! print output results
  call printresults(lkval)

! deallocate the arrays
  call modeldeallocat()

END PROGRAM ring_polymer_surface_hopping

