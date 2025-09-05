  module runtraj_module
!**********************************************************************
!     SHARP Pack module for running trajectories
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
  use global_module
  use modelvar_module
  use models_module
  use initial_module
  use rpmd_module
  use nonadiabatic_module
  use propagation_module
  use print_module
  use dboc_module

  implicit none

  contains

  subroutine runTraj()
!**********************************************************************
!     SHARP Pack routine to run trajectories
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
  implicit none

  integer  :: i,j,k,l,itime,itraj,ibd,ip
  integer  :: istate, inext
  integer  :: nlit
  integer  :: iprint
  integer  :: ia
  
  real*8   :: rnlit
  real*8   :: RANF
  real*8   :: crit,length,currentlength
  real*8   :: KE,PE,Ering,TotE,Vn
  real*8   :: t2au=41340.d0
  real*8,allocatable :: db2pop(:,:)
  
  integer,parameter  :: nrite_dtl1=1001
  integer,parameter  :: nrite_dtl2=1002

  real*8  :: systmp,syctmp,rgr

  real*8   :: rp_pimd(np,nb,ntraj),vp_pimd(np,nb,ntraj)
  real*8   :: rp0samp(np,nb,ntraj),vp0samp(np,nb,ntraj)
  real*8   :: eva_avg

! bead approximation 
  real*8   :: eva_bead(nstates,2,nb)!, psi_bead(nstates,nstates,2,nb)
  real*8   :: dhel_bead(nstates,nstates,np,nb)
  real*8   :: vdotd_bead(nstates,nstates,nb)
  real*8   :: norm

  real*8   :: rc0

! diagonal born-oppenheimer corretion JCP 144, 154103 (2016)
  real*8   :: dboc(nstates)
  real*8   :: enr_dboc(nstates)
  logical  :: lpsi

  ntrajR = 0
 
  ia = floor((nsteps*dt)/t2au)+1
  allocate(db2pop(nstates,ia))
  db2pop = 0.d0
 
  nlit=int(dt/dtq)
  rnlit = 1.0/REAL(nlit,8)
  !nequil = int(0.2*nsample)

  rp0samp = 0.d0
  vp0samp = 0.d0

  if(dlevel .eq. 1)then
     open(nrite_dtl1,file='HISTORYc',status='unknown')
     write(nrite_dtl1,'(1x,A)') '# nTraj       Time(au)    rc     vc '
  endif
  if(dlevel .eq. 2)then
     open(nrite_dtl2,file='HISTORYb',status='unknown')
     write(nrite_dtl2,'(1x,A)') '# nTraj       Time(au)    rp(np,nb)     vp(np,nb) '
  endif

  if(lnorm .and. R0key.ne.1)then
     write(0,*) 'rmode should be direct!!'
     stop
  endif

! main trajectory loop
  DO itraj=1,ntraj
     
!  sample initial distribution of position and momenta
    if(.not.lpimd)then
      call sample_init(vp_samp,rp_samp)

!  transform postions normal to coordinates system
      if(lnorm .and. R0key.eq.1) call realft(rp_samp,np,nb,-1)
 
      if((keymodel==12).or.(keymodel==13))call sample_init_lchain(vp_samp,rp_samp)

!  sample initial distribution pimd run for position and momenta
    elseif(lpimd)then

!  pimd sampling from single run trajectory
      if(itraj .eq. 1)then

        call sample_pimd(itraj,rp_pimd,vp_pimd)

      endif

      rp_samp = rp_pimd(:,:,itraj)
      vp_samp = vp_pimd(:,:,itraj)

    endif

! reset centroid of ring polymer to dividing surface
    if(ldsurf)then
      do ip = 1, np
        rc0 = sum(rp_samp(ip,:))/real(nb)
        rp_samp(ip,:) = rp_samp(ip,:) - rc0
      enddo
    endif

    rp = rp_samp
    vp = vp_samp

!   store sampling pos,vel to compute initial distribution
    rp0samp(:,:,itraj) = rp_samp
    vp0samp(:,:,itraj) = vp_samp

! Caculating the initial centroid variables

!   centroid approximation
!    if(apkey == 1)then

    rc=0.d0
    vc=0.d0

    do ibd=1,nb
      rc(:)=rc(:)+rp(:,ibd)/real(nb)
      vc(:)=vc(:)+vp(:,ibd)/real(nb)
    enddo

! Initialize position and momentum of solvent

! get the intial wavefunction 

      call gethel(rc(:),hel(:,:),dhel_rc(:,:,:))
      call Diag(eva(:,1),psi(:,:,1),hel(:,:))

    if(apkey == 1)then
      call compute_vdotd_old(vdotd_old,psi(:,:,1),vc,eva(:,1),dhel_rc)

!   bead approximation
    elseif(apkey == 2)then

      eva(:,1) = 0.d0
      !psi(:,:,1) = 0.d0
      vdotd_old = 0.d0

      do ibd=1,nb

        call gethel(rp(:,ibd),hel(:,:),dhel_bead(:,:,:,ibd))
        call Diag(eva_bead(:,1,ibd),psi_bead(:,:,1,ibd),hel(:,:))

        call compute_vdotd_old(vdotd_bead(:,:,ibd),psi_bead(:,:,1,ibd),vp(:,ibd),eva_bead(:,1,ibd),dhel_bead(:,:,:,ibd))

        eva(:,1) = eva(:,1) + eva_bead(:,1,ibd)
       ! psi(:,:,1) = psi(:,:,1) + psi_bead(:,:,1,ibd)
        vdotd_old = vdotd_old + vdotd_bead(:,:,ibd)

      enddo

      eva(:,1) = eva(:,1)/real(nb)
      vdotd_old = vdotd_old/real(nb)

    endif

    adia(:)=CMPLX(psi(:,1,1),0.)

    CALL RANDOM_NUMBER(RANF)

    crit=RANF

    currentlength=0.d0
    do i=1,nstates
      currentlength=currentlength+psi(i,1,1)**2
      if (crit.le.currentlength) then
        exit
      else
        continue
      endif
    enddo
    istate=i
    inext=i
    nIniStat(istate) = nIniStat(istate) + 1

!=================================================================================
! Initialize the forces

    do ibd=1,nb
      call gethel(rp(:,ibd),hel(:,:),dhel_bead(:,:,:,ibd))
      call Diag(eva_bead(:,1,ibd),psi_bead(:,:,1,ibd),hel(:,:))

      psi_dboc(:,:,1,ibd) = psi_bead(:,:,1,ibd)

      if((keymodel==12).or.(keymodel==13))then
        CALL FORCE_Lchain(rp(1:np,ibd), fp(1:np,ibd))

      else
        CALL FORCE(psi_bead(:,:,1,ibd),istate,dhel_bead(:,:,1:np,ibd), fp(1:np,ibd))

      endif
    enddo

!==================================================================================

! start a trajectory for n-steps     
    DO itime = 1, NSTEPS
      
      d_12 = 0.d0  
! main nuclear propagation step
      CALL ADVANCE_MD(istate,rp,vp,fp,eva_bead,dhel_bead)
       
! Caculating the centroid variables at each time step

!      if(apkey == 1)then
      rc=0.d0
      vc=0.d0
   
      do ibd=1,nb
        rc(:)=rc(:)+rp(:,ibd) !/real(nb)
        vc(:)=vc(:)+vp(:,ibd) !/real(nb)
      enddo
      rc(:)=rc(:)/real(nb)
      vc(:)=vc(:)/real(nb)

      CALL gethel(rc(:),hel(:,:),dhel_rc(:,:,:))
      ! get the intial wavefunction 
      CALL Diag(eva(:,2),psi(:,:,2),hel(:,:))

      do i = 1, nstates
        psi(:,i,2) = psi(:,i,2)*dot_product(psi(:,i,1),psi(:,i,2))/ &
                     abs(dot_product(psi(:,i,1),psi(:,i,2)))
      end do

!     centroid approximation
      if(apkey == 1)then

        call compute_vdotd(vdotd_new,psi)

!     bead approximation
      elseif(apkey == 2)then

        eva(:,2) = 0.d0
        vdotd_new = 0.d0

        do ibd = 1, nb

          do i = 1, nstates
            psi_bead(:,i,2,ibd) = psi_bead(:,i,2,ibd) * &
                    dot_product(psi_bead(:,i,1,ibd),psi_bead(:,i,2,ibd))/ &
                    abs(dot_product(psi_bead(:,i,1,ibd),psi_bead(:,i,2,ibd)))
          end do

          call compute_vdotd(vdotd_bead(:,:,ibd),psi_bead(:,:,:,ibd))

          eva(:,2) = eva(:,2) + eva_bead(:,2,ibd)
          vdotd_new = vdotd_new + vdotd_bead(:,:,ibd)


        end do

        eva(:,2) = eva(:,2)/real(nb)
        vdotd_new = vdotd_new/real(nb)

      endif

!  dboc implementation
     if(ldboc)then

      if(dvkey .eq. 2) d_12 = d_12/real(nb)

      dboc = 0.d0
      do i = 1, nstates
        do j = 1, nstates
          if(i .ne. j)then
              
             dboc(i) = dboc(i) + sum(0.5d0*hbar*hbar/mp(:,1)*d_12(i,j,:)*d_12(i,j,:))

          endif

        enddo
      enddo

      if(((ntraj .lt.10) .or. (mod(itraj,istrj).eq.0)).and.(mod(itime,iskip).eq.0))then
        write(nrite_dboc,101) itraj,itime*dt,rc(1),vc(1),eva(:,2),dboc(:),eva(:,2)+dboc(:)
      endif

     endif

      deva(:) = eva(:,2) - eva(:,1)  
       
      vdotd = vdotd_new - vdotd_old

      ! electronic time steps
      DO j = 1, nlit
        eva(:,2) = eva(:,1) + deva * rnlit
        vdotd_new = vdotd_old + vdotd * rnlit 

        ! advance the adiabatic wave fucntion
        CALL ADVANCE_ADIA(eva(:,:),vdotd_old(:,:),vdotd_new(:,:),adia(:))

        CALL Hopping(istate,inext,vdotd_new(:,:),adia(:))         
           
        ! give the current to previous
        eva(:,1) = eva(:,2)      !eva(now) give to eva(old)
        vdotd_old = vdotd_new
      END DO

      psi(:,:,1) = psi(:,:,2)
      psi_bead(:,:,1,:) = psi_bead(:,:,2,:)
      eva_bead(:,1,:) = eva_bead(:,2,:)

!     add decoherence scheme(s)
      if(dckey == 1)call decoh_damp(adia,istate,eva(:,2))

      IF (istate.NE.inext) THEN
      ! total potential energy difference for the whole ring polymer hamiltonian 

!     rpsh-ca method of energy/velocity rescaling
        if(apkey == 1)then 
          CALL EKINRESCALE(istate,inext,vc,vp,rp,eva(:,:),eva_bead(:,2,:),dhel_rc(:,:,:),psi(:,:,:),itime)    
!     rpsh-ba method of energy/velocity rescaling      
        elseif(apkey == 2)then
          CALL EKINRESCALE_BEAD(istate,inext,vp,rp,eva(:,:),eva_bead(:,2,:),dhel_bead,psi_bead(:,:,2,:),itime)    
        endif 
      END IF

!=============store the adiabatic/diabatic density matrix by differnet ways============
      if((mod(itime,iskip).eq.0).or.(itime==1))then
        iprint = int(itime/iskip)
        CALL  pop_estimator(istate,iprint,diabat1,diabat2,diabat3,redmat,redmat_ec)
      
! reflected probability calculation for tully model
        if(keymodel <= 3)then
          if((vc(1)<0).and.(rc(1) < 0))then
            CALL  pop_estimator(istate,iprint,diabat1R,diabat2R,diabat3R,redmatR,redmat_ecR)
          endif
        endif
      endif

! db2pop calculation
      ! bined coefficient every 1 ps 
     if((keymodel==12).or.(keymodel==13))then
       ia = floor((itime*dt)/t2au)+1
       db2pop(istate,ia)=db2pop(istate,ia)+1.
     endif

     ! printdetail()
     IF(ldtl)THEN
       ! rc=0.d0
       vc=0.d0
       do ibd=1,nb
       !  rc(:)=rc(:)+rp(:,ibd) !/real(nb)
         vc(:)=vc(:)+vp(:,ibd) !/real(nb)
       enddo
       !rc(:)=rc(:)/real(nb)
       vc(:)=vc(:)/real(nb)

       if(((ntraj .lt.10) .or. (mod(itraj,istrj).eq.0)).and.(mod(itime,iskip).eq.0))then
         call calEnergy(rp,vp,KE,Ering,Vn,TotE,istate,systmp,eva_bead(:,2,:))
         syctmp = sum(vc*vc*mp(:,1))/np*temp
!  compute radius of gyration (debugging mode)
         call calcRg(itraj,itime*dt,rp,rgr)

         write(nrite_dcoup,101) itraj,itime*dt,rc(1),vc(1),KE,Ering,Vn,TotE, &
                  real(istate),eva(istate,1),vdotd_new,systmp,syctmp,rgr

         if(dlevel .eq. 1) write(nrite_dtl1,102) itraj,itime*dt,rc,vc
         if(dlevel .eq. 2) write(nrite_dtl2,102) itraj,itime*dt,(rp(ip,:),ip=1,np),(vp(ip,:),ip=1,np)

         write(nrite_psi,601) itraj,itime*dt,istate,rc(1),psi(:,:,2),adia,eva(:,2),vdotd_new(1,2)

       endif
     ENDIF

! nuc step loop
    ENDDO
    if((vc(1)<0).and.(rc(1)<0)) ntrajR = ntrajR + 1

   ! traj loop
  ENDDO

! print initial sampling distribution
  call printInitSamp(rp0samp,vp0samp)

 101 format(i10,1x,f15.3,1x,150(e18.8E3,2x))
 102 format(i10,1x,f15.3,1x,350(f15.6,2x))
 601 format(i10,1x,f15.3,1x,i6,20(1x,e20.6))
 700 format(i8,1x,f12.3,1x,1000(f12.4,2x))
   
  if(dlevel .eq. 1) close(nrite_dtl1)
  if(dlevel .eq. 2) close(nrite_dtl2)
  
  if((keymodel==12).or.(keymodel==13))then
    open(1,file='db2pop.out')
    write(1,*) '# db2pop'
    write(1,*) '#time(ps) populations'
    do i = 1, ia-1
      write(1,'(5(e13.5E3,2x))') real(i), (db2pop(j,i)/(ntraj*t2au/dt),j=1,nstates)
    enddo
   close(1)
  endif

  deallocate(db2pop)

  end subroutine


  subroutine compute_vdotd_old(vdotd,psi,vp,eva,delH)
!**********************************************************************
!     SHARP Pack routine to compute vdotd_old
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
  implicit none

  integer            :: i,j,k,l,ip
  real*8,intent(out) :: vdotd(nstates,nstates)
  real*8,intent(in)  :: psi(nstates,nstates)
  real*8,intent(in)  :: vp(np)
  real*8,intent(in)  :: eva(nstates)
  real*8,intent(in)  :: delH(nstates,nstates,np)

  vdotd=0.d0

  if((keymodel==12).or.(keymodel==13))then
    do k=1,nstates
      do l=1,nstates
        if(k.ne.l)then
          do ip=1,1  !np
            vdotd(k,l) = vdotd(k,l)+vp(ip)*d_ab(k,l)
          enddo
        endif
      enddo
    enddo

  else
    do k=1,nstates
      do l=1,nstates
        if (k.ne.l) then
          do i=1,nstates
            do j=1,nstates
              do ip=1,np
                vdotd(k,l) = vdotd(k,l)+psi(i,k)* &
                        delH(i,j,ip)*psi(j,l)*vp(ip)/(eva(l)-eva(k))
              enddo
            enddo
          end do
        endif
      end do
    end do
  endif

  end subroutine compute_vdotd_old


  subroutine compute_vdotd(vdotd,psi)
!**********************************************************************
!     SHARP Pack routine to compute vdotd_new
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
  implicit none

  integer            :: i,j,k,l,ip
  real*8,intent(out) :: vdotd(nstates,nstates)
  real*8,intent(in)  :: psi(nstates,nstates,2) ! old,current or new electronic step

  ! middle decoup based on SHS-Tully 94

  vdotd=0.d0

! CHECK using eq 32\times velocity Tully94 paper
  do k=1,nstates
    do l=1,nstates
      if(k.ne.l) then
        if((keymodel==12).or.(keymodel==13))then
        !!v*d_ij = p1*d_ij/m as for LinearChainModel
          vdotd(k,l) = vc(1)*d_ab(k,l)
        else
          vdotd(k,l) = sum(psi(:,k,1)*psi(:,l,2)-psi(:,k,2)*psi(:,l,1))/(2.d0*dt)
        endif
      end if
    enddo
  enddo

  end subroutine compute_vdotd
        

  subroutine calEnergy(rp,vp,KE,Erng,Vn,Hn,istate,systmp,ener_bead)
!**********************************************************************
!     SHARP Pack routine to calculate energies at particular time
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
!  use global_module
  implicit none

  integer            :: ip,ibd
  integer,intent(in) :: istate

  real*8,intent(in)  :: rp(np,nb),vp(np,nb)
  real*8,intent(in)  :: ener_bead(nstates,nb)
  real*8,intent(out) :: systmp
  real*8             :: Wn2
  real*8             :: KE,Vn,Erng,Hn

  wn2 = nb/(beta*hbar)
  wn2 = wn2*wn2

  KE = 0.5d0*sum(vp*vp*mp)

! get system temperature for relation: KE = 0.5*kB*T/Ndof  
  systmp = 2.d0*KE/(np*nb)*temp
  !syctmp = sum(vc*vc*mp(:,1))/np*temp

!  PE = nb*E_i
  ! calculate active energy of each bead
  Vn = 0.d0
  Vn = sum(ener_bead(istate,:))

  Erng = 0.d0
  do ip = 1, np
    do ibd = 1, nb
      if(ibd .eq. 1)then
        Erng = Erng + mp(ip,ibd)*(rp(ip,ibd) - rp(ip,nb))**2
      else  
        Erng = Erng + mp(ip,ibd)*(rp(ip,ibd) - rp(ip,ibd-1))**2
      endif
    enddo
  enddo

  Erng = 0.5d0*wn2*Erng

  Hn = KE+Erng+Vn

  end subroutine calEnergy
 

  subroutine printHel()
!**********************************************************************
!     SHARP Pack routine to print analytical energy surface
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
!  use global_module, only : np,nstates
  use modelvar_module, only : dhel_rc,eva
   
  implicit none
  
  integer :: i,j,k
  integer :: keymodel

  real*8  :: x(np)
  real*8  :: hel(nstates,nstates)!,eva(nstates,2)
  real*8  :: d_ab(nstates,nstates),psi(nstates,nstates,2)!,dhel_rc(nstates,nstates,np)

  open(1,file='energy_surface.out',status='unknown')
  write(1,'(A13,2x,A)') '# MODEL:    ', modelname
  write(1,'(A)') ' # R[-12:12], E_dia(n,n), E_adia(n), nondaibatic_coupling_d_ab(n,n)'

  vc(:) = 1.d0

  do k = 1,961
     x(:) = -12+(k-1)*0.025

     call gethel(x(:),hel(:,:),dhel_rc(:,:,:))
     call Diag(eva(:,1),psi(:,:,2),hel(:,:))

     if(k == 1)psi(:,:,1)=psi(:,:,2)

      do i = 1, nstates
        psi(:,i,2) = psi(:,i,2)*dot_product(psi(:,i,1),psi(:,i,2))/ &
                     abs(dot_product(psi(:,i,1),psi(:,i,2)))
      end do

     !call compute_vdotd_old(d_ab,psi)
      call compute_vdotd_old(d_ab,psi(:,:,2),vc,eva(:,1),dhel_rc)

     write(1,'(150e15.6e3)') x(1),hel,eva(:,1),d_ab

     psi(:,:,1)=psi(:,:,2)
   enddo

   close(1)

   end subroutine printHel

  subroutine calcRg(itraj,tstep,rx,rgr)
!**********************************************************************
!     SHARP Pack routine to comppute and print radius of gyration 
!     of ring polymer
!
!     date:: 03/20/2025

!     input::
!       rx : ring polymer coordinates
!
!     output::
!       rg : radius of gyration of ring polymer
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
  use global_module, only : np,nb
   
  implicit none
  
  integer :: ip,ib

  integer,intent(in) :: itraj
  real*8,intent(in)  :: tstep
  real*8,intent(in)  :: rx(np,nb)
  real*8 :: rav, rg(np),rgr

! compute center of mass of each ring polymer
  rgr = 0.d0
  do ip = 1, np
    rav = 0.d0
!    rav = sum(mp(ip,1:nb)*rx(ip,1:nb))/sum(mp(ip,1:nb))
    rav = sum(rx(ip,1:nb))/nb

!    rg(ip) = sum(mp(ip,1:nb)*(rx(ip,1:nb)-rav)**2)/sum(mp(ip,1:nb))
    rg(ip) = sum((rx(ip,1:nb)-rav)**2)/nb
    rg(ip) = sqrt(rg(ip))

    rgr = rgr + rg(ip)
  enddo

  rgr = rgr/np

  end subroutine calcRg

  subroutine printInitSamp(rp0samp,vp0samp)
!**********************************************************************
!     SHARP Pack routine to compute and print initial sampling 
!     distribution position and velocity 
!
!     date:: 04/11/2025

!     input::
!       rpsamp : coordinates of ring polymer
!       vpsamp : velocity or ring polymer
!     output::
!       histrogram distribution of ring polymer
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
  use global_module, only : np,nb,ntraj,nrite_samp
  implicit none

  real*8,intent(in) :: rp0samp(np,nb,ntraj)
  real*8,intent(in) :: vp0samp(np,nb,ntraj)

  integer :: i,j,k,bin
  integer :: ip,ib,half_window
  integer,parameter :: nbin=50,window_size=5
  real*8  :: x,xmin,xmax,dx,zero
  real*8  :: xbin(nbin),hist(nbin),shist(nbin)

  half_window = (window_size - 1)/2
  zero = 0.0

  ! write position distribution for n trajectories
  if(ntraj .gt. 1)then
  if((nb .gt. 1).or.((nb .eq. 1).and.(.not.ldsurf)))then
  open(nrite_samp,file='sampling_pos.out',status='unknown')
  do ip = 1, np
    do ib = 1, nb
       xmin = minval(rp0samp(ip,ib,1:ntraj))
       xmax = maxval(rp0samp(ip,ib,1:ntraj))
       dx = (xmax-xmin)/nbin
      
       hist(:) = 0.0
       do i = 1, ntraj
         bin = int((rp0samp(ip,ib,i)-xmin)/dx)+1
         if(bin .le. nbin)hist(bin)=hist(bin)+1.0
       enddo

       hist(:) = hist(:)/ntraj
       shist(:) = 0.0

       do i = 1, nbin
         xbin(i) = xmin+(i-0.5)*dx
         k = 0
         do j = max(1,i-half_window), min(nbin,i+half_window)
           shist(i) = shist(i) + hist(j)
           k = k +1
         enddo
         if(k > 0) shist(i) = shist(i)/real(k)
       enddo

       write(nrite_samp,102)'# X0sample:np-nb :',xmin,xmax,ip,ib
       write(nrite_samp,101) xbin(1)-0.5*dx,zero,zero,ip,ib
       do i = 1, nbin
         write(nrite_samp,101) xbin(i),hist(i),shist(i),ip,ib
       enddo
       write(nrite_samp,101) xbin(nbin)+0.5*dx,zero,zero,ip,ib
       write(nrite_samp,*)

    enddo
  enddo
  close(nrite_samp)
  endif

  ! write velocity distribution if vsamp is not fixed
  if(V0key .ne. 0)then
  open(nrite_samp,file='sampling_vel.out',status='unknown')
  do ip = 1, np
    do ib = 1, nb
       xmin = minval(vp0samp(ip,ib,:))
       xmax = maxval(vp0samp(ip,ib,:))
       dx = (xmax-xmin)/nbin
       
       hist(:) = 0.0
       do i = 1, ntraj
         bin = int((vp0samp(ip,ib,i)-xmin)/dx)+1
         if(bin .le. nbin)hist(bin)=hist(bin)+1.0
       enddo

       hist(:) = hist(:)/ntraj
       shist(:) = 0.0

       do i = 1, nbin
         xbin(i) = xmin+(i-0.5)*dx
         k = 0
         do j = max(1,i-half_window), min(nbin,i+half_window)
           shist(i) = shist(i) + hist(j)
           k = k +1
         enddo
         if(k > 0) shist(i) = shist(i)/real(k)
       enddo

       write(nrite_samp,102)'# X0sample:np-nb :',xmin,xmax,ip,ib
       write(nrite_samp,101) xbin(1)-0.5*dx,zero,zero,ip,ib
       do i = 1, nbin
         write(nrite_samp,101) xbin(i),hist(i),shist(i),ip,ib
       enddo
       write(nrite_samp,101) xbin(nbin)+0.5*dx,zero,zero,ip,ib
       write(nrite_samp,*)

    enddo
  enddo
  close(nrite_samp)

! Maxwell-Boltzmann distribution of speeds  
  open(nrite_samp,file='sampling_dis.out',status='unknown')
  xmin = 0.0 
  xmax = maxval(abs(vp0samp(:,:,:)))
  dx = (xmax-xmin)/nbin
       
  hist(:) = 0.0
  do ip = 1, np
    do ib = 1, nb
       do i = 1, ntraj
         bin = int((abs(vp0samp(ip,ib,i))-xmin)/dx)+1
         if(bin .le. nbin)hist(bin)=hist(bin)+1.0
       enddo
    enddo
  enddo

  hist(:) = hist(:)/nb
  shist(:) = 0.0

  do i = 1, nbin
    xbin(i) = xmin+(i-0.5)*dx
    k = 0
    do j = max(1,i-half_window), min(nbin,i+half_window)
      shist(i) = shist(i) + hist(j)
      k = k +1
    enddo
    if(k > 0) shist(i) = shist(i)/real(k)
  enddo

  write(nrite_samp,102)'# X0sample:min_max :',xmin,xmax
  write(nrite_samp,101) xbin(1)-0.5*dx,zero,zero
  do i = 1, nbin
    write(nrite_samp,101) xbin(i),xbin(i)**2*hist(i),xbin(i)**2*shist(i),ip,ib
  enddo
  write(nrite_samp,101) xbin(nbin)+0.5*dx,zero,zero
  close(nrite_samp)
  endif

  endif

 101 format(f15.4,2(1x,f12.5),2(1x,i6))
 102 format(A,2(1x,f15.6),2(1x,i6))

  end subroutine printInitSamp


  subroutine decoh_damp(adia,istate,ener)
!**********************************************************************
!     SHARP Pack routine to implement decoherence of exopnential dampling
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
  implicit none

  integer,intent(in)       :: istate
  real*8,intent(in)        :: ener(nstates)
  complex*16,intent(inout) :: adia(nstates)

  integer  :: j
  real*8   :: tau_ji
  real*8   :: C,E_0,T_j
  complex*16   :: norm

  C = 1
  E_0 = 0.d0
  T_j = 1.d0
  norm = 0.d0

  do j = 1, nstates
    if(istate .ne. j)then
      tau_ji = hbar/(abs(ener(j) - ener(istate))) * (C + E_0/T_j)
      adia(j) = adia(j)*dexp(-dt/tau_ji)

!      sum c_j^2 for renormalization      
      norm = norm + CONJG(adia(j))*adia(j)
    endif
  enddo

  if(abs(adia(istate)) .gt. 1.0d-12)then
    norm = sqrt((1.d0-norm)/(CONJG(adia(istate))*adia(istate)))
    adia(istate) = adia(istate) * norm
  endif
  
  end subroutine decoh_damp

  subroutine sample_pimd(itraj,rp_pimd,vp_pimd)
!**********************************************************************
!     SHARP Pack routine to run pimd sampling
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
  implicit none


  integer,intent(in) :: itraj
  real*8,intent(out) :: rp_pimd(np,nb,ntraj),vp_pimd(np,nb,ntraj)

  integer :: i,j,ismpl
  real*8  :: systmp,rgr

  if((nsample .gt. 0).and.(nsample-nequil) .lt. ntraj)then
    write(0,*) 'Increase number of pimd nsample step'
    stop
  endif

  call sample_init(vp_samp,rp_samp)

!  transform postions normal to coordinates system
  if(lnorm .and. R0key.eq.1) call realft(rp_samp,np,nb,-1)
 
  open(nrite_pimd,file='pimd_sample.out',status='unknown',access='append')
  write(nrite_pimd,'(1x,A)') '# Traj  Time(au)   Temp   Rgr   rp(np,nb)     vp(np,nb) '

  ismpl = 0
  !nskip = int((nsample-nequil)/ntraj)

  do i=1,nsample

    call run_traj(rp_samp,vp_samp)

!   store sampling pos,vel to compute initial distribution
    if((i .gt. nequil).and.(mod(i,nskip).eq.0).and.(ismpl.lt.ntraj))then
      ismpl = ismpl + 1
      rp_pimd(:,:,ismpl) = rp_samp
      vp_pimd(:,:,ismpl) = vp_samp
    endif

!  compute radius of gyration (debugging mode)
    if(mod(i,iskip).eq.0)then
      call calcRg(itraj,i*dt,rp_samp,rgr)
    
! write pimd sampling trajectory for 1st-trajectory only        
      systmp = sum(vp_samp*vp_samp*mp)/(np*nb)*temp

      write(nrite_pimd,'(i7,1x,f15.3,1x,1000f15.5)') itraj,i*dt,systmp,rgr,(rp_samp(j,:),j=1,min(5,np)),(vp_samp(j,:),j=1,min(5,np))
    endif

  enddo
  close(nrite_pimd)

  if((nsample.gt.0).and.(ntraj.ne.ismpl))write(0,*) 'pimd sample NOT matching with ntraj!!'

  end subroutine sample_pimd
!**********************************************************************
  end module runtraj_module
