      module rpmd_module
      implicit none

      contains

      subroutine run_traj(rp,vp)
!**********************************************************************
!     
!     SHARP PACK subroutine to evolve the free ring polymer according 
!     to its ground state potential for initial bead postion sampling

!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     
!**********************************************************************
      use global_module
      use modelvar_module, only : mp
      implicit none

      real*8             :: rp(np,nb),vp(np,nb)
      real*8             :: fp(np,nb),p(np,nb)
      real*8             :: dt2

      integer  :: istart = 1 

      dt2 = dt*0.5d0

!     apply PILE thermostat - 1/2 step
      if(lpimd)call pile_thermo(istart,vp)      

      call grhel(rp,fp)

      vp = vp + dt2*fp/mp

      p = vp*mp

! Evolve the ring polymer according to its ground state potential

      call freerp(np,p,rp)

      vp = p/mp

      call grhel(rp,fp)

      vp = vp + dt2*fp/mp

!     apply PILE thermostat - 1/2 step
      if(lpimd)call pile_thermo(istart,vp)      

      end subroutine run_traj


      subroutine grhel(rp,fp)
!**********************************************************************
!     SHARP PACK subroutine to calculate force based on ground state 
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      use global_module
      use modelvar_module, only : mp,Wb,cb
      implicit none

      integer             :: ip,ibd
      real*8              :: rp(np,nb)
      real*8, intent(out) :: fp(np,nb)
      real*8              :: hel(nstates,nstates)
      real*8              :: dhel(nstates,nstates,np)
      real*8              :: Wn

      do ip=1,np
         Wn = Wb(ip)
         do ibd=1,nb
           fp(ip,ibd) = -mp(ip,ibd)*Wn*Wn*(rp(ip,ibd)-R0)
         enddo
      enddo
   
      return

      end subroutine grhel


      subroutine freerp (nf,p,q)
!**********************************************************************
!     SHARP PACK routine for 
!     -----------------------------------------------------------------
!     Free harmonic ring-polymer evolution through a time interval dt.
!     -----------------------------------------------------------------
!
!**********************************************************************
      use global_module, only: nb
      use modelvar_module, only: mp

      implicit none
      integer,parameter :: nbmax=128   !!1024
      integer      :: j,k,nf
      integer,save :: init = 0

      real*8       :: p(nf,nb),q(nf,nb)
      real*8       :: pjknew
      real*8,save  :: poly(4,nbmax)
!
      if (init .eq. 0) then
         if (nb .gt. nbmax) stop 'freerp 1'
         call ring(nf,poly)
         init = 1
      endif

      if (nb .eq. 1) then
         do j = 1,nf
            q(j,1) = q(j,1)+p(j,1)*poly(3,1)/mp(j,1)
         enddo
      else
         call realft(p,nf,nb,+1)
         call realft(q,nf,nb,+1)
         do k = 1,nb
            do j = 1,nf
               pjknew = p(j,k)*poly(1,k)+q(j,k)*poly(2,k)*mp(j,k)
               q(j,k) = p(j,k)*poly(3,k)/mp(j,k)+q(j,k)*poly(4,k)
               p(j,k) = pjknew
            enddo
         enddo
         call realft(p,nf,nb,-1)
         call realft(q,nf,nb,-1)
      endif

      return
      end subroutine

      subroutine ring(nf,poly)
!**********************************************************************
!     SHARP PACK routine to calculate 
!     -----------------------------------------------------------------
!     Monodromy matrix elements for free ring-polymer evolution.
!     -----------------------------------------------------------------
!     
!**********************************************************************
      use global_module, only: nb,dt,beta,hbar
      use modelvar_module, only : nfreq

      implicit none

      integer  :: k,nf
      real*8   :: poly(4,nb)
      real*8   :: betan,twown,pibyn
      real*8   :: wk,wt,wm,cwt,swt
!
      poly(1,1) = 1.d0
      poly(2,1) = 0.d0
      poly(3,1) = dt !/mp
      poly(4,1) = 1.d0

      if (nb .gt. 1) then
!         betan = beta/nb
!         twown = 2.d0/(betan*hbar)
!         pibyn = dacos(-1.d0)/nb

         do k = 1,nb/2
!            wk = twown*dsin(k*pibyn)
            wk = nfreq(k+1)
            wt = wk*dt
!            wm = wk*mp
            cwt = dcos(wt)
            swt = dsin(wt)
            poly(1,k+1) = cwt
            poly(2,k+1) = -wk*swt
            poly(3,k+1) = swt/wk
            poly(4,k+1) = cwt
         enddo

         do k = 1,(nb-1)/2
            poly(1,nb-k+1) = poly(1,k+1)
            poly(2,nb-k+1) = poly(2,k+1)
            poly(3,nb-k+1) = poly(3,k+1)
            poly(4,nb-k+1) = poly(4,k+1)
         enddo
      endif

      return
      end subroutine

      subroutine realft(data,m,n,mode)
!**********************************************************************
!     SHARP PACK routine to calculate 
!     -----------------------------------------------------------------
!     FFT of m real arrays (if mode = 1) or complex Hermitian
!     arrays in real storage (if mode = -1), using -lfftw3.
!     Works equally well with f77 and ifc.
!     -----------------------------------------------------------------
!     
!**********************************************************************
      use modelvar_module, only : cmat
      use global_module,   only : lfft

      implicit none
      integer,parameter    :: nmax=1024

      integer              :: m,n,mode
      integer              :: j,k,np
      integer*8            :: plana,planb

      real*8               :: scale
      real*8               :: copy(nmax)
      real*8,intent(inout) :: data(m,n)

      data np /0/
      save copy,scale,plana,planb,np
!
      if(lfft)then
        if (n .ne. np) then
          if (n .gt. nmax) stop 'realft 1'
          scale = dsqrt(1.d0/n)
          call dfftw_plan_r2r_1d(plana,n,copy,copy,0,64)
          call dfftw_plan_r2r_1d(planb,n,copy,copy,1,64)
          np = n
        endif
      endif

      do k = 1,m

         if(lfft)then
            do j = 1,n
               copy(j) = data(k,j)
            enddo
         else
            copy(1:n) = 0.d0
         endif

         if (mode .eq. 1) then
            if(lfft)then
               call dfftw_execute(plana)
            else
               do j=1,n
                  copy(j)=dot_product(cmat(1:n,j),data(k,1:n))    
               enddo
            endif

         else if (mode .eq. -1) then
            if(lfft)then
               call dfftw_execute(planb)
            else
               do j=1,n
                  copy(j)=dot_product(cmat(j,1:n),data(k,1:n))    
               enddo
            endif

         else
            stop 'realft 2'
         endif

         if(lfft)then
            do j = 1,n
               data(k,j) = scale*copy(j)
            enddo
         else
            do j = 1,n
               data(k,j) = copy(j)
            enddo
         endif

      enddo

      return
      end subroutine realft


      subroutine cmat_init()
!**********************************************************************
!     SHARP PACK routine to initialize transform matrix 
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      use global_module, only : pi,nb,hbar,beta,pi,mdkey
      use modelvar_module, only : cmat,nfreq

      implicit none
      
      integer  :: j,k
      real*8   :: freq

!     normal mode frequencies
      freq = real(nb)/(beta*hbar)
      !freq = sqrt(real(nb))/(beta*hbar)
      do k=1,nb
        nfreq(k) = 2.d0*freq*dsin((k-1)*pi/nb)
      enddo

!     transformation matrix
      do j=1,nb

        cmat(j,1) = 1.d0/sqrt(real(nb))
      
        do k=1,nb/2

          cmat(j,k+1) = sqrt(2.d0/nb)*cos(2.d0*pi*j*k/nb)

        enddo

        do k=nb/2+1,nb-1

          cmat(j,k+1) = sqrt(2.d0/nb)*sin(2.d0*pi*j*k/nb)

        enddo

        if(mod(nb,2).eq.0) then
                
          cmat(j,nb/2+1)=1.d0/sqrt(real(nb))*(-1.d0)**j

        endif

      enddo

      end subroutine cmat_init

!**********************************************************************

      subroutine pile_thermo(istart,vp)
!**********************************************************************
!     SHARP PACK subroutine to implement PILE thermostat as in
!     PIMD/TRPMD in normal mode
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      use global_module, only : np,nb,dt,beta,taut,pi
      use modelvar_module, only : mp,pileC1,pileC2
      use initial_module, only : gaussn

      implicit none

      integer :: ip, ib, istart

      real*8, intent(inout)  :: vp(np,nb)
      real*8       :: ranf(nb)
      real*8       :: gama

      logical, save :: newjob = .true.

      if(newjob)then

        do ib = 1, nb
          gama = 4.d0*nb/beta * sin((ib-1)*pi/nb)   ! 2*wk

          if((nb.gt.1).and.(ib.eq.1))gama = 1.d0/taut

          pileC1(ib) = dexp(-0.5d0*dt*gama)
          pileC2(ib) = dsqrt(1.d0 - pileC1(ib)*pileC1(ib))

        enddo

!        write(*,*) newjob,pileC1,pileC2
!        pileC1(1) = 1, pileC2(1) = 0 for t-rpmd

        newjob = .false.

      endif

!     velocity to momenta      
      vp = vp*mp

!     normal mode tranformation 
      call realft(vp,np,nb,+1)

!     apply PILE thermostat to bead momenta - 1/2 step      
      do ip = 1, np
         do ib = 1, nb
            ranf(ib) = gaussn()
         enddo

!     istart = 1 (for pimd) and = 2 (for t-rpmd)
         do ib = istart,nb
           vp(ip,ib) = pileC1(ib)*vp(ip,ib) + &
                       dsqrt(mp(ip,ib)*nb/beta)*pileC2(ib)*ranf(ib)
         enddo
      enddo

!     back to 
      call realft(vp,np,nb,-1)

!      momenta to velocity
      vp = vp/mp

      end subroutine pile_thermo

!**********************************************************************
      end module rpmd_module
