      module initial_module
      use global_module
      use modelvar_module, only : mp,Wb,rp0
      use modelvar_module, only : nfreq

      implicit none

      contains

      subroutine modelParam(keymodel)
!**********************************************************************
!     SHARP PACK routine to set some specific  model parameters
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none
      
      integer  :: ib,keymodel

!mass, beta parameters for Tully's Models
      if(keymodel .le. 3)then  
        mp = 2000.0d0
        beta = mp(1,1)/(P0*P0)   !mp/(P0**2)
!        beta = 10000

        R0 = -15.d0

! Parameters for Morse1 Model
      elseif(keymodel .eq. 4)then
        mp = 20000.0d0
        beta = 1052.58d0   ! beta corresponding to 300K
        omega = 0.005d0

        R0 = 2.9d0

! Parameters for Morse2 Model
      elseif(keymodel .eq. 5)then
        mp = 20000.0d0
        beta = 1052.58d0   ! beta corresponding to 300K
        omega = 0.005d0

        R0 = 3.3d0
         
! Parameters for Morse3 Model
      elseif(keymodel .eq. 6)then
        mp = 20000.0d0
        beta = 1052.58d0   ! beta corresponding to 300K
        omega = 0.005d0

        R0 = 2.1d0
         
! Mass, beta parameters for Spin-Boson Model
      elseif(keyModel .eq. 7)then  
        mp = 1.0d0
        if(keybath == 1 .and. enu .gt. 1.d0)then
          enu = enu/freq
        elseif(keybath == 1 .and. enu .eq. 1.0)then
          enu = 1.d0
        elseif(keybath == 2)then
          enu = 1.d0
        endif

        eps= eps * enu
        delta= delta * enu
        beta = 1.d0/(KT*enu)
        wc = wc * enu ! * (2.d0*pi)
        E_r = E_r * enu
!       wmax = wmax * wc
        R0 = 0

!Mass, beta parameters for 2-state linear chain Model
!mass, beta parameters for 2-state linear chain Model
      elseif(keymodel .eq. 12)then  
        mp = 12.d0 * amu2au   !! amu --> a.u.
        beta = 1.d0/kT*temp    !! a.u.
        v11 = 0.d0 * kJ_mol2au  !! kJ/mol --> a.u.
        v22 = 8.0d0 * kJ_mol2au  !! kJ/mol --> a.u.
        d_ab = 0.d0
        d_ab(1,2)= -6.d0 * 0.52918d0  !-6.0 A^(-1)
        d_ab(2,1)=-d_ab(1,2)
        gamaLC(1) = 0.002418d0   !! 10^14 S^-1 --> a.u.
        sigmaLC(1) = sqrt(2.d0*gamaLC(1)*mp(1,1)*nb/beta/dt) 
        do ib = 2, nb  !! ib = 2, nb
          gamaLC(ib) = 2.d0*(2.d0*nb/(beta*hbar)*sin((ib-1)*pi/nb)) 
          sigmaLC(ib) = sqrt(2.d0*gamaLC(ib)*mp(1,ib)*nb/beta/dt)
        enddo 

!mass, beta parameters for 3-state SuperExchange linear chain Model
      elseif(keymodel .eq. 13)then  
        mp = 12.d0 * amu2au   !! amu --> a.u.
        beta = 1.d0/kT*temp    !! a.u.
        v11 = 0.0   !!a.u.
        v22 = 39.d0*kJ_mol2au  !!kJ/mol --> a.u.
        v33 = 13.d0*kJ_mol2au  !!kJ/mol --> aa.u.
        d_ab = 0.d0
        d_ab(1,2)= -6.0d0 * 0.52918d0 !!a.u.
        d_ab(2,1)=-d_ab(1,2) 
        d_ab(2,3)= 8.d0 * 0.52918d0 !!a.u.
        d_ab(3,2)=-d_ab(2,3)
        gamaLC(1) = 0.002418d0   !! 10^14 S^-1 --> a.u.
        sigmaLC(1) = sqrt(2.d0*gamaLC(1)*mp(1,1)*nb/beta/dt) 
        do ib = 2, nb
          gamaLC(ib) = 2.d0*(2.d0*nb/(beta*hbar)*sin((ib-1)*pi/nb)) 
          sigmaLC(ib) = sqrt(2.d0*gamaLC(ib)*mp(1,ib)*nb/beta/dt)
        enddo

      elseif(keymodel .eq. 14)then  
        mp = 2000.0d0
        beta = mp(1,1)/(P0*P0)   !mp/(P0**2)

        R0 = -10.d0

!mass, beta parameters for Tully's Models
      elseif(keymodel .eq. 21)then  
        mp = 2000.0d0
        beta = mp(1,1)/(P0*P0)   !mp/(P0**2)

        R0 = -10.d0

      endif

      return
      end subroutine modelParam


      subroutine sample_init(vp,rp)
!**********************************************************************
!     SHARP PACK routine to sample initial velocity and position 
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none

      integer             :: ibd,ip
      real*8, intent(out) :: vp(np,nb),rp(np,nb)
      real*8              :: sigR,sigP,alpha,sigRw,sigPw
      real*8              :: R0c

      if(keymodel .le. 3)then
         alpha = 0.25d0
         sigR = dsqrt(1.d0/(2.d0*alpha))
         sigP = dsqrt(hbar*hbar*alpha/2.d0)
         sigRw = dsqrt(1.d0/(4.d0*alpha))
         sigPw = dsqrt(hbar*hbar*alpha)
      elseif(keymodel==4)then
         sigR = dsqrt(hbar/(2.d0*mp(1,1)*omega))
         sigP = dsqrt(mp(1,1)*hbar*omega/2.d0)
         sigRw = sigR
         sigPw = sigP
         V0key = 1
      elseif(keymodel==5)then
         sigR = dsqrt(hbar/(2.d0*mp(1,1)*omega))
         sigP = dsqrt(mp(1,1)*hbar*omega/2.d0)
         sigRw = sigR
         sigPw = sigP
         V0key = 1
      elseif(keymodel==6)then
         sigR = dsqrt(hbar/(2.d0*mp(1,1)*omega))
         sigP = dsqrt(mp(1,1)*hbar*omega/2.d0)
         sigRw = sigR
         sigPw = sigP
         V0key = 1
      elseif(keymodel==14)then
         alpha = 0.25d0
         sigR = dsqrt(1.d0/(2.d0*alpha))
         sigP = dsqrt(hbar*hbar*alpha/2.d0)
         sigRw = dsqrt(1.d0/(4.d0*alpha))
         sigPw = dsqrt(hbar*hbar*alpha)
      elseif(keymodel==21)then
         alpha = (P0/20.d0)*(P0/20.d0)   !0.25d0
         sigR = dsqrt(1.d0/(2.d0*alpha))
         sigP = dsqrt(hbar*hbar*alpha/2.d0)
         sigRw = dsqrt(1.d0/(4.d0*alpha))
         sigPw = dsqrt(hbar*hbar*alpha)
      endif

      do ip=1,np
         
         if(keymodel == 7)then

           if(v0key == 1)then
              ! initialize Gaussian distributed bead momenta 
              sigP = dsqrt(mp(ip,1)*real(nb)/beta)
           elseif(v0key == 2)then
              ! initialize Wigner distributed bead momenta 
              sigPw = dsqrt(mp(ip,1)*hbar*Wb(ip)/(2.d0*tanh(0.5d0*beta/nb*Wb(ip))))
           else
              write(0,*) 'INCORRECT MOMENTUM INIT SCHEME SELECTED!!!'
              stop
           endif

           if(R0key .lt. 2)then
              ! initialize Gaussian distributed bead position 
              sigR = dsqrt(mp(ip,1)*real(nb)/beta)/(mp(ip,1)*Wb(ip))
           elseif(R0key == 2)then
              ! initialize Wigner distributed bead position 
              sigRw = dsqrt(mp(ip,1)*hbar*Wb(ip)/(2.d0*tanh(0.5d0*beta/nb*Wb(ip))))/(mp(ip,1)*Wb(ip))
           endif

           if(R0key==0)sigR=gaussn()*sigR + R0
         endif
         
         do ibd=1,nb
            ! initialize bead positions
            if(R0key==0)then
              ! assign each bead same position with fraction of de-Broglie length
              rp(ip,ibd) = sigR + gaussn()*(0.001)*(hbar*sqrt(beta/(mp(ip,ibd)*nb)))
            elseif(R0key==1)then
              if(lnorm.and.(keymodel.eq.7))sigR = sigP/mp(ip,1)/dsqrt(Wb(ip)**2+nfreq(ibd)**2)
              ! assign each bead gaussian distributed position
              rp(ip,ibd) = gaussn()*sigR + R0 

            elseif(R0key==2)then
              ! assign each bead wigner distributed position
              rp(ip,ibd) = gaussn()*sigRw + R0  
            endif

            ! initialize bead momenta
            if(v0key == 0)then
               ! initialize deterministic bead momenta 
               vp(ip,ibd) = P0
            elseif(v0key == 1)then
               ! initialize Gaussian distributed bead momenta 
               vp(ip,ibd) = gaussn()*sigP + P0
            elseif(v0key == 2)then
               ! initialize Wigner distributed bead momenta 
               vp(ip,ibd) = gaussn()*sigPw + P0
            else
               write(0,*) 
               write(0,*) ' ERROR!! '
               write(0,*) ' INCORRECT MOMENTUM INIT SCHEME SELECTED!!!'
               write(0,*) 
               stop
            endif
         end do

      end do
      
      vp = vp/mp

      return
      end subroutine sample_init


      subroutine sample_init_lchain(vp,rp)
!**********************************************************************
!     SHARP PACK routine to sample initial velocity and position for
!     N-Linear Chain Model
!
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none

      integer             :: ibd,ip
      real*8, intent(out) :: vp(np,nb),rp(np,nb)
      real*8              :: sigR,sigP

      sigP = dsqrt(mp(1,1)*real(nb)/beta)
      
      do ip=1,np

         do ibd=1,nb
           ! initialize bead positions
           rp(ip,ibd) = (ip-1)*2.d0/np + 0.01d0*gaussn()
           !rp(ip,ibd) = gaussn()

           ! initialize Gaussian distributed bead momenta 
           vp(ip,ibd) = gaussn()*sigP
         enddo
      enddo

      vp = vp/mp

      end subroutine sample_init_lchain
      

      FUNCTION gaussn()
!**********************************************************************
!     SHARP PACK function to calculate normal distribution center at
!     origin with unit standard deviation 
!     
!**********************************************************************
      IMPLICIT NONE

      REAL(8), PARAMETER  :: PI2=2.0*3.141592654
      REAL(8)             :: Z1,Z2,gaussn

      CALL RANDOM_NUMBER(z1)
      CALL RANDOM_NUMBER(z2)
      gaussn=dSQRT(-2.d0*dLOG(Z1))*dCOS(PI2*Z2)
!      gaussn=0.00193
      RETURN
      END FUNCTION gaussn


!**********************************************************************
!     
!    bath initialization for spin-boson models 
! 
!**********************************************************************
      subroutine bath_init(keybath)
      use modelvar_module, only: Wb,cb
      implicit none
      
      integer            :: i
      integer,intent(in) :: keybath
      real*8             :: delw(np),Jw(np),Jwcal(np)
      real*8             :: wmax


      open(nrite_bath,file='bathfrequency.out',status='unknown')
      write(nrite_bath,'(A3,1x,5(A10,6x))')'#n ','w_i','c_i','delw_i','J(w)','Jcal(w)'
      write(nrite_bath,'(A,1x,f15.8)')'#enu: ',enu

      if(keybath == 1)then
!     Bath discretization: Debye Spectral form
        wmax = 20.d0*wc

        do i=1,np
          Wb(i) = tan((real(i)*atan(wmax/wc))/real(np))*wc
          cb(i) = Wb(i)*dsqrt(E_r*atan(wmax/wc)/(pi*real(np)))
        enddo

        delw(1) = wb(2)-wb(1)
        do i=2,np-1
          delw(i) = (wb(i+1)-wb(i-1))/2.d0
        enddo
        !two points linear interpolation
        delw(np) = delw(np-1)+(wb(np)-wb(np-1))* &
                (delw(np-1)-delw(np-2))/(wb(np-1)-wb(np-2))

        do i=1,np
          Jw(i) = 0.5d0*E_r*wb(i)*wc/(wb(i)*wb(i) + wc*wc)
          Jwcal(i) = pi/2.d0 * cb(i)*cb(i)/wb(i)/delw(i)
        enddo

      elseif(keybath == 2)then
!     Bath discretization: Ohmic Spectral form
!     J. Chem. Phys. 151, 024105 (2019)
!        wmax = wc/nb*(1-exp(-3.d0))
        do i=1,np
!          Wb(i) = -wc * log(1.0 - real(i)*wmax/wc)
!          cb(i) = Wb(i) * dsqrt(E_r*wmax)
          Wb(i) = -wc * log(1.0 - real(i)/(1+np))
          cb(i) = Wb(i) * dsqrt(E_r*wc/(1+np))
        enddo

        delw(1) = wb(2)-wb(1)
        do i=2,np-1
          delw(i) = (wb(i+1)-wb(i-1))/2.d0
        enddo
        !two points linear interpolation
        delw(np) = delw(np-1)+(wb(np)-wb(np-1))* &
                (delw(np-1)-delw(np-2))/(wb(np-1)-wb(np-2))

        do i=1,np
          Jw(i) = pi/2.d0 * E_r * wb(i) * dexp(-wb(i)/wc)
          Jwcal(i) = pi/2.d0 * cb(i)*cb(i)/wb(i)/delw(i)
        enddo

      endif

      do i=1,np
        write(nrite_bath,'(i4,5e15.6)') i, wb(i)/enu,cb(i)/(enu*sqrt(enu)),&
                delw(i)/enu, Jw(i)/enu, Jwcal(i)/enu
      enddo
      
      close(nrite_bath)

      end subroutine bath_init

!**********************************************************************

    end module initial_module
