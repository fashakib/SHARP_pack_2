      module propagation_module
!**********************************************************************
!     SHARP PACK module that contains propragation schemes     
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
      use global_module
      use models_module
      use rpmd_module
      use initial_module
      use dboc_module
!      use nhc_module

      contains


      SUBROUTINE ADVANCE_ADIA(eva,vdotd_old,vdotd_new,adia)
!**********************************************************************
!     SHARP PACK subroutine to advance the electronic coefficients
!     by 4th order Runge-Kutta (RK) method
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none

      integer    :: i
      real*8,intent(in) :: eva(NSTATES,2)
      real*8,intent(in) :: vdotd_old(NSTATES,NSTATES),vdotd_new(NSTATES,NSTATES)
      real*8     :: eva_half(NSTATES)
      real*8     :: vdotd_half(NSTATES,NSTATES)
      complex*16 :: adia(NSTATES)
      complex*16 :: k1(NSTATES),k2(NSTATES),k3(NSTATES),k4(NSTATES)

!      real*8 :: avg
      !testing with eva_avg removing

!      avg = sum(eva(:,1))/nstates
!      eva(:,1) = eva(:,1)-avg
!      avg = sum(eva(:,2))/nstates
!      eva(:,2) = eva(:,2)-avg

      eva_half = 0.5d0*(eva(:,1)+eva(:,2))

      vdotd_half = 0.5d0*(vdotd_old+vdotd_new)
      
! 4th order RK scheme

      do i=1,nstates
        k1(i) = adia(i)*eva(i,1)/eye - sum(vdotd_old(i,:)*adia(:))
      end do

      do i=1,nstates
        k2(i) = (adia(i)+k1(i)*0.5d0*dtq)*eva_half(i)/eye- &
                sum(vdotd_half(i,:)*(adia(:)+k1(:)*0.5d0*dtq))
      end do

      do i=1,nstates
        k3(i) = (adia(i)+k2(i)*0.5d0*dtq)*eva_half(i)/eye- &
                sum(vdotd_half(i,:)*(adia(:)+k2(:)*0.5d0*dtq))
      end do

      do i=1,nstates
        k4(i) = (adia(i)+k3(i)*dtq)*eva(i,2)/eye- &
                sum(vdotd_new(i,:)*(adia(:)+k3(:)*dtq))
      end do

      adia = adia+dtq/6.d0*(k1+2.d0*k2+2.d0*k3+k4)

      return
      END SUBROUTINE ADVANCE_ADIA


      SUBROUTINE ADVANCE_MD(istate,rp,vp,fp,eva_bead,dhel_bead)
!**********************************************************************
!     SHARP PACK subroutine to advance the classical positions
!     and velocities by velocity verlet
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
!      use initial_module, only : gaussn
      use modelvar_module, only : psi_bead
      implicit none

      integer              :: i,ip,ibd
      integer              :: istart = 2
      integer, intent(in)  :: istate
      real*8, intent(inout):: rp(np,nb),vp(np,nb),fp(np,nb)
      real*8               :: fran(np,nb)
      real*8               :: frc_dboc(np)
      real*8               :: hel(NSTATES,NSTATES)
      real*8               :: p(np,nb)

      real*8               :: aa(np,nb)
      real*8               :: dt2
      real*8               :: rc(np)

      
      real*8   :: eva_bead(nstates,2,nb)
      real*8   :: dhel_bead(nstates,nstates,np,nb)

! for LinearChainModel
      real*8    :: fR(nb)

      frc_dboc = 0.d0

! RP propogation, velocity-Verlet
      dt2 = 0.5d0*dt

!     apply PILE thermostat - 1/2 step for t-rpmd dynamics
      if(ltrpmd)call pile_thermo(istart,vp)      

      !!Langevin dynamic for just N-th particle in LinearChain Model
      if((keymodel==12).or.(keymodel==13))then

        call Langevin(vp)

        aa = fp/mp

      else

        aa = fp/mp

      endif

     !!! Velocity-Verlet
      vp=vp+dt2*aa     !! fp/mp

      p=vp*mp

      CALL freerp(np,p,rp)
      !rp = rp + vp*dt

      vp=p/mp

      do ibd=1,nb
        if((keymodel==12).or.(keymodel==13))then
        ! LinearChainModel classical force calculation
          CALL FORCE_Lchain(rp(1:np,ibd), fp(1:np,ibd))

        else
         CALL gethel(rp(:,ibd),hel(:,:),dhel_bead(:,:,:,ibd))
         CALL Diag(eva_bead(:,2,ibd),psi_bead(:,:,2,ibd),hel(:,:))

         do i = 1, nstates
           psi_bead(:,i,2,ibd) = psi_bead(:,i,2,ibd) * &
                    dot_product(psi_bead(:,i,1,ibd),psi_bead(:,i,2,ibd))/ &
                    abs(dot_product(psi_bead(:,i,1,ibd),psi_bead(:,i,2,ibd)))
         end do

         CALL FORCE(psi_bead(:,:,2,ibd),istate,dhel_bead(:,:,:,ibd),fp(:,ibd))

!         psi_bead(:,:,1,ibd) = psi_bead(:,:,2,ibd)

!        force from dboc (2nd derivative nac) - bead approximation force         
         if(dvkey .eq. 2)then
           if(ldboc)then
             call compute_grad_dboc(np,ibd,rp(:,ibd),mp(:,ibd),istate,frc_dboc)

             fp(:,ibd) =  fp(:,ibd) + frc_dboc(:)
           endif
         endif

        endif
      enddo

!     force from dboc (2nd derivative nac) - centroid approximation force
      if(dvkey .eq. 1)then
        if(ldboc)then
          rc=0.d0
          do ibd=1,nb
            rc(:)=rc(:)+rp(:,ibd)
          enddo
          rc(:)=rc(:)/real(nb)

          call compute_grad_dboc(np,1,rc(:),mp(:,1),istate,frc_dboc)
          do ibd = 1, nb
            fp(:,ibd) =  fp(:,ibd) + frc_dboc(:)
          enddo

        endif
      endif

! second half of RP propogation

      if((keymodel==12).or.(keymodel==13))then
        aa = fp/mp

        vp=vp+dt2*aa    !!fp/mp

        call Langevin(vp)

      else
        aa = fp/mp
        vp=vp+dt2*aa    !!fp/mp
      endif

!      vp=vp+dt2*aa    !!fp/mp

!     apply PILE thermostat - 1/2 step for t-rpmd dynamics
      if(ltrpmd)call pile_thermo(istart,vp)      

      return
      END SUBROUTINE ADVANCE_MD
      
      SUBROUTINE FORCE(psi,istate,dhel,fp)
!**********************************************************************
!     SHARP PACK subroutine to calculate force 
!     by Hellmann–Feynman theorem on active surface
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************

      use global_module, only : np, NSTATES
      implicit none

      integer             :: i,j,ip
      integer, intent(in) :: istate
      real*8, intent(in)  :: psi(NSTATES,NSTATES)
      real*8, intent(in)  :: dhel(NSTATES,NSTATES,NP)
      real*8, intent(out) :: fp(np)

      fp=0.d0
      
      do ip=1,np
        do i=1,nstates   
          do j=1,nstates
            fp(ip)=fp(ip)-psi(i,istate)*dhel(i,j,ip)*psi(j,istate)
          enddo
        end do
      end do

      END SUBROUTINE FORCE


      SUBROUTINE FORCE_Lchain(r,ff)
!**********************************************************************
!     SHARP PACK subroutine to calculate force of LinearChainModel
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************

      use global_module, only : np
      implicit none

      integer             :: i,j,iL,iR
      real*8, intent(in)  :: r(np)
      real*8, intent(out) :: ff(np)
      real*8              :: rij,rr
      real*8              :: Vo,a

      Vo = 175.d0 * kJ_mol2au  !!0.00038125d0  !! kJ/mol --> a.u.
      a = 4.d0 * 0.5292d0 !! A^(-1) --> a.u^(-1)

      ff=0.d0

      do i=1,np
        iL = i-1
        iR = i+1
        if(i==1)iL=1
       ! if(i==np)iR=np

        do j=iL,iR
          if(i .ne. j)then
            if((i==np).and.(j.gt.np))then
              rij = r(i) - 2.0d0 !N+1 particle is fixed, !!10.d0
              rr = sqrt(rij*rij)
            else
              rij = r(i)-r(j)
              rr = sqrt(rij*rij)
            endif

            ff(i) = ff(i)-Vo*(2.d0*a*a*rr-3.d0*a*a*a*rr*rr+ &
                    2.32d0*a*a*a*a*rr*rr*rr)*rij/rr
          endif
        enddo

      enddo

      END SUBROUTINE FORCE_Lchain


      subroutine Langevin(vp)
!**********************************************************************
!     SHARP PACK subroutine to implement Langevin thermostat as 
!     implemented in Ceriotti2010, to LinearChain model in normal mode
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      use global_module, only : np,nb,dt,beta,dt,gamaLC ! ,pi,tau0
      use modelvar_module, only : mp

      implicit none

      integer :: ip, ib

      real*8,intent(inout)  :: vp(np,nb)
      real*8   :: C1(nb), C2(nb), ranf(nb)
      real*8  :: dt2

!      ph = pi/nb
!      betan = beta/nb
      dt2 = dt/2.d0

      do ib = 1, nb
!        wk(ib) = 2.d0/betan * sin((ib-1)*ph)
!        gama(ib) = 2.d0 * wk(ib)
!        if(ib==1) gama(ib) = 1.0/tau0
        C1(ib) = exp(-dt2*gamaLC(ib))
        C2(ib) = dsqrt(1.d0 - C1(ib)*C1(ib))
      enddo

      vp = vp*mp
!  Langevin start here
      call realft(vp,np,nb,+1)

      do ip = np, np
         do ib = 1, nb
            ranf(ib) = gaussn()
         enddo

         vp(ip,:) = C1*vp(ip,:) + dsqrt(mp(ip,:)*nb/beta)*C2 * ranf
      enddo

      call realft(vp,np,nb,-1)

      vp = vp/mp

      end subroutine Langevin


!**********************************************************************
      end module propagation_module
