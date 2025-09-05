      module nonadiabatic_module
!**********************************************************************
!     SHARP PACK module that contains nonadiabatic      
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************

      implicit none

      contains

      SUBROUTINE HOPPING(istate,inext,vdotd_new,adia)
!**********************************************************************
!     SHARP PACK subroutine to calculate hopping of PES      
!     IMPLEMENTING TULLY'S FEWEST SWITCH ALGORITHM.  
!     WHEN DISCONTINUITY IN THE WAVE FUNCTION IS DETECTED, 
!     DIABATIC MOTION IS ASSUMED
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************

      use global_module, only : NSTATES,dtq
      implicit none

!#------------- dummy arguments ------------
      integer, intent(in)    :: istate
      integer, intent(inout) :: inext
      real*8, intent(in)     :: vdotd_new(NSTATES,NSTATES)
      complex*16, intent(in) :: adia(NSTATES)

!#------------- local variables ------------
      integer    :: i,j,k
      real*8     :: acumulator,dice,rholl,gik(NSTATES)
      real*8     :: RANF 
      complex*16 :: aux1(NSTATES)  !#  TIME DERIVATIVES OF PSI I
      complex*16 :: aux2(NSTATES)  !#  TIME DERIVATIVES OF PSI II
      complex*16 :: rho(NSTATES)

!     numerically integrate gik
      rholl = CONJG(adia(istate))*adia(istate)
      CALL RANDOM_NUMBER(ranf)
      dice = RANF
      acumulator = 0.0
      rho(:) = (adia(:)*CONJG(adia(istate)))

      aux2 = (0.0, 0.0)
      DO j = 1, NSTATES
         aux2(j) = -2.0 * dtq * vdotd_new(j,istate)
      END DO

       gik(:) = REAL(rho(:)*aux2(:), 8) / rholl

! downward transitions (if any)

      DO i=1,istate-1             !replaced inext with istate
           IF (gik(i).GT.0) THEN
            acumulator=acumulator+gik(i)
             IF(acumulator.GT.dice) THEN
               inext=i         ! we jumped!! (still have to check energy)
               GOTO 10
             END IF
          END IF
      END DO

! upward transitions (if any)

      DO i=istate+1,NSTATES
         IF (gik(i).GT.0)THEN
           acumulator=acumulator+gik(i)
            IF(acumulator.GT.dice)THEN
                inext=i         ! we jumped!! (still have to check energy)
               GOTO 10
            END IF
         END IF
      END DO

  10 CONTINUE
             
      END SUBROUTINE HOPPING


      SUBROUTINE EKINRESCALE(istate,inext,vc,vp,rp,eva,eva_b,dhel_rc,psi,ctime)
!**********************************************************************
!     SHARP PACK subroutine to scale velocity (KE) after successful
!     hopping in between two adiabatic energy states
!
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      use global_module
      use modelvar_module, only : nJump, nJumpFail, nFrust, nFrustR,mp
      use models_module
      use propagation_module 
      implicit none

      integer                :: i,j,ip,ibd
      integer, intent(inout) :: istate,inext  
      real*8, intent(inout)  :: vp(np,nb)
      real*8, intent(in)     :: eva(NSTATES,2),vc(np)
      real*8, intent(in)     :: eva_b(NSTATES,nb)
      real*8, intent(in)     :: dhel_rc(NSTATES,NSTATES,np)
      real*8, intent(in)     :: psi(NSTATES,NSTATES,2)
      real*8, intent(in)     :: rp(np,nb)
       ! derivative coupling from k to k' state
      real*8                 :: dcoup(np),a_tot,b_tot
      real*8                 :: ediff,gama,det
      real*8                 :: fp1(np),fp2(np)
      real*8                 :: f1,f2
      real*8                 :: avetim
      integer                :: npscale
      integer                :: ctime

      avetim = 0.d0

      if(vckey == 1)then 
        ediff=eva(istate,2)-eva(inext,2)

      elseif(vckey == 2)then
      !!bead_level
        ediff = (sum(eva_b(istate,:))-sum(eva_b(inext,:)))/nb

      endif
!!!
      dcoup=0.0

      a_tot=0.d0
      b_tot=0.d0

      if((keymodel==12).or.(keymodel==13))then

        avetim = 30.d0

        do ip=1,1   !np !!! VELOCITY SCALING JUST P-1 !!!  
          
          dcoup(ip) = d_ab(istate,inext) 

          a_tot = a_tot + (0.5d0*dcoup(ip)*dcoup(ip))/mp(ip,1)
          b_tot = b_tot + dcoup(ip)*vc(ip)     
        enddo

      else

        do i=1,nstates
          do j=1,nstates
            do ip=1,np
              dcoup(ip) = dcoup(ip)+psi(i,istate,2)*dhel_rc(i,j,ip)* &
                          psi(j,inext,2)/(eva(inext,2)-eva(istate,2))
            enddo
          enddo
        enddo
 
        do ip=1,np 
          a_tot = a_tot + (0.5d0*dcoup(ip)*dcoup(ip))/mp(ip,1)
          b_tot = b_tot + dcoup(ip)*vc(ip)     
        enddo

      endif
    
      if(ldtl)then
        WRITE(nrite_hopp,*)
        WRITE(nrite_hopp,'(" Trying to hop from ",I3," to ",I3)') &
                           istate,inext
      endif

      det = b_tot*b_tot + 4.d0*a_tot*ediff
      
      IF (det>=0) THEN
        if(ldtl)WRITE(nrite_hopp,'(" Enough Kinetic energy to hop, &
                                   accepted")')
        if((ctime*dt/41340.d0) > avetim)then
           nJump(istate,inext) = nJump(istate,inext) + 1
        endif

         IF (b_tot<0) THEN
            gama=(b_tot+dsqrt(det))/(2.d0*a_tot)
         ELSE
            gama=(b_tot-dsqrt(det))/(2.d0*a_tot)
         ENDIF
         istate=inext

      !! FRUNSTRATED HOPPING
      ELSEIF (det<0) THEN
        if((ctime*dt/41340.d0)>avetim)then
        nfrust_hop = nfrust_hop + 1
        nFrust(istate,inext) = nFrust(istate,inext) + 1
        endif

        ! Never reverse velocity
        if(vrkey == 0)then
          gama = 0.d0
          if((ctime*dt/41340.d0)>avetim)then
          nfrust_hop2 = 0
          nFrustR(istate,inext) = nFrustR(istate,inext) + 0
          endif
        ! Always reverse velocity based on Hammes-Tully1994
        elseif(vrkey == 1)then
          gama = b_tot/a_tot
          if((ctime*dt/41340.d0)>avetim)then
          nfrust_hop2 = nfrust_hop2 + 1
          nFrustR(istate,inext) = nFrustR(istate,inext) + 1
          endif

        ! Truhlar scheme of velocity reversal, Chem. Phys. Lett. 369, 60 (2003)
        elseif(vrkey == 2)then
          gama = 0.d0

          CALL FORCE(psi(:,:,2),inext,dhel_rc(:,:,1:np), fp2(1:np))
         
          f2 = sum(dcoup*fp2)

          !! Truhlar condtion as in J.Chem.Phys.147,214112(2017) 
          if((b_tot*f2) < 0.d0)then
            gama = b_tot/a_tot
            if((ctime*dt/41340.d0)>avetim)then
            nfrust_hop2 = nfrust_hop2 + 1
            nFrustR(istate,inext) = nFrustR(istate,inext) + 1
            endif
          endif
        
        !! Jain & Subotnik condtion as 
        elseif(vrkey == 3)then
          gama = 0.d0

          CALL FORCE(psi(:,:,2),istate,dhel_rc(:,:,1:np), fp1(1:np))
          CALL FORCE(psi(:,:,2),inext,dhel_rc(:,:,1:np), fp2(1:np))
         
          f1 = sum(dcoup*fp1)
          f2 = sum(dcoup*fp2)

          !! Truhlar condtion as in J.Chem.Phys.147,214112(2017) 
          if(((b_tot*f2) < 0.d0).and.((f1*f2) < 0.d0))then
            gama = b_tot/a_tot
            if((ctime*dt/41340.d0)>avetim)then
            nfrust_hop2 = nfrust_hop2 + 1
            nFrustR(istate,inext) = nFrustR(istate,inext) + 1
            endif
          endif

        endif

        if(ldtl)WRITE(nrite_hopp,'(" Not enough kinetic energy, &
                                    rejected")')
        if((ctime*dt/41340.d0)>avetim)then
          nJumpFail(istate,inext) = nJumpFail(istate,inext) + 1
        endif
        inext=istate
      ENDIF

!     RESCALING VELOCITY ALONG DIRECTION OF dij TO CONSERVE ENERGY     
      !! velocity rescaling for P-1 only
      if((keymodel==12).or.(keymodel==13))then
        DO ip=1,1
          DO ibd=1,nb
            vp(ip,ibd)=vp(ip,ibd) - gama*dcoup(ip)/mp(ip,ibd)
          ENDDO
        ENDDO 

      else
        DO ip=1,np
          DO ibd=1,nb
            vp(ip,ibd)=vp(ip,ibd) - gama*dcoup(ip)/mp(ip,ibd)
          ENDDO
        ENDDO 

      endif

      if(ldtl)WRITE(nrite_hopp,*)
      
      END SUBROUTINE EKINRESCALE

      
      SUBROUTINE pop_estimator(istate,itime,diabat1,diabat2,diabat3, &
                 adiabat1,adiabat2)
!**********************************************************************
!     SHARP PACK subroutine to calculate diabatic/adiabatic populations
!
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
      use global_module, only : nstates, nprint, nb
      use modelvar_module, only : psi, adia
      implicit none

      integer                  :: i,j,k,l
      integer, intent(in)      :: istate,itime
      real*8, intent(inout)    :: diabat1(NSTATES,0:nprint)
      real*8, intent(inout)    :: diabat3(NSTATES,0:nprint)
      real*8, intent(inout)    :: adiabat1(NSTATES,0:nprint)
      complex*16, intent(inout):: diabat2(NSTATES,0:nprint)
      complex*16, intent(inout):: adiabat2(NSTATES,NSTATES,0:nprint)
      complex*16               :: diabat(nstates)

!=============store the adiabatic density matrix by differnet ways============
      ! bined coefficient  
      adiabat1(istate,itime)=adiabat1(istate,itime)+1.

      ! electronic coeffcienct
      do i=1,NSTATES
        do j=1,NSTATES
          adiabat2(i,j,itime)=adiabat2(i,j,itime)+adia(i)*CONJG(adia(j))
        end do
      end do

!=============store the diabatic density matrix by differnet ways============
! Method 1 from Subotink, JCP, 139, 211101, 2013

      do i=1,NSTATES
        diabat1(i,itime)=diabat1(i,itime)+psi(i,istate,1)**2
      enddo

! Method 2 from Subotink, JCP, 139, 211101, 2013

      diabat = cmplx(0.,0.)
      do i=1,NSTATES
        do k=1,NSTATES
          diabat(i)=diabat(i)+psi(i,k,1)*adia(k)
        enddo
      enddo

      do i=1,NSTATES
        diabat2(i,itime)=diabat2(i,itime)+abs(diabat(i))**2
      end do

! Method 3 from Subotink, JCP, 139, 211101, 2013

      do i=1,NSTATES
        diabat3(i,itime)=diabat3(i,itime)+psi(i,istate,1)**2
      enddo

      do k=1,NSTATES
        do l=1,NSTATES
          if (k.lt.l) then
            do i=1,NSTATES
              diabat3(i,itime)=diabat3(i,itime)+2*Real(psi(i,k,1)* &
                               adia(k)*conjg(adia(l))*psi(i,l,1))
            enddo
          endif
        enddo
      enddo

      END SUBROUTINE

      SUBROUTINE EKINRESCALE_BEAD(istate,inext,vp,rp,eva,eva_b,delH,psi,ctime)
!**********************************************************************
!     SHARP PACK subroutine to scale velocity (KE) after successful
!     hopping in between two adiabatic energy states
!     using bead approximation
!
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      use global_module
      use modelvar_module, only : nJump, nJumpFail, nFrust, nFrustR,mp
      use models_module
      use propagation_module 
      implicit none

      integer                :: i,j,ip,ib
      integer, intent(inout) :: istate,inext  
      real*8, intent(inout)  :: vp(np,nb)
      real*8, intent(in)     :: rp(np,nb)
      real*8, intent(in)     :: eva(NSTATES,2)
      real*8, intent(in)     :: eva_b(NSTATES,nb)
      real*8, intent(in)     :: delH(NSTATES,NSTATES,np,nb)
      real*8, intent(in)     :: psi(NSTATES,NSTATES,nb)
       ! derivative coupling from k to k' state
      real*8                 :: dcoup(np,nb),a_tot,b_tot
      real*8                 :: ediff,gama,det
      real*8                 :: fp1(np),fp2(np)
      real*8                 :: f1,f2
      real*8                 :: avetim
      integer                :: npscale
      integer                :: ctime

      avetim = 0.d0
      if(vckey == 1)then 
        ediff=eva(istate,2)-eva(inext,2)

      elseif(vckey == 2)then
      !!bead_level
        ediff = (sum(eva_b(istate,:))-sum(eva_b(inext,:)))/nb

      endif

      dcoup=0.0
      a_tot=0.d0
      b_tot=0.d0

      if((keymodel==12).or.(keymodel==13))then

        avetim = 30.d0

        do ip=1,1   !np !!! VELOCITY SCALING JUST P-1 !!! 
          do ib=1,nb 
          
            dcoup(ip,ib) = d_ab(istate,inext) 

            a_tot = a_tot + (0.5d0*dcoup(ip,ib)*dcoup(ip,ib))/mp(ip,ib)
            b_tot = b_tot + dcoup(ip,ib)*vp(ip,ib)     
          enddo
        enddo

      else

        do i=1,nstates
          do j=1,nstates
            do ip=1,np
              do ib=1,nb
                dcoup(ip,ib) = dcoup(ip,ib)+psi(i,istate,ib) * &
                     delH(i,j,ip,ib)* psi(j,inext,ib) / &
                     (eva_b(inext,ib)-eva_b(istate,ib))
              enddo
            enddo
          enddo
        enddo
 
        do ip=1,np 
          do ib=1,nb
            a_tot = a_tot + (0.5d0*dcoup(ip,ib)*dcoup(ip,ib))/mp(ip,ib)
            b_tot = b_tot + vp(ip,ib)*dcoup(ip,ib)     
          enddo
        enddo

      endif

      a_tot = a_tot/real(nb)
      b_tot = b_tot/real(nb)

    
      if(ldtl)then
        WRITE(nrite_hopp,*)
        WRITE(nrite_hopp,'(" Trying to hop from ",I3," to ",I3)') &
                           istate,inext
      endif

      det = b_tot*b_tot + 4.d0*a_tot*ediff
      
      IF (det>=0) THEN
        if(ldtl)WRITE(nrite_hopp,'(" Enough Kinetic energy to hop, &
                                   accepted")')
        if((ctime*dt/41340.d0) > avetim)then
           nJump(istate,inext) = nJump(istate,inext) + 1
        endif

         IF (b_tot<0) THEN
            gama=(b_tot+dsqrt(det))/(2.d0*a_tot)
         ELSE
            gama=(b_tot-dsqrt(det))/(2.d0*a_tot)
         ENDIF
         istate=inext

      !! FRUNSTRATED HOPPING
      ELSEIF (det<0) THEN
        if((ctime*dt/41340.d0)>avetim)then
        nfrust_hop = nfrust_hop + 1
        nFrust(istate,inext) = nFrust(istate,inext) + 1
        endif

        ! Never reverse velocity
        if(vrkey == 0)then
          gama = 0.d0
          if((ctime*dt/41340.d0)>avetim)then
          nfrust_hop2 = 0
          nFrustR(istate,inext) = nFrustR(istate,inext) + 0
          endif
        ! Always reverse velocity based on Hammes-Tully1994
        elseif(vrkey == 1)then
          gama = b_tot/a_tot
          if((ctime*dt/41340.d0)>avetim)then
          nfrust_hop2 = nfrust_hop2 + 1
          nFrustR(istate,inext) = nFrustR(istate,inext) + 1
          endif

        ! Truhlar scheme of velocity reversal, Chem. Phys. Lett. 369, 60 (2003)
        elseif(vrkey == 2)then
          gama = 0.d0

!          CALL FORCE(psi(:,:,2),inext,dhel_rc(:,:,1:np), fp2(1:np))
         
!          f2 = sum(dcoup*fp2)

          !! Truhlar condtion as in J.Chem.Phys.147,214112(2017) 
!          if((b_tot*f2) < 0.d0)then
!            gama = b_tot/a_tot
!            if((ctime*dt/41340.d0)>avetim)then
!            nfrust_hop2 = nfrust_hop2 + 1
!            nFrustR(istate,inext) = nFrustR(istate,inext) + 1
!            endif
!          endif
        
        !! Jain & Subotnik condtion as 
        elseif(vrkey == 3)then
          gama = 0.d0

!          CALL FORCE(psi(:,:,2),istate,dhel_rc(:,:,1:np), fp1(1:np))
!          CALL FORCE(psi(:,:,2),inext,dhel_rc(:,:,1:np), fp2(1:np))
         
!          f1 = sum(dcoup*fp1)
!          f2 = sum(dcoup*fp2)

          !! Truhlar condtion as in J.Chem.Phys.147,214112(2017) 
!          if(((b_tot*f2) < 0.d0).and.((f1*f2) < 0.d0))then
            gama = b_tot/a_tot
!            if((ctime*dt/41340.d0)>avetim)then
!            nfrust_hop2 = nfrust_hop2 + 1
!            nFrustR(istate,inext) = nFrustR(istate,inext) + 1
!            endif
!          endif

        endif

        if(ldtl)WRITE(nrite_hopp,'(" Not enough kinetic energy, &
                                    rejected")')
        if((ctime*dt/41340.d0)>avetim)then
          nJumpFail(istate,inext) = nJumpFail(istate,inext) + 1
        endif
        inext=istate
      ENDIF

!     RESCALING VELOCITY ALONG DIRECTION OF dij TO CONSERVE ENERGY     
      !! velocity rescaling for P-1 only
      if((keymodel==12).or.(keymodel==13))then
        DO ip=1,1
          DO ib=1,nb
            vp(ip,ib)=vp(ip,ib) - gama*dcoup(ip,ib)/mp(ip,ib)
          ENDDO
        ENDDO 

      else
        DO ip=1,np
          DO ib=1,nb
            vp(ip,ib)=vp(ip,ib) - gama*dcoup(ip,ib)/mp(ip,ib)
          ENDDO
        ENDDO 

      endif

      if(ldtl)WRITE(nrite_hopp,*)
      
      END SUBROUTINE EKINRESCALE_BEAD

!**********************************************************************
      end module nonadiabatic_module

