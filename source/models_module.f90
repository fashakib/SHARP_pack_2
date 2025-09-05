      module models_module
!**********************************************************************
!     SHARP Pack module for defining system models
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
      use global_module

      implicit none

      real*8,private :: A, B, C, D, Eo
      real*8,private :: Dm(3),bm(3),Re(3), cm(3)
      real*8,private :: Aij(3), Rij(3), amj(3)
      real*8,private :: aa(3),bb(3)  ! for 3-state super exchange model

      contains

      subroutine gethel(rxyz,hel, dhel)
!**********************************************************************
!     
!     SHARP PACK subroutine to define Hamiltonian of models
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!
!           | V11   V12   ... |
!      H  = | V21   V22   ... |
!           | ...   ...   ... |
!  
!     keymodel = 1-3 (Tully Model 1,2,3)
!     keymodel = 4-6 (Morse Model 1,2,3)     
!     keymodel = 12 (2-state coupled with N-linear chain)     
!     keymodel = 13 (2-state coupled with N-linear chain)     
!     
!**********************************************************************
      use modelvar_module, only : Wb, cb

      implicit none

      real*8,intent(in) :: rxyz(np)
      real*8,intent(out) :: hel(NSTATES,NSTATES)
      real*8,intent(out) :: dhel(NSTATES,NSTATES,np)

      real*8,parameter :: m=1.d0
      real*8  :: Ho(nstates,nstates)
      real*8  :: Hba, Hsb

      integer :: i, j, ip

      real*8  :: ff(np),Vq

      hel = 0.d0
      dhel = 0.d0

! if(keymodel == 1) <= Tully Model I
      if(keymodel == 1)then
         A = 0.01d0
         B = 1.6d0
         C = 0.005d0
         D = 1.0d0
         if(rxyz(1) > 0.0)then
             hel(1,1) = A*(1.0-dexp(-B*rxyz(1)))
             dhel(1,1,1) = A*B*dexp(-B*rxyz(1))
         else
             hel(1,1) = -A*(1.0-dexp(B*rxyz(1)))
             dhel(1,1,1) = A*B*dexp(B*rxyz(1))
         endif

         hel(2,2) = -1.d0*hel(1,1)
         hel(1,2) = C*dexp(-D*rxyz(1)*rxyz(1))
         hel(2,1) = hel(1,2)

         dhel(2,2,1) = -1.d0*dhel(1,1,1)
         dhel(1,2,1) = hel(1,2)*(-2.*D*rxyz(1))
         dhel(2,1,1) = dhel(1,2,1)

! if(keymodel == 2) <= Tully Model II
      elseif(keymodel == 2)then
         A = 0.10d0
         B = 0.28d0
         C = 0.015d0
         D = 0.06d0
         Eo = 0.05d0

         hel(2,2) = -A*dexp(-B*rxyz(1)*rxyz(1)) + Eo
         hel(1,2) = C*dexp(-D*rxyz(1)*rxyz(1))
         hel(2,1) = hel(1,2)

         dhel(2,2,1) = (hel(2,2)-Eo) * (-2.*B*rxyz(1))
         dhel(1,2,1) = hel(1,2) * (-2.*D*rxyz(1))
         dhel(2,1,1) = dhel(1,2,1)

! if(keymodel == 3) <= Tully Model III
      elseif(keymodel == 3)then
         A = 0.0006d0
         B = 0.10d0
         C = 0.90d0
      ! interchange sign of A here for currect initialazation of wavepacket
      ! into state1

         hel(1,1) = -A  
         hel(2,2) = A

         if(rxyz(1) < 0.0)then
            hel(1,2) = B*dexp(C*rxyz(1))
            dhel(1,2,1) = B*C*dexp(C*rxyz(1))
         else
            hel(1,2) = B*(2.d0-dexp(-C*rxyz(1)))
            dhel(1,2,1) = B*C*dexp(-C*rxyz(1))
         endif

         hel(2,1) = hel(1,2)
         dhel(2,1,1) = dhel(1,2,1)

! if(keymodel == 4) <= Morse Model I
      elseif(keymodel == 4)then
         Dm = (/0.003d0, 0.004d0, 0.003d0/)
         bm = (/0.650d0, 0.600d0, 0.650d0/)
         Re = (/5.000d0, 4.000d0, 6.000d0/)
         cm = (/0.000d0, 0.010d0, 0.006d0/)

         Aij = (/0.002d0, 0.000d0, 0.002d0/)
         Rij = (/3.400d0, 0.000d0, 4.800d0/)
         amj = (/16.00d0, 0.000d0, 16.00d0/)

         !Diagonal elements
         do i = 1,3
            hel(i,i) = Dm(i)*(1.0d0-dexp(-bm(i)*(rxyz(1)-Re(i))))**2 &
                       +cm(i)

            dhel(i,i,1) =2.d0*Dm(i)*bm(i)* &
                         (1.0d0-dexp(-bm(i)*(rxyz(1)-Re(i))))* &
                         dexp(-bm(i)*(rxyz(1)-Re(i)))
         enddo

         !Off-Diagonal elements
         hel(1,2) = Aij(1)*dexp(-amj(1)*(rxyz(1)-Rij(1))**2)
         hel(2,3) = Aij(3)*dexp(-amj(3)*(rxyz(1)-Rij(3))**2)

         hel(2,1) = hel(1,2)
         hel(3,2) = hel(2,3)

         dhel(1,2,1) = -2*amj(1)*(rxyz(1)-Rij(1)) * hel(1,2)
         dhel(2,3,1) = -2*amj(3)*(rxyz(1)-Rij(3)) * hel(2,3)

         dhel(2,1,1) = dhel(1,2,1)
         dhel(3,2,1) = dhel(2,3,1)

! if(keymodel == 5) <=  Morse Model II      
      elseif(keymodel == 5)then
         Dm = (/0.020d0, 0.010d0, 0.003d0/)
         bm = (/0.650d0, 0.400d0, 0.650d0/)
         Re = (/4.500d0, 4.000d0, 4.400d0/)
         cm = (/0.000d0, 0.010d0, 0.020d0/)

         Aij = (/0.005d0, 0.005d0, 0.000d0/)
         Rij = (/3.660d0, 3.340d0, 0.000d0/)
         amj = (/32.00d0, 32.00d0, 0.000d0/)

        !Diagonal elements
         do i = 1,3
            hel(i,i) = Dm(i)*(1.0d0-dexp(-bm(i)*(rxyz(1)-Re(i))))**2 &
                       +cm(i)

            dhel(i,i,1) = 2*Dm(i)*bm(i)*(1.0d0-dexp(-bm(i)*(rxyz(1)- &
                          Re(i))))*dexp(-bm(i)*(rxyz(1)-Re(i)))
         enddo

         !Off-Diagonal elements
         hel(1,2) = Aij(1)*dexp(-amj(1)*(rxyz(1)-Rij(1))**2)
         hel(1,3) = Aij(2)*dexp(-amj(2)*(rxyz(1)-Rij(2))**2)

         hel(2,1) = hel(1,2)
         hel(3,1) = hel(1,3)

         dhel(1,2,1) = -2*amj(1)*(rxyz(1)-Rij(1)) * hel(1,2)
         dhel(1,3,1) = -2*amj(2)*(rxyz(1)-Rij(2)) * hel(1,3)

         dhel(2,1,1) = dhel(1,2,1)
         dhel(3,1,1) = dhel(1,3,1)

! if(keymodel == 6) <=  Morse Model III      
      elseif(keymodel == 6)then
         Dm = (/0.020d0, 0.020d0, 0.003d0/)
         bm = (/0.400d0, 0.650d0, 0.650d0/)
         Re = (/4.000d0, 4.500d0, 6.000d0/)
         cm = (/0.020d0, 0.000d0, 0.020d0/)

         Aij = (/0.005d0, 0.005d0, 0.000d0/)
         Rij = (/3.400d0, 4.970d0, 0.000d0/)
         amj = (/32.00d0, 32.00d0, 0.000d0/)

         !Diagonal elements
         do i = 1,3
            hel(i,i) = Dm(i)*(1.0d0-dexp(-bm(i)*(rxyz(1)-Re(i))))**2 &
                       +cm(i)

            dhel(i,i,1) = 2*Dm(i)*bm(i)*(1.0d0-dexp(-bm(i)*(rxyz(1)- &
                          Re(i))))*dexp(-bm(i)*(rxyz(1)-Re(i)))
         enddo

         !Off-Diagonal elements
         hel(1,2) = Aij(1)*dexp(-amj(1)*(rxyz(1)-Rij(1))**2)
         hel(1,3) = Aij(2)*dexp(-amj(2)*(rxyz(1)-Rij(2))**2)

         hel(2,1) = hel(1,2)
         hel(3,1) = hel(1,3)

         dhel(1,2,1) = -2*amj(1)*(rxyz(1)-Rij(1)) * hel(1,2)
         dhel(1,3,1) = -2*amj(2)*(rxyz(1)-Rij(2)) * hel(1,3)

         dhel(2,1,1) = dhel(1,2,1)
         dhel(3,1,1) = dhel(1,3,1)

! if(model == 7) <= Spin-Boson Model, Debye  (J. Chem. Phys. 144, 094104 (2016)
      elseif(keymodel == 7)then

         hel = 0.d0     
         dhel = 0.d0

         Hba = 0.5d0*m*sum(Wb(:)*Wb(:)*rxyz(:)*rxyz(:))
         Hsb = sum(cb(:)*rxyz(:))

         hel(1,1) =  eps + Hba + Hsb
         hel(2,2) = -eps + Hba - Hsb
         !hel(1,1) =  eps + Hsb
         !hel(2,2) = -eps - Hsb
         hel(1,2) = delta
         hel(2,1) = delta

         do ip=1,np
            dhel(1,1,ip) = m * Wb(ip)**2 * rxyz(ip) + cb(ip)
            dhel(2,2,ip) = m * Wb(ip)**2 * rxyz(ip) - cb(ip)
         enddo


! if(keymodel == 12) <= detailed balance  2-State with N-Chain Model      
      elseif(keymodel == 12)then

         hel(1,1) = v11
         hel(2,2) = v22

! if(keymodel == 13) <= detailed balance 3-State SuperExchange with N-Chain Model      
      elseif(keymodel == 13)then

         hel(1,1) = v11 
         hel(2,2) = v22 
         hel(3,3) = v33

! if(model == 14) <=  3-State Super Exchange Model      
      elseif(keymodel == 14)then
         aa = (/0.d0,0.01d0,0.005d0/)
         bb = (/0.001d0,0.01d0,0.d0/)

         !Diagonal elements
         do i = 1,3
             hel(i,i) = aa(i)

             dhel(i,i,1) = 0.d0
         enddo

         !Off-Diagonal elements
         hel(1,2) = bb(1)*exp(-rxyz(1)*rxyz(1)*0.5d0)
         hel(2,3) = bb(2)*exp(-rxyz(1)*rxyz(1)*0.5d0)
         hel(1,3) = 0.d0

         hel(2,1) = hel(1,2)
         hel(3,2) = hel(2,3)
         hel(3,1) = hel(1,3)

         dhel(1,2,1) = -rxyz(1) * hel(1,2)
         dhel(2,3,1) = -rxyz(1) * hel(2,3)
         dhel(1,3,1) = -rxyz(1) * hel(1,3)

         dhel(2,1,1) = dhel(1,2,1)
         dhel(3,2,1) = dhel(2,3,1)
         dhel(3,1,1) = dhel(1,3,1)

! if(keymodel == 21) <= Flat BO model
      elseif(keymodel == 21)then
         A = 0.005d0
         C = 5.5d0
         D = 0.8d0
         B = C*pi*(tanh(D*rxyz(1))+1.d0)

         hel(1,1) = -A*dcos(B)
         hel(1,2) = A*dsin(B)
         hel(2,1) = hel(1,2)
         hel(2,2) = -hel(1,1)

         dhel(1,1,1) = pi*C*D* 1.d0/(dcosh(D*rxyz(1))*dcosh(D*rxyz(1)))*hel(1,2)
         dhel(1,2,1) = pi*C*D*1.d0/(dcosh(D*rxyz(1))*dcosh(D*rxyz(1)))*hel(2,2)
         dhel(2,1,1) = dhel(1,2,1)
         dhel(2,2,1) = -dhel(1,1,1)

      else
         write(0,*) 'Wrong Model System'
         write(0,*) 'Choose Proper Model System!!!'
         stop
      endif

      return
      end subroutine gethel
     

      SUBROUTINE DIAG(EVALUES,EVECT,CRV)
!**********************************************************************
!     
!     CRV: HERMITIAN MATRIX (INPUT)
!     EVECT: EIGENVECTORS (OUTPUT)
!     EVALUES: EIGENVALUES (OUTPUT)
! 
!**********************************************************************

      use global_module, only: nstates
      implicit none

      character :: JOBZ,UPLO
      integer   :: INFO
      integer   :: N,NAP,LDZ
      integer   :: I,J,IND
      real*8    :: AP(nstates*(nstates+1)/2),WORK(3*nstates)
      real*8    :: EVALUES(nstates)
      real*8    :: CRV(nstates,nstates),EVECT(nstates,nstates)
       
      N=nstates
      NAP=N*(N+1)/2
      LDZ=N

      EVALUES=0.
      EVECT=0.
      JOBZ='V' ! calculate both eigenvalue and eigenvector

      UPLO='L' ! lower diagonal matrix

      IND=0

      DO J=1,N
         DO I=J,N
            IND=IND+1
            AP(IND)=CRV(I,J)
         END DO
      END DO

      CALL DSPEV(JOBZ,UPLO,N,AP,EVALUES,EVECT,LDZ,WORK,INFO)

      RETURN
      
      END SUBROUTINE DIAG


!**********************************************************************
      end module models_module
