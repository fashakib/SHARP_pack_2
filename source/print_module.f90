      module print_module
!**********************************************************************
!     SHARP PACK module that contains printing variables     
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!     New Jersey Institute of Technology
!**********************************************************************
      use global_module
      use modelvar_module
      implicit none

      contains

      subroutine printin()
!**********************************************************************
!     SHARP PACK subroutine to print system model parameters     
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************

      implicit none

      integer, parameter :: nrite=2
      character(len=24)  :: datentime

      call fdate(datentime)

      open(nrite,file='param.out',status='unknown')

      write(nrite,*) '                                                                             '
      write(nrite,*) '*****************************************************************************'
      write(nrite,*) '                                                                             '
      write(nrite,*) '    *******      ***      ***        **         ********       ******        '
      write(nrite,*) '   **********    ***      ***       ****        **********     *********     '
      write(nrite,*) '  ***      ***   ***      ***      ******       ***     ***    ***     ***   '
      write(nrite,*) '  ***            ***      ***      *** ***      ***      ***   ***      ***  '
      write(nrite,*) '   *********     ************     ***   ***     ***    ****    ***     ***   '
      write(nrite,*) '     *********   ************     *********     *** *****      *********     '
      write(nrite,*) '            ***  ***      ***    ***********    *** ****       ********      '
      write(nrite,*) '  ***      ***   ***      ***    ***     ***    ***    ***     ***           '
      write(nrite,*) '   **********    ***      ***   ***       ***   ***     ***    ***           '
      write(nrite,*) '     *******     ***      ***  ****       ****  ***      ****  ***    Pack 2 '
      write(nrite,*) '                                                                             '
      write(nrite,*) '*****************************************************************************'
      write(nrite,*) '                                                                             '
      write(nrite,'(1x,A,2x,A)') '                              ',trim(method)                             
      if(ldboc)then
        if(dvkey.eq.1)then
          write(nrite,*)'     with diagonal Born-Oppenheimer correction in centroid level'
        elseif(dvkey.eq.2)then
          write(nrite,*)'       with diagonal Born-Oppenheimer correction in bead level'
        endif
        if(apkey.eq.1)then
          write(nrite,*)'        and probablity of transition also in centroid level'
        elseif(apkey.eq.2)then
          write(nrite,*)'          and probablity of transition also in bead level'
        endif
      endif
      if((apkey.eq.1).and.(nb.gt.1))then
      write(nrite,*) '                        Centroid Approximation                               '
      elseif((apkey.eq.2).and.(nb.gt.1))then
      write(nrite,*) '                          Bead Approximation                                 '
      endif
      write(nrite,*) '                    SIMULATION CONTROL PARAMETERS                            '
      write(nrite,*) '                                                                             '
      if(lfft)then
      write(nrite,*) '                              with FFT                                       '
      endif
      write(nrite,*) '*****************************************************************************'
      write(nrite,*) '                                                 '
      if(ncpu == 1)then
       write(nrite,*)'   Serial Execution on 1 CPUs                    '
      elseif(ncpu > 1)then
       write(nrite,*)'   Parallel Execution: Running on',ncpu,' core(s)   '
      endif
      write(nrite,*) '                                                 '
      write(nrite,*) "JOB STARTED AT: ", datentime   
      write(nrite,*) '                                                 '

      write(nrite,'(/,1x,A,2x,A)')     'MODEL                               ', modelname

      write(nrite,'(/,1x,A,2x,I8)')    'number of states                    ', nstates

      write(nrite,'(/,1x,A,2x,I8)')    'number of particle(s)               ', np

      write(nrite,'(/,1x,A,2x,I8)')    'number of nbead(s)                  ', nb

      write(nrite,'(/,1x,A,2x,I8)')    'number of ntraj (in each CPU)       ', ntraj

      if(nsample .ne. 0)then
         write(nrite,'(/,1x,A,2x,I8)') 'number of sampling simulation steps ', nsample
         write(nrite,'(/,1x,A,I8,A)') 'with',nequil,' pimd equlibrium steps'
         write(nrite,'(/,1x,A,I7,A,I5,A)') 'and',ntraj,' sample(s) at every ',nskip,' steps'
         if(lpimd)then
           write(nrite,'(/,1x,A,2x,f12.3)') 'pimd sampling simulation with taut   ', taut
           write(nrite,'(/,1x,A,2x,f8.2,1x,A)') 'and PILE thermostat at temperature  ',temp/beta,'K'
         else
           write(nrite,'(/,1x,A,2x,f8.2,1x,A)') 'pimd sampling w/o thermostat at temperature  ',temp/beta,'K'
         endif
      endif

      if(ltrpmd)then
        write(nrite,'(/,1x,A,2x,f8.3)') 'thermostated RPMD simulation for RPSH dynamics '
      endif

      write(nrite,'(/,1x,A,2x,I8)')    'number of nsteps                    ', nsteps

      write(nrite,'(/,1x,A,2x,f8.2,A)')'R0 at which wave-function excited   ', R0,' (a.u.)'

      write(nrite,'(/,1x,A,2x,f8.2,A)')'P0, initial momentum, k             ', P0,' (a.u.)'

      if(R0key == 0)then
         write(nrite,'(/,1x,A,4x,A)')  'R0 beads for specific particle      ', "'fraction of de-Broglie length'"
      elseif(R0key == 1)then
         write(nrite,'(/,1x,A,4x,A)')  'R0 beads for specific particle      ', "'Gaussian'"
      elseif(R0key == 2)then
         write(nrite,'(/,1x,A,4x,A)')  'R0 beads for specific particle      ', "'Wigner'"
      endif

      if(R0key == 1)then
        if(lnorm .and. ldsurf)then
          write(nrite,'(/,1x,A)') "normal mode ring polymer sampling, and centroid mapping to dividing surface"
        elseif(lnorm .and. (.not. ldsurf))then
          write(nrite,'(/,1x,A)')"normal mode ring polymer sampling, and centroid NOT map to dividing surface"
        elseif((.not.lnorm) .and. ldsurf)then
          write(nrite,'(/,1x,A)') "direct sampling, and centroid map to dividing surface"
        elseif((.not.lnorm) .and. (.not.ldsurf))then
          write(nrite,'(/,1x,A)') "direct sampling, and centroid NOT map to dividing surface"
        endif
      else
        if(ldsurf)then
          write(nrite,'(/,1x,A)') "direct sampling, and centroid map to dividing surface"
        elseif(.not.ldsurf)then
          write(nrite,'(/,1x,A)') "direct sampling, and centroid NOT map to dividing surface"
        endif
      endif

      if(v0key == 0)then
         write(nrite,'(/,1x,A,4x,A)')  'P0 initializatoin scheme            ', "'Deterministic'"
      elseif(v0key == 1)then
         write(nrite,'(/,1x,A,4x,A)')  'P0 initialization scheme            ', "'Gaussian'"
      elseif(v0key == 2)then
         write(nrite,'(/,1x,A,4x,A)')  'P0 initialization scheme            ', "'Wigner'"
      endif

      if(vrkey == 0)then
         write(nrite,'(/,1x,A,4x,A)')  "Velocity Reversal for Frustrated Hop    'Never' "
      elseif(vrkey == 1)then
         write(nrite,'(/,1x,A,4x,A)')  "Velocity Reversal for Frustrated Hop    'Always' "
      elseif(vrkey == 2)then
         write(nrite,'(/,1x,A,4x,A)')  "Velocity Reversal for Frustrated Hop    'delV1' "
      elseif(vrkey == 3)then
         write(nrite,'(/,1x,A,4x,A)')  "Velocity Reversal for Frustrated Hop    'delV2' "
      endif

      if(vckey == 1)then
         write(nrite,'(/,1x,A,4x,A)')  "Velocity Rescaling Scheme               'Centroid' "
      elseif(vckey == 2)then
         write(nrite,'(/,1x,A,4x,A)')  "Velocity Rescaling Scheme               'Beads-Average' "
      endif

      if(dckey == 0)then
         write(nrite,'(/,1x,A,4x,A)')  "Decoherence Scheme                      'None' "
      elseif(dckey == 1)then
         write(nrite,'(/,1x,A,4x,A)')  "Decoherence Scheme                      'Exponential damping' "
      endif

      write(nrite,'(/,1x,A,2x,64f15.4)')  'mass of particle, mp (a.u.)         ', mp(1,1:nb)
      write(nrite,'(/,1x,A,/,(4e15.6))')  'normal mode freqeuncies,wk (a.u.): ', nfreq

      write(nrite,'(/,1x,A,2x,e13.6)') 'Inverse temperature of system, β    ', beta

      if(keymodel >=4 .and. keymodel <=6)then
         write(nrite,'(/,1x,A,2x,e13.6)') 'Vibrational frequency, ω            ', omega
      endif

      write(nrite,'(/,1x,A,2x,f8.4)')  'Simulation time step, dt (a.u.)     ', dt

      write(nrite,'(/,1x,A,2x,f8.4)')  'Electronic time step, dtE (a.u.)    ', dtq

      write(nrite,'(/,1x,A,2x,I8)')    'print skip                          ', iskip
      if(iseed.ne.0)then
        write(nrite,'(/,1x,A,2x,I12)') 'Random seed                         ', iseed
      endif

      if(ldboc)then
        write(nrite,'(/,1x,A,2x,f8.4)')  'DBOC force from finite difference method with eps (a.u.)', eps_dboc
      endif

      if(keymodel==7)then
       if(keybath==1)write(nrite,'(/,1x,A,2x,e13.6)')   'SPIN-BOSON (DEBYE) PARAMATERS:'
       if(keybath==2)write(nrite,'(/,1x,A,2x,e13.6)')   'SPIN-BOSON (OHMIC) PARAMATERS:'
         write(nrite,'(/,1x,A,2x,f8.4)')    '   ϵ             ', eps/enu
         write(nrite,'(/,1x,A,2x,f8.4)')    '   Δ             ', delta/enu
         write(nrite,'(/,1x,A,2x,f8.4)')    '   E_r           ', E_r/enu
         write(nrite,'(/,1x,A,2x,f8.4)')    '   ωc            ', wc/enu
!         write(nrite,'(/,1x,A,2x,f8.4)')    '   ωmax          ', wmax/enu
         write(nrite,'(/,1x,A,2x,f8.4)')    '   KT(1/β)       ',1.d0/(beta*enu)
         write(nrite,'(/,1x,A,2x,e13.6,2x,e13.6)')    '   enu           ', enu*freq, enu
         write(nrite,'(/,1x,A,2x,e13.6)') 'FINISH SPIN-BOSON PARAMATERS.'
      endif

      if(keymodel == 12)then
         write(nrite,'(/,1x,A,2x,e13.6)') '2-State LinearChain MODEL PARAMATERS:'
         write(nrite,'(/,1x,A,2x,e13.6)') '   T (in Kelvin)   ', kT
         write(nrite,'(/,1x,A,2x,e13.6)') '   V11 (kJ/mol)    ', v11/kJ_mol2au
         write(nrite,'(/,1x,A,2x,e13.6)') '   V22 (kJ/mol)    ', v22/kJ_mol2au
         write(nrite,'(/,1x,A,2x,e13.6)') '   d_12 (A^-1))    ', d_ab(1,2)/0.52918d0
         write(nrite,'(/,1x,A,2x,I4)')    '   N               ', np
         write(nrite,'(/,1x,A,2x,e13.6)') '   m (amu)         ', mp(1,1)/amu2au
         write(nrite,'(/,1x,A,2x,e13.6)') '   Vo (kJ/mol)     ', 175.d0
         write(nrite,'(/,1x,A,2x,e13.6)') '   a (A^-1)        ', 4.d0
         write(nrite,'(/,1x,A,2x,e13.6)') '   γ (s^-1)        ', 1.0e14
         write(nrite,'(/,1x,A,2x,e13.6)') 'FINISH 2-State LinearChain PARAMATERS.'
      endif
      if(keymodel == 13)then
         write(nrite,'(/,1x,A,2x,e13.6)') '3-State LinearChain MODEL PARAMATERS:'
         write(nrite,'(/,1x,A,2x,e13.6)') '   T (in Kelvin)   ', kT
         write(nrite,'(/,1x,A,2x,e13.6)') '   V11 (a.u.)      ', v11
         write(nrite,'(/,1x,A,2x,e13.6)') '   V22 (a.u.)      ', v22
         write(nrite,'(/,1x,A,2x,e13.6)') '   V33 (a.u.)      ', v33
         write(nrite,'(/,1x,A,2x,e13.6)') '   d_12 (a.u.)     ', d_ab(1,2)
         write(nrite,'(/,1x,A,2x,e13.6)') '   d_23 (a.u.)     ', d_ab(2,3)
         write(nrite,'(/,1x,A,2x,e13.6)') '   d_13 (a.u.)     ', d_ab(1,3)
         write(nrite,'(/,1x,A,2x,I4)')    '   N               ', np
         write(nrite,'(/,1x,A,2x,e13.6)') '   m (amu)         ', mp(1,1)/amu2au
         write(nrite,'(/,1x,A,2x,e13.6)') '   Vo (kJ/mol      ', 175.d0
         write(nrite,'(/,1x,A,2x,e13.6)') '   a (A^-1)        ', 4.d0
         write(nrite,'(/,1x,A,2x,e13.6)') '   γ (s^-1)        ', 1.0e14
         write(nrite,'(/,1x,A,2x,e13.6)') 'FINISH 3-State LinearChain PARAMATERS.'
      endif

      write(nrite,*) '                                                 '
      write(nrite,*) '*************************************************************************'

      close(nrite)

      return

      end subroutine


      subroutine printresults(lkval)
!**********************************************************************
!     SHARP PACK subroutine to print detail results of populations 
!     by different methods
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************

      implicit none
      
      integer :: i, j, itime, istart
      logical :: lkval,lfile

        open(56,file='pop_adiabat1.out',status='unknown')
        open(66,file='pop_adiabat2.out',status='unknown')
        open(76,file='pop_diabat1.out',status='unknown')
        open(77,file='pop_diabat2.out',status='unknown')
        open(78,file='pop_diabat3.out',status='unknown')

        write(56,'(A13,2x,A)') '# MODEL:    ', modelname
        write(66,'(A13,2x,A)') '# MODEL:    ', modelname
        write(76,'(A13,2x,A)') '# MODEL:    ', modelname
        write(77,'(A13,2x,A)') '# MODEL:    ', modelname
        write(78,'(A13,2x,A)') '# MODEL:    ', modelname

        write(56,'(A13,2x,A)') '#TIME (a.u.)', 'ADIABATIC POPULATIONS: METHOD I'
        write(66,'(A13,2x,A)') '#TIME (a.u.)', 'ADIABATIC POPULATIONS: METHOD II'
        write(76,'(A13,2x,A)') '#TIME (a.u.)', 'DIABATIC POPULATIONS: METHOD I'
        write(77,'(A13,2x,A)') '#TIME (a.u.)', 'DIABATIC POPULATIONS: METHOD II'
        write(78,'(A13,2x,A)') '#TIME (a.u.)', 'DIABATIC POPULATIONS: METHOD III'
        
!       if(ntrajR > 0)then 
!        open(86,file='pop_adiabatR1.out',status='unknown')
!        open(87,file='pop_adiabatR2.out',status='unknown')
!        open(88,file='pop_diabatR3.out',status='unknown')

!        write(86,'(A13,2x,A)') '# MODEL:    ', modelname
!        write(87,'(A13,2x,A)') '# MODEL:    ', modelname
!        write(88,'(A13,2x,A)') '# MODEL:    ', modelname
        
!        write(86,'(A13,2x,A)') '#TIME (a.u.)', 'ADIABATIC REFLECTED POPULATIONS: METHOD I'
!        write(87,'(A13,2x,A)') '#TIME (a.u.)', 'ADIABATIC REFLECTED POPULATIONS: METHOD II'
!        write(88,'(A13,2x,A)') '#TIME (a.u.)', 'DIABATIC REFLECTED POPULATIONS: METHOD III'
!       endif

        istart=1
        if(iskip.gt.1)istart=0
        do itime=istart,nprint
           write(56,222) itime*iskip*dt,(redmat(i,itime)/real(ntraj),i=1,nstates)
           write(66,222) itime*iskip*dt,((redmat_ec(i,j,itime)/real(ntraj),i=1,nstates),j=1,nstates)
           write(76,222) itime*iskip*dt,(diabat1(i,itime)/real(ntraj),i=1,nstates)
           write(77,222) itime*iskip*dt,(diabat2(i,itime)/real(ntraj),i=1,nstates)
           write(78,222) itime*iskip*dt,(diabat3(i,itime)/real(ntraj),i=1,nstates)
         
!          if(ntrajR > 0)then 
!           write(86,222) itime*iskip*dt,(redmatR(i,itime)/real(ntraj),i=1,nstates)
!           write(87,222) itime*iskip*dt,((redmat_ecR(i,j,itime)/real(ntraj),i=1,nstates),j=1,nstates)
!           write(88,222) itime*iskip*dt,(diabat3R(i,itime)/real(ntraj),i=1,nstates)
!          endif
        end do

        close(56)
        close(66)
        close(76)
        close(77)
        close(78)

        if(ntrajR > 0)then 
        close(86)
        close(87)
        close(88)
        endif

      if(lkval)then
        
        inquire(file='pop_branch.out', exist=lfile)
        if(lfile)then
           open(3,file='pop_branch.out',status='old',access='append')
           open(4,file='pop_branch_diabat.out',status='old',access='append')
        else
           open(3,file='pop_branch.out',status='new')
           write(3,'(A36,2x,A)') '# BRANCHING PROBABILTY OF MODEL:    ', modelname
           write(3,'(A,A12,9(A15))') '#','k-val','T1-pop','T2-pop','R1-pop','R2-pop',&
                   'nRevTraj','nFrust','nFrustRev','nAttempted','nSuccess'

           open(4,file='pop_branch_diabat.out',status='new')
           write(4,'(A36,2x,A)') '# BRANCHING PROBABILTY OF MODEL:    ', modelname
           write(4,'(A,A12,9(A15))') '#','k-val','T1-pop','T2-pop','R1-pop','R2-pop',&
                   'nRevTraj','nFrust','nFrustRev','nAttempted','nSuccess'
        endif

        write(3,222) P0, ((redmat(i,nprint)-redmatR(i,nprint))/(ntraj), i=1,nstates),&
                (redmatR(i,nprint)/(ntraj),i=1,nstates),&
                real(ntrajR)/ntraj,real(nfrust_hop)/ntraj,real(nfrust_hop2)/ntraj,&
                real(sum(nJump)+sum(nJumpFail))/ntraj,real(sum(nJump))/ntraj

        write(4,222) P0, ((diabat1(i,nprint)-diabat1R(i,nprint))/(ntraj), i=1,nstates),&
                (diabat1R(i,nprint)/(ntraj),i=1,nstates),&
                real(ntrajR)/ntraj,real(nfrust_hop)/ntraj,real(nfrust_hop2)/ntraj,&
                real(sum(nJump)+sum(nJumpFail))/ntraj,real(sum(nJump))/ntraj

        close(3)
        close(4)
      endif

222 format(120(e13.5E3,2x))
!222 format(f15.6,100(2x,f13.8))
      return

      end subroutine


      subroutine printlogo(nrite)
!**********************************************************************
!     SHARP PACK subroutine to print logo
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************

      implicit none

      integer            :: nrite
      !character(len=*)  :: ofile

      !open(nrite,file=ofile,status='unknown')

      write(nrite,*) '                                                                             '
      write(nrite,*) '*****************************************************************************'
      write(nrite,*) '                                                                             '
      write(nrite,*) '    *******      ***      ***        **         ********       ******        '
      write(nrite,*) '   **********    ***      ***       ****        **********     *********     '
      write(nrite,*) '  ***      ***   ***      ***      ******       ***     ***    ***     ***   '
      write(nrite,*) '  ***            ***      ***      *** ***      ***      ***   ***      ***  '
      write(nrite,*) '   *********     ************     ***   ***     ***    ****    ***     ***   '
      write(nrite,*) '     *********   ************     *********     *** *****      *********     '
      write(nrite,*) '            ***  ***      ***    ***********    *** ****       ********      '
      write(nrite,*) '  ***      ***   ***      ***    ***     ***    ***    ***     ***           '
      write(nrite,*) '   **********    ***      ***   ***       ***   ***     ***    ***           '
      write(nrite,*) '     *******     ***      ***  ****       ****  ***      ****  ***           '
      write(nrite,*) '                                                                             '
      write(nrite,*) '*****************************************************************************'
      write(nrite,*) '                                                                             '

      return 
      end subroutine


      subroutine openfile()
!**********************************************************************
!     SHARP PACK subroutine to print open files
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none
  
      character(len=24) :: file_hopp, file_dcoup
      character(len=24) :: file_psi, file_dboc

      file_hopp = 'hoppinghist.out'
      file_dcoup = 'dcoupling.out'
      file_psi = 'psi.out'
      file_dboc = 'energy_dboc.out'

      if(ldtl)then
         open(nrite_hopp,file=file_hopp,status='unknown')
         call printlogo(nrite_hopp)
         write(nrite_hopp,'(A13,2x,A)') '# MODEL:    ', modelname
      
         open(nrite_dcoup,file=file_dcoup,status='unknown')
         write(nrite_dcoup,'(A13,2x,A)') '# MODEL:    ', modelname
         write(nrite_dcoup,'(1x,A)') '# nTraj      nSteps      R_1 (a.u.)    Velocity_1 (a.u.)     &
           KE    Ering     PE    TotalEnergy    istate    Energy(istate)   dCoupling (nstates,nstates)   Temp(n*T)  Temp(cent)   Rgr'
  
         open(nrite_psi,file=file_psi,status='unknown')
          write(nrite_psi,'(1x,A)') '# Traj  Time(au)  istate   rc(1)   psi(n,n)   c_i(n)   E_i(n)     d12'

      endif

      if(ldboc)then
        open(nrite_dboc,file=file_dboc,status='unknown')
        write(nrite_dboc,'(1x,A)') '# Traj   Time(au)    Xc(1)    Vc(1)   BO_Energy(n)   DBOC(n)   Adia_Energy(n)'
      endif

      return

      end subroutine openfile


      subroutine hopping_stat()
!**********************************************************************
!     SHARP PACK subroutine for writing hopping statistics
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none

      integer :: i,j

      write(nrite_hopp,'("============= Hopping Statistics ==============")')
      do i = 1, nstates
        do j = 1, nstates
          if(i .ne. j)then
            write(nrite_hopp,'(" # Accepted jump from     ",I3," to ",I3," : ",I8)') i,j,nJump(i,j)
            write(nrite_hopp,'(" # Rejected jump from     ",I3," to ",I3," : ",I8)') i,j,nJumpFail(i,j)
          endif
        enddo
      enddo

      write(nrite_hopp,*)
      do i = 1, nstates
        do j = 1, nstates
          if(i .ne. j)then
            write(nrite_hopp,'(" # Frustrated Hop from    ",I3," to ",I3," : ",I8)') i,j,nFrust(i,j)
            write(nrite_hopp,'(" # FrustratedRev Hop from ",I3," to ",I3," : ",I8)') i,j,nFrustR(i,j)
          endif
        enddo
      enddo

      write(nrite_hopp,*)
      write(nrite_hopp,'(" #Total Successful Jump:  ",I8)') sum(nJump)
      write(nrite_hopp,'(" #Total Attempted Jump:   ",I8)') sum(nJump)+sum(nJumpFail)
      write(nrite_hopp,*)
      write(nrite_hopp,*) "#Total Initial states:   ", nIniStat
      write(nrite_hopp,*)
      write(nrite_hopp,*) "#Total Frustrated Hop:   ", nfrust_hop
      write(nrite_hopp,*) "#Reversed Frustrated Hop:", nfrust_hop2
      write(nrite_hopp,'("===============================================")')

      return

      end subroutine hopping_stat


      subroutine hopping_stat2()
!**********************************************************************
!     SHARP PACK subroutine for writing hopping statistics
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none

      integer :: i,j,nrite_hop2=130
      character(len=24)  :: datentime

      open(unit=nrite_hop2,file='param.out',status='old',access='append')

      write(nrite_hop2,'("============= Hopping Statistics ==============")')
      do i = 1, nstates
        do j = 1, nstates
          if(i .ne. j)then
            write(nrite_hop2,'(" # Accepted jump from     ",I3," to ",I3," : ",I8)') i,j,nJump(i,j)
            write(nrite_hop2,'(" # Rejected jump from     ",I3," to ",I3," : ",I8)') i,j,nJumpFail(i,j)
          endif
        enddo
      enddo

      write(nrite_hop2,*)
      do i = 1, nstates
        do j = 1, nstates
          if(i .ne. j)then
            write(nrite_hop2,'(" # Frustrated Hop from    ",I3," to ",I3," : ",I8)') i,j,nFrust(i,j)
            write(nrite_hop2,'(" # FrustratedRev Hop from ",I3," to ",I3," : ",I8)') i,j,nFrustR(i,j)
          endif
        enddo
      enddo

      write(nrite_hop2,*)
      write(nrite_hop2,'(" #Total Successful Jump:  ",I8)') sum(nJump)
      write(nrite_hop2,'(" #Total Attempted Jump:   ",I8)') sum(nJump)+sum(nJumpFail)
      write(nrite_hop2,*)
      write(nrite_hop2,*) "#Total Initial states:   ", nIniStat
      write(nrite_hop2,*)
      write(nrite_hop2,*) "#Total Frustrated Hop:   ", nfrust_hop
      write(nrite_hop2,*) "#Reversed Frustrated Hop:", nfrust_hop2
      write(nrite_hop2,'("===============================================")')
      write(nrite_hop2,*) '                                                 '
      write(nrite_hop2,*) "When using SHARP Pack, please cite the following papers:"
      write(nrite_hop2,*) "D. K. Limbu and F. A. Shakib, Software Impacts 19, 100604 (2024)."
      write(nrite_hop2,*) "D. K. Limbu and F. A. Shakib, J. Phys. Chem. Lett. 14, 8658–8666 (2023)."

      call fdate(datentime)

      write(nrite_hop2,*) '                                                 '
      write(nrite_hop2,*) "JOB FINISHED AT: ", datentime   
      write(nrite_hop2,*) '                                                 '
      write(nrite_hop2,'("===============================================")')

      close(nrite_hop2)

      return

      end subroutine hopping_stat2


      subroutine close_file()
!**********************************************************************
!     SHARP PACK subroutine to close files
!     
!     authors    - D.K. Limbu & F.A. Shakib     
!     copyright  - D.K. Limbu & F.A. Shakib
!
!     Method Development and Materials Simulation Laboratory
!**********************************************************************
      implicit none

      if(ldtl)close(nrite_hopp)
      if(ldtl)close(nrite_dcoup)
      if(ldtl)close(nrite_psi)
      if(ldboc)close(nrite_dboc)
      
      return 
      end subroutine close_file

!***********************************************************************
     subroutine printdetail()

     implicit none

     integer  :: ibd
     integer  :: itraj,itime


     end subroutine


!***********************************************************************

      end module print_module
