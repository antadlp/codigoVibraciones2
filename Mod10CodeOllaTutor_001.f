

* ********************************************************
*           Interactive Molecular Dynamics                *
*               Module: pressure cooker                   *
*                         by                              *
*                Dr. Ruben Santamaria                     *
*                 rso@fisica.unam.mx                      *
*            Dept. of Theoretical Physics                 *
*                Institute of Physics                     *
*               Univ. of Mexico, UNAM                     *
* ********************************************************
      program pressure_cooker
      implicit double precision (a-h,o-z)


      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)

      CHARACTER atsymbol(nmax)*5
      dimension atomass(nmax)

* positions

      dimension x(nmax),    y(nmax),    z(nmax)
      dimension xnew(nmax), ynew(nmax), znew(nmax)
      dimension Xn(nmax),   Yn(nmax),   Zn(nmax)
      dimension x0(nmax),   y0(nmax),   z0(nmax)

* speeds

      dimension vx(nmax),    vy(nmax),    vz(nmax)
      dimension vxnew(nmax), vynew(nmax), vznew(nmax)
      dimension Vxn(nmax),   Vyn(nmax),   Vzn(nmax)

* ab-initio gradients

      dimension Gradx(nmax),    Grady(nmax),    Gradz(nmax)
      dimension Gradxnew(nmax), Gradynew(nmax), Gradznew(nmax)

* kinetic & potential energies

      dimension enekin(-1:maxt), enepot(-1:maxt)

* temperatures, pressure, force constants

      dimension tCage(maxt), tConf(maxt), press(maxt), sK(nmax)

      dimension ftest(nmax), sum(nmax), knt(m1)

      CHARACTER imodule*40

      integer seed
      real*8 zbqlu01

* Otros

        CHARACTER (LEN = 400) RutaCarpetaDinamica 
        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaInp
        CHARACTER (LEN = 400) NombreOllaOut
        CHARACTER (LEN = 400) MovieOut
        CHARACTER (LEN = 400) NwchemScript
        CHARACTER (LEN = 400) Gradientes

* ------------------------------------------
* identify initial cpu processes & directories
* ------------------------------------------
c       Pido la ruta de la carpeta dinamica
      
      open (unit = 1, file ="DinamicaACorrer", status='old')
       
      read(1,*) RutaCarpetaDinamica

      close (1)


c       Abro y escribo los DatosGrlsDinamica
      open (unit=88, file ='.ckout')

      write(88,*) trim(RutaCarpetaDinamica)
      write(88,*) " "

      CLOSE (88)

c       LLamo a la rutina para los archivos de 
c       entrada y salida

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)


       NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

       Gradientes = trim(RutaArchivosSalida)
     & // '/Gradientes_' // trim(NombreDinamica)
     & // '.out'
 

       NwchemScript = trim(RutaArchivosEntrada)
     & // '/NwchemScript_' // trim(NombreDinamica)
       
        print*, NombreOllaOut
        open (unit = 77, file = trim(NombreOllaOut))
        open (unit = 59, file = trim(Gradientes))


      call infoUser


      write(77,*) "* ---------------------------------------------- *"
      write(77,*) "          INTERACTIVE MOLECULAR DYNAMICS          "
      write(77,*) "             Module:  Pressure Cooker"
      write(77,*)
      write(77,*) "  CITATION:"
      write(77,*) "  IMD program by Ruben Santamaria, IFUNAM, 2014"
      write(77,*) "* ---------------------------------------------- *"

      write(77,*) "internal modules (alphabetical order)"
      write(77,*) "analytic"
      write(77,*) "classify"
      write(77,*) "positions"
      write(77,*) "powerSer"
      write(77,*) "printGauss"
      write(77,*) "speeds"
      write(77,*)
      write(77,*) "external modules"
      write(77,*) "averageT"
      write(77,*) "breakBond"
      write(77,*) "cartesian"
      write(77,*) "citation"
      write(77,*) "distances"
      write(77,*) "eneNwchem"
      write(77,*) "eneSpring"
      write(77,*) "estimateCPU"
      write(77,*) "infoUser"
      write(77,*) "initConditOlla"
      write(77,*) "inptNwchem"
      write(77,*) "kinEne"
      write(77,*) "kparse"
      write(77,*) "largest"
      write(77,*) "makeMovie"
      write(77,*) "massCenter"
      write(77,*) "molInfo"
      write(77,*) "periodicTable"
      write(77,*) "perturb"
      write(77,*) "polar"
      write(77,*) "restartOlla"
      write(77,*) "scriptNwchem"
      write(77,*) "shrink"
      write(77,*) "slowSpeeds"
      write(77,*) "taylor"
      write(77,*) "units"
      write(77,*) "workDir"
      write(77,*) "writeMol"
      write(77,*) "zRndGenerator"

* -------------------------
* collect molecular info
* -------------------------

c      itest= 1  ! run a debug test
c      itest= 2  ! run the ab-initio calculations

      imodule= "Pressure Cooker (nwchem version)"

      call molInfo (npar,ncage, atsymbol,atomass, x,y,z,
     &   dt,Time, rfac,vfac, tempEquil, sK,xi,
     &   radiusf, ifreqShr,ifreqStor,irestDyn,iwavFunct,nsteps,
     &   itest)

c     call units
      write(77,*) "DESPUES DE MOLINFO"
      call inptNwchem (3,npar,atsymbol,x,y,z)  ! print input to the output file
      call scriptNwchem          ! create the script to run nwchem

      call powerSer (dt,xi, c0,c1,c2,c3, aa,bb,cc) ! parameters of the distrib funct

      call zbqlini(seed)

* restart the dynamics

      if (irestDyn.eq.2) then ! <----
                                     !
      call restartOlla (2,trest,npar, x,y,z, vx,vy,vz,
     &       Gradx,Grady,Gradz, x0,y0,z0, refPE)
                                     !
      write(77,*)                         !
      write(77,18) "restarting the dynamics from time step", trest
                                     !
                                     ! reuse the wave function
      call inptNwchem (2,npar,atsymbol,x,y,z)
                                     !
      goto 120                       !
                                     !
      end if ! <---------------------

* ----------------------------------------------
* impose initial conditions on positions & speeds
* ----------------------------------------------

      imzero= 0

      call massCenter (ncage,imzero,atomass,x,y,z,cmx,cmy,cmz)

      call initConditOlla (npar,ncage, atomass, x,y,z, vx,vy,vz,
     &    x0,y0,z0, cmx,cmy,cmz, tempEquil,
     &    nsteps,ibTime1,ibTime2,ibTime3)

* ----------------------------------------------
* write the molec, submit a job & get new forces
* ----------------------------------------------

      call perturb (npar,ncage, x,y,z, vx,vy,vz,
     &    xnew,ynew,znew, vxnew,vynew,vznew, rfac,vfac)

      t= -dt  ! present time

      write(77,*) "--------------------"
      write(77,"(1x,a,i6,1x,f9.2)") "step & time", -1, t

      call inptNwchem (iwavFunct,npar,atsymbol,
     & xnew,ynew,znew)  ! create the external input file

      if (itest.eq.2) then ! <----------
c       call system ("./script.nwchem") ! submit a job
        call system ("./"//trim(NwchemScript))  !
        call sleep (20)                 ! seconds to close files
                                        ! get new forces
        call eneNwchem (t,npar,Gradxnew,Gradynew,Gradznew,Epot,refPE)
      end if ! <------------------------

      call largest (t,npar,ncage, xnew,ynew,znew, radConf,radCage, Vol)
      call eneSpring (t,npar,ncage,xnew,ynew,znew,x0,y0,z0,sK,Espring)
      call kinEne (npar,ncage, atomass, vxnew,vynew,vznew,
     &     EkinConf,EkinCage, tempConf,tempCage,Vol,dynPress)

* 1 nano-metre= 10 Angst

      Volnano= Vol/1000.
      Etot= EkinConf + EkinCage + Epot + Espring

      write(77,11) "  energies [Hartree] t=  ", t
      write(77,12) "kinetic energy confined  ", EkinConf
      write(77,12) "kinetic energy cage      ", EkinCage
      write(77,12) "potential ab energy      ", Epot
      write(77,12) "spring energy            ", Espring
      write(77,12) "total energy             ", Etot
      write(77,13) "temperature cage [K]     ", tempCage
      write(77,13) "temperature conf [K]     ", tempConf
      write(77,13) "dynamic pressure [atm]   ", dynPress
      write(77,13) "confinement rad [Ang]    ", radConf
      write(77,13) "cage radius [Ang]        ", radCage
      write(77,13) "confinement vol [nanom^3]", Volnano
      write(77,*)

      do i= 1, npar ! <--
                         !
         x(i)= xnew(i)   ! shuffle variables
         y(i)= ynew(i)   !
         z(i)= znew(i)   !
                         !
        vx(i)= vxnew(i)  !
        vy(i)= vynew(i)  !
        vz(i)= vznew(i)  !
                         !
        Gradx(i)= Gradxnew(i)
        Grady(i)= Gradynew(i)
        Gradz(i)= Gradznew(i)

        write(59, *) Gradx(i)
                         !
      end do ! <---------

      close (59)

      cpu1= 0
      call estimateCPU (cpu1)

* ------------------------------------------
* evolve positions & speeds with Taylor from t=-dt to t= 0
* ------------------------------------------

      call taylor (npar, ncage, atomass,
     &   x,y,z, vx,vy,vz, Gradx,Grady,Gradz,
     &   xnew,ynew,znew, vxnew,vynew,vznew, dt)

      t= 0  ! present time

      write(77,*) "--------------------"
      write (77,"(1x,a,i6,1x,f9.2)") "step & time", 0, t

      call inptNwchem (2,npar,atsymbol,
     & xnew,ynew,znew)  ! reuse the wave function

      if (itest.eq.2) then ! <----------
c       call system ("./script.nwchem") ! submit a job
        call system ("./"//trim(NwchemScript))  !
        call sleep (09)                 ! seconds to close files
                                        ! get new forces
        call eneNwchem (t,npar,Gradxnew,Gradynew,Gradznew,Epot,refPE)
      end if ! <------------------------

      call largest (t,npar,ncage, xnew,ynew,znew, radConf,radCage, Vol)
      call eneSpring (t,npar,ncage,xnew,ynew,znew,x0,y0,z0,sK,Espring)
      call kinEne (npar,ncage, atomass, vxnew,vynew,vznew,
     &     EkinConf,EkinCage, tempConf,tempCage,Vol,dynPress)

      Volnano= Vol/1000.
      Etot= EkinConf + EkinCage + Epot + Espring

      write(77,11) "  energies [Hartree] t=  ", t
      write(77,12) "kinetic energy confined  ", EkinConf
      write(77,12) "kinetic energy cage      ", EkinCage
      write(77,12) "potential ab energy      ", Epot
      write(77,12) "spring energy            ", Espring
      write(77,12) "total energy             ", Etot
      write(77,13) "temperature cage [K]     ", tempCage
      write(77,13) "temperature conf [K]     ", tempConf
      write(77,13) "dynamic pressure [atm]   ", dynPress
      write(77,13) "confinement rad [Ang]    ", radConf
      write(77,13) "cage radius [Ang]        ", radCage
      write(77,13) "confinement vol [nanom^3]", Volnano
      write(77,*)

      do i= 1, npar ! <--
                         !
        x(i)= xnew(i)    ! shuffle variables
        y(i)= ynew(i)    !
        z(i)= znew(i)    !
                         !
        vx(i)= vxnew(i)  !
        vy(i)= vynew(i)  !
        vz(i)= vznew(i)  !
                         !
        Gradx(i)= Gradxnew(i)
        Grady(i)= Gradynew(i)
        Gradz(i)= Gradznew(i)
                         !
      end do ! <---------

c     goto 100

* ---------------
* estimate the cpu
* ---------------

      cpu2= 0
      call estimateCPU (cpu2)

* cpu time per step

      cputime= (cpu1 + cpu2)/ 2.

* ab-initio computations= nsteps + 2
* sleeping time= nsteps*2 + 30
* radiation + charge uptake + sleeping time computations= 6*cputime+6*2

      cpu= cputime*(nsteps+2.) + nsteps*2+30 + 6*cputime+6*2  ! in sec

      ndays=   int(cpu/ (3600.*24.))
      nhours=  int(cpu/3600. -ndays*24.)
      minutes= int(cpu/60 -ndays*24.*60.-nhours*60.)
      nsec=    int(mod(cpu,60.))

      write(77,14) "number of time steps to evaluate", nsteps
      write(77,20) ndays, nhours, minutes, nsec
      write(77,*)

* --------------------------
* initialization of parameters
* --------------------------

  120 sEkin=  0.d0   ! sum of kinetic energies
      sEpot=  0.d0   ! sum of pot energies
      sEtot=  0.d0   ! sum of total energies
      sEtot2= 0.d0   ! sum of squares of kin energies
      stCage= 0.d0   ! sum of cage temperatures
      stConf= 0.d0   ! sum of confined atoms temperatures
      sPres=  0.d0   ! sum of pressures
      sVol=   0.d0   ! sum of volumes

      do k= 1, nsteps
       enekin(k)= 0.d0
       enepot(k)= 0.d0
       tConf(k)=  0.d0
       tCage(k)=  0.d0
       press(k)=  0.d0
      end do

      do i= 1, ncage
       sum(i)= 0.d0  ! summations for every atom speed
      end do

      icount= 0  ! counter of points for the histogram
      kount=  0  ! counter of frames of the simulation
      seed= 15

c     xl=  0.01       ! scan speeds in = [0,xl]
      xl=  0.03       ! an ordinary Maxwell distrib. funct.
      wid= xl/1000.   ! sub-intervals to classify speeds

      nintervals= int(xl/wid) ! total num of sub-intervals
      kdim= 0

      do i= 1, nintervals  ! initialize counters
       kdim= kdim +1
       knt(i)= 0
      end do

      write(77,18) "interval to scan atom speeds [0,xl]", 0.0, xl
      write(77,16) "resolution to classify speeds", wid
      write(77,14) "number of intervals in the histogram", kdim

* open the file for the animation frames

        MovieOut = trim(RutaArchivosSalida)
     & // '/MovieOut_' // trim(NombreDinamica)
     & // '.out'


c     open (unit=12, file= 'movie.xyz',status= 'unknown')
      open (unit = 92, file = trim(MovieOut), status = 'unknown')

      t= 0.d0
      iAft= 10  ! time after breaking bond

* ----------------------
      write(77,*) "evolving coords and speeds with velVerlet"
* ----------------------

      do k= 1, nsteps, 1 ! <-------------
                                         !
      t= t + dt                          !
      kount= kount + 1                   !
                                         !
      write(77,*) "--------------------" !
      write(77,"(1x,a,i6,1x,f9.2)") "step & time", k, t
                                         !
                                         !
                                ! -----------------
                                ! get new positions
                                         !
      call positions (t,npar,ncage, atomass,
     &  x,y,z, vx,vy,vz, Gradx,Grady,Gradz,
     &  xnew,ynew,znew, Xn,Yn,Zn, x0,y0,z0, sK,
     &  tempEquil, c1,c2,cc, dt)         !
                                         !
                                         !
                                 ! --------------
                                 ! get new forces
                                         !
                                         !
c     if  (k.le.ibTime1)                 !                       ! initial state
               call inptNwchem (2,npar,atsymbol,
     & xnew,ynew,znew)
c     if ((k.ge.ibTime1+1).and.(k.le.ibTime1+iAft))              ! excited state
c    &         call inptNwchem (5,npar,atsymbol,xnew,ynew,znew)
c     if ((k.ge.ibTime1+iAft+1).and.(k.le.ibTime2))              ! charge recovered
c    &         call inptNwchem (2,npar,atsymbol,xnew,ynew,znew)
c                                        !
c                                        !
c     if ((k.ge.ibTime2+1).and.(k.le.ibTime2+iAft))
c    &         call inptNwchem (5,npar,atsymbol,xnew,ynew,znew)
c     if ((k.ge.ibTime2+iAft+1).and.(k.le.ibTime3))
c    &         call inptNwchem (2,npar,atsymbol,xnew,ynew,znew)
c                                        !
c                                        !
c     if ((k.ge.ibTime3+1).and.(k.le.ibTime3+iAft))
c    &         call inptNwchem (5,npar,atsymbol,xnew,ynew,znew)
c     if  (k.ge.ibTime3+iAft+1)          !
c    &         call inptNwchem (2,npar,atsymbol,xnew,ynew,znew)
                                         !
                                         !
      if (itest.eq.2) then               !
c       call system ("./script.nwchem")  !
        call system ("./"//trim(NwchemScript))  !
        call sleep (02)                  !
        call eneNwchem (t,npar,Gradxnew,Gradynew,Gradznew,Epot,refPE)
      end if                             !
                                         !
                                         !
                                 ! --------------
                                 ! get new speeds
                                         !
      call speeds (t, npar,ncage, atomass,
     & x,y,z, vx,vy,vz, Gradx,Grady,Gradz, xnew,ynew,znew,
     & vxnew,vynew,vznew, Gradxnew,Gradynew,Gradznew,
     & Xn,Yn,Zn, Vxn,Vyn,Vzn, sK, x0,y0,z0,
     & tempEquil, c0,c1,c2, aa,bb,cc, dt, ftest)
                                         !
                                         !
      do i= 1, ncage                     ! count speeds for histogram
       icount= icount + 1                !
       sum(i)= sum(i) + ftest(i)         ! every cage atom has its
       point= ftest(i)                   ! own summation
       call classify (nintervals,point,knt,xl,wid)
      end do                             !
                                         !
                                         !
                                  ! ------------
                                  ! get energies
                                         !
      call largest (t,npar,ncage, xnew,ynew,znew, radConf,radCage, Vol)
      call eneSpring (t,npar,ncage,xnew,ynew,znew,x0,y0,z0,sK,Espring)
      call kinEne (npar,ncage, atomass, vxnew,vynew,vznew,
     &     EkinConf,EkinCage, tempConf,tempCage,Vol,dynPress)
                                         !
                                         !
      Etot= EkinConf + EkinCage + Epot + Espring
                                         !
                                         !
                               ! ------------------
                               ! store dynamic data
                                         !
      enekin(k)= EkinCage + EkinConf     !
      enepot(k)= Epot + Espring          !
      tCage(k)= tempCage                 !
      tConf(k)= tempConf                 !
      press(k)= dynPress                 !
                                         !
      Volnano= Vol/1000.                 ! vol in nano^3
                                         !
      write(77,11) "  energies [Hartree] t=  ", t
      write(77,12) "kinetic energy confined  ", EkinConf
      write(77,12) "kinetic energy cage      ", EkinCage
      write(77,12) "potential ab energy      ", Epot
      write(77,12) "spring energy            ", Espring
      write(77,12) "total energy             ", Etot
      write(77,13) "temperature cage [K]     ", tempCage
      write(77,13) "temperature conf [K]     ", tempConf
      write(77,13) "dynamic pressure [atm]   ", dynPress
      write(77,13) "confinement rad [Ang]    ", radConf
      write(77,13) "cage radius [Ang]        ", radCage
      write(77,13) "confinement vol [nanom^3]", Volnano
                                         !
                                         !
      sEkin=  sEkin + EkinCage +EkinConf ! sum kin energies
      sEpot=  sEpot + Epot + Espring     ! sum pot energies
      sEtot=  sEtot + Etot               ! sum tot energies
      sEtot2= sEtot2 + Etot**2           ! sum (tot energies)^2
      stCage= stCage + tempCage          ! sum cage temperatures
      stConf= stConf + tempConf          ! sum confined atoms temps
      sPres=  sPres + dynPress           ! sum pressures
      sVol=   sVol + Vol                 ! sum volumes
                                         !
                                         !
                               ! -----------------
                               ! shuffle variables
      do i= 1, npar                      !
         x(i)= xnew(i)                   !
         y(i)= ynew(i)                   !
         z(i)= znew(i)                   !
                                         !
        vx(i)= vxnew(i)                  !
        vy(i)= vynew(i)                  !
        vz(i)= vznew(i)                  !
                                         !
        Gradx(i)= Gradxnew(i)            !
        Grady(i)= Gradynew(i)            !
        Gradz(i)= Gradznew(i)            !
      end do                             !
                                         !
      call makeMovie (k,npar, atsymbol,x,y,z)
                                         !
                                         !
                                 ! --------------
                                 ! cage shrinking
                                         !
c     call distances (t,npar, x,y,z)     !
                                         !
      if (mod(k,ifreqShr).eq.0) then     !
      write(77,*) "shrinking the cage"       !
      call shrink (npar, ncage, x0,y0,z0, x,y,z, radiusf)
      end if                             !
                                         !
                               ! ------------------
                               ! implicit radiation
                                         !
                                         !
c     if ((k.eq.ibTime1).or.(k.eq.ibTime2).or.(k.eq.ibTime3)) then
c     write(77,"(/,a)") "ionizing radiation: bond stochastically broken"
c                                        !
c     call breakBond (npar,ncage,atsymbol,x,y,z)
c     call inptNwchem (4,npar,atsymbol,x,y,z) ! write new input file
c                                        !
c     ene1= Epot                         ! pot ene from previous step
c                                        !
c     if (itest.eq.2) then               ! submit a job
c     call system ("./script.nwchem")    !
c     call sleep (02)                    !
c     call eneNwchem (t,npar,Gradx,Grady,Gradz,Epot,refPE)
c     end if                             !
c                                        !
c     ene2= Epot                         ! in hartree
c     radEne= (ene2-ene1)* 27.2107       ! in eV
c     wLength= (4.55619348)/(ene2-ene1) *d-8 ! in meters
c                                        !
c     write(77,"(/,a,f9.5)") "ionization energy [eV]", radEne
c     write(77,"(a,f9.5)")   "energy wavelength [m] ", wLength
c     call makeMovie (k,npar, atsymbol,x,y,z)
c     end if                             !
                                         !
                                 ! -------------
                                 ! charge uptake
                                         !
                                         !
c     if ((k.eq.ibTime1+iAft).or.(k.eq.ibTime2+iAft).or.
c    &    (k.eq.ibTime3+iAft)) then      !
c     write(77,"(/,a)") "charge uptake from environment"
c                                        !
c     call inptNwchem (1,npar,atsymbol,x,y,z) ! write initial input file
c                                        !
c     if (itest.eq.2) then               ! submit a job
c     call system ("./script.nwchem")    !
c     call sleep (02)                    !
c     call eneNwchem (t,npar,Gradx,Grady,Gradz,Epot,refPE)
c     end if                             !
c     end if                             !
                                         !
                              ! -------------------
                              ! store relevant data
                                         !
      if (mod(k,ifreqStor).eq.0) then    !
      write(77,"(/,a)") "storing data for restart"
                                         !
      call restartOlla (1,t,npar, x,y,z, vx,vy,vz,
     &       Gradx,Grady,Gradz, x0,y0,z0, refPE)
                                         !
                                         !
      write(77,*)                             ! average temps
      call averageT (k,tCage,ifreqStor,averTemp)
      write(77,"(1x,a,f9.2)") "average cage temp [K]", averTemp
                                         !
      call averageT (k,tConf,ifreqStor,averTemp)
      write(77,"(1x,a,f9.2)") "average confined temp [K]", averTemp
      end if                             !
                                         !
      end do ! <-------------------------

c     close (12) ! close movie file
        close (92)

* -------------------------
* averages of main energies
* -------------------------

      Ekin=  sEkin  /dble(kount)    ! kin energy
      Epot=  sEpot  /dble(kount)    ! pot energy
      Etot=  sEtot  /dble(kount)    ! tot energy
      Etot2= sEtot2 /dble(kount)    ! (tot energy)^2
      TempCage= stCage /dble(kount) ! cage temp
      TempConf= stConf /dble(kount) ! confined atoms temp
      Pres=  sPres  /dble(kount)    ! pressure
      Vol=   sVol   /dble(kount)    ! volume

      Volnano= Vol/1000.

* energy fluctuation

      Efluct= dsqrt(dabs(Etot2 - Etot**2))/ Etot

* average speed of cage particles (~Maxwell
* average speed for the given temperature)

      sum2= 0.d0

      do i= 1, ncage
       v1= sum(i)/ dble(k)
       sum2= sum2 + v1
      end do

      avSpd= sum2/dble(ncage)

* info for the Vprom speed:
* http://hyperphysics.phy-astr.gsu.edu/hbase/kinetic/kintem.html

* [sqrt(RT/m)]= sqrt [J/(mol*Kelvin) * Kelvin/amu]= fac *Angst/fs

      fac= 1.d-5* sqrt(1000.)

      Vprom= fac* sqrt(2.*8.3144621 *tempEquil/atomass(1)) ! in Angst/fs

* create the speeds histogram and the analytic Maxwell speed distrib. funct.

      call printGauss (nintervals,icount,knt,xl,wid)
      call analytic (atomass,tempEquil,xl,wid)

      write(77,*)  "--------------------"
      write(77,14) "number of intervals in the histogram", kdim
      write(77,14) "number of classified speeds", icount
      write(77,16) "average rnd speed of cage particles [Angst/fs]",
     & avSpd
      write(77,16) "maxwellian speed [Angst/fs]", Vprom
      write(77,11) "maxwellian speed [m/s]", Vprom *1.d+5

      write(77,*)
      write(77,*) "the speed histogram of the cage particles was" 
     & // "created"
      write(77,*) "the analytic Maxwell speed distrib. funct."
     & //" was created"

* ------------------------------
* printing of physical variables
* ------------------------------

* print energies in terms of t

      write(77,*) "--------------------------------------------"
      write(77,*) "time   Ekin [Hartree] Epot [Hartree] Tcage [K]  "
     &//"Tconf [K] Pdyn [GPa]"

      t= 0.d0
      sum1= 0.
      sum2= 0.

      do k= 1, nsteps, 1 ! <----
                                !
       t= t + dt                !
       sum1= sum1 + tCage(k)    !
       sum2= sum2 + tCage(k)**2 !
                                !
       write(77,10) t, enekin(k),enepot(k), tCage(k),tConf(k),press(k)
                                !
      end do ! <----------------

      tempfluct= dsqrt(dabs(sum2 - sum1**2))/ sum1

      write(77,*)
      write(77,*) "number of computed time steps= ", kount
      write(77,*)

      write(77,*) "  average values [Hartree]"
      write(77,12) "kinetic ene             ", Ekin
      write(77,12) "potential ene           ", Epot
      write(77,12) "total ene               ", Etot
      write(77,12) "(total ene)^2           ", Etot2
      write(77,12) "ene fluctuat            ", Efluct
      write(77,13) "cage temperature [K]    ", TempCage
      write(77,13) "confined atoms temp [K] ", TempConf
      write(77,13) "cage temp fluctuat      ", tempfluct
      write(77,13) "dynamic pressure [atm]  ", Pres
      write(77,13) "volume of conf [nanom^3]", Volnano

      write(77,*) "--------------------------------------------"
      write(77,*) "final parameters of the dynamics"
      write(77,16) "time step [fs] ", dt
      write(77,16) "total simulation time [fs] ", Time
      write(77,16) "rfac      ", rfac
      write(77,16) "vfac      ", vfac
      write(77,14) "npar      ", npar
      write(77,14) "ncage     ", ncage
      write(77,16) "tempEquil ", tempEquil
      write(77,16) "spring K  ", sK(1)  ! all are the same
      write(77,16) "friction  ", xi
      write(77,16) "radiusf   ", radiusf
      write(77,14) "ifreqShr  ", ifreqShr
      write(77,14) "ifreqStor ", ifreqStor
      write(77,14) "irestDyn  ", irestDyn
      write(77,14) "iwavFunct ", iwavFunct
      write(77,*)

      write(77,*) "created files: movie.xyz, analytic, rndgauss"

      

  100 call citation (imodule)

   10 format (1x,f8.2,2(1x,d12.5), 3(1x,f10.3))
   11 format (1x,a,2x,f10.3)
   12 format (1x,a,2x,d16.8)
   13 format (1x,a,2x,f14.4)
   14 format (1x,a,2x,i10)
   16 format (1x,a,2x,f16.8)
   18 format (1x,a,2f9.5)
   20 format (1x,"shortest estimated cpu time for the dynamics",/,
     &     i3," days:", i3," hrs:", i3," min:", i3," sec")

* delete non-required files

c     call system ("rm mol.xyz a.out script.nwchem nwchem.inp")

      CLOSE (77)
      stop
      end
* *****************************
* create an analytic gaussian function

* info on the normal distribution function:
* http://hyperphysics.phy-astr.gsu.edu/hbase/kinetic/maxspe.html#c3
* http://en.wikipedia.org/wiki/Normal_distribution
* variance= sigma^2
* standard deviation= sigma
* *****************************
      subroutine analytic (atomass,tempEquil,xl,wid)
      implicit double precision (a-h,o-z)
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension atomass(nmax)
      real*8 kB

        
        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) RutaArchivosEntrada 
        CHARACTER (LEN = 400) AnalyticOut
        
        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
      
         AnalyticOut  = trim(RutaArchivosSalida)
     & // '/AnalyticOut_' // trim(NombreDinamica)
     & // '.out'
      
      open (unit = 89, file = trim(AnalyticOut), status = 'unknown')

c     open (16,file="analytic",status="unknown")

      pi= 3.1415926535898

* Boltzmann constant

      kB= 1.3806504d-23    ! in Joule/Kelvin

* facKT *(Angst/fs)^2= (Joule/Kelvin)*(Kelvin/amu)= J/amu

      facKT= 1.0d17/ 1.6605402

* normalization factor

      v1= atomass(1)/ (2.*kB*tempEquil) *1./facKT ! in (fs/Angst)^2
      fac= 4.*pi*(v1/pi)**(3./2.) ! in (fs/Angst)^3

* define intervals

      nn= 125
      dy= 0.0002
      y= -dy
      ik= 0

* loop over y-values

      do ki= 1, nn ! <---
                         !
      ik= ik + 1         !
      y= y + dy          !
                         ! in fs/Angst
      distrib= fac*y**2 *exp(-v1*y**2)
                         !
c     write (16,"(2f9.4)") y, distrib
      write (89,"(2f9.4)") y, distrib
                         !
      end do ! <---------

c     close (16)
      CLOSE (89)

      return
      end
* *****************************
* classify speeds in sub-intervals to build a histogram
* *****************************
      subroutine classify (nintervals,point,knt,xl,wid)
      implicit double precision (a-h,o-z)
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension knt(m1)  ! m1= ncage*Time~ 180*1000= 180 mil

      j= 0
      sum= 0.d0  ! values of point >= 0

      do ki= 1, nintervals ! <--
                                ! creat the j^th
      j=j+1                     ! sub-interval with
                                ! limits [sum,v1]
      sum= sum + wid            !
      v1=  sum + wid            !
                                !
c     write(77,"(i5,2f10.5)") j, sum, v1
                                !
                                ! is the point in the
                                ! interval [sum,v1)?
                                !
      if ((point.ge.sum).and.(point.lt.v1)) knt(j)= knt(j)+1
                                !
      end do ! <----------------

      return
      end
* *******************************
* evolve positions using vVerlet integrators
* *******************************
      subroutine positions (time, npar,ncage, atomass,
     &  x,y,z, vx,vy,vz, Gradx,Grady,Gradz,
     &  xnew,ynew,znew, Xn,Yn,Zn, x0,y0,z0, sK,
     &  tempEquil, c1,c2,cc, dt)

      implicit double precision (a-h,o-z)
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      external gauss

      dimension atomass(nmax)

      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)
      dimension Gradx(nmax), Grady(nmax), Gradz(nmax)

      dimension xnew(nmax), ynew(nmax), znew(nmax)
      dimension Xn(nmax), Yn(nmax), Zn(nmax)

      dimension x0(nmax), y0(nmax), z0(nmax), sK(nmax)

      real*8 kB

        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida      
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)

        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')


* Boltzmann constant

      kB= 1.3806504d-23    ! in Joule/Kelvin= Kg (m/s)^2 /Kelvin

* faccl * Angst= (1/amu)*(Hartree/Bohr)*(fs^2)

      faccl= 0.49616184    ! factors to change units

* fspring *(amu/fs^2) = Newton/cm

      fspring= 1.d0/ 16.605402

* ------------------
* evolve positions of cage atoms
* ------------------

      do i= 1, ncage ! <------------!----
                                    !
      amas= atomass(i)              !
                          ! -------------------
                          ! systematic forces
                                    !
      a1x= -faccl *Gradx(i) /amas   !
      a1y= -faccl *Grady(i) /amas   ! ab-initio forces
      a1z= -faccl *Gradz(i) /amas   !
                                    !
                                    ! F= -k*dX
      a2x= -fspring *sK(i) *(x(i)-x0(i)) /amas
      a2y= -fspring *sK(i) *(y(i)-y0(i)) /amas ! in Angst/fs^2
      a2z= -fspring *sK(i) *(z(i)-z0(i)) /amas
                                    !
      ax= a1x + a2x                 ! total
      ay= a1y + a2y                 !
      az= a1z + a2z                 !
                                    !
                          ! -------------------
                          ! stochastic forces
                                    !
      pkT= kB*tempEquil/amas *cc    ! in Angst^2
                                    !
      etax= gauss()                 ! random numbers
      etay= gauss()                 !
      etaz= gauss()                 ! stochastic positions
                                    ! Eqs 1.13.15 of chapter Langevin
      Xn(i)= etax * sqrt(pkT)       !
      Yn(i)= etay * sqrt(pkT)       ! in Angst
      Zn(i)= etaz * sqrt(pkT)       !
                                    !
                                    !
                          ! -------------------
                          ! position variables
                                    !
* x(t)= x(t-dt) + c1*v(t-dt)*dt + c2*F(t-dt)/m *dt^2 + Xn(dt)
                                    !
                                    ! in Angst
      xnew(i)= x(i) + c1*vx(i)*dt + c2*ax * dt**2 + Xn(i)
      ynew(i)= y(i) + c1*vy(i)*dt + c2*ay * dt**2 + Yn(i)
      znew(i)= z(i) + c1*vz(i)*dt + c2*az * dt**2 + Zn(i)
                                    !
                                    !
                          ! -------------------
                          ! debugging positions
                                    !
c     if ((time.le.1.0).and.(i.le.5)) then
c     if (i.le.3) then              !
c                                   !
c     write(77,15) i, sqrt(pkT)      !
c     write(77,51) "c1,c2, dt", c1,c2,dt
c     write(77,51) "aki posicion", x(i),c1*vx(i)*dt,c2*ax*dt**2,Xn(i)
c     write(77,51) "aki posicion", vx(i),vy(i),vz(i)
c                                   !
c  51 format (a,4f10.4)             !
c     end if                        !
                                    !
      end do   ! <------------------!-

* ------------------
* evolve positions of confined atoms
* ------------------

      do i= ncage +1, npar ! <------!----
                                    !
      amas= atomass(i)              !
                                    !
                          ! -------------------
                          ! systematic forces
                                    !
      ax= -faccl *Gradx(i) /amas    ! ab-initio forces
      ay= -faccl *Grady(i) /amas    !
      az= -faccl *Gradz(i) /amas    !
                                    !
                          ! -------------------
                          ! position variables
                                    !
* x(t)= x(t-dt) + v(t-dt)*dt + F(t-dt)/m *dt^2/2
                                    !
                                    ! in Angst
      xnew(i)= x(i) + vx(i)*dt + ax/2.d0 *dt**2
      ynew(i)= y(i) + vy(i)*dt + ay/2.d0 *dt**2
      znew(i)= z(i) + vz(i)*dt + az/2.d0 *dt**2
                                    !
      end do   ! <------------------!-
      return
      end
* ************************************
* power series of coefficients
* this module is called once
* ************************************
      subroutine powerSer (dt,xi, c0,c1,c2,c3, aa,bb,cc)
      implicit double precision (a-h,o-z)
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)

      CHARACTER (LEN = 400) NombreDinamica
      CHARACTER (LEN = 400) RutaArchivosEntrada
      CHARACTER (LEN = 400) RutaArchivosSalida
      CHARACTER (LEN = 400) NombreOllaOut

        
      call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
      NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out' 
c     open (unit = 77, file = trim(NombreOllaOut), status = 'old')

* xi in 1/ps,  1 ps= 1000 fs

      hxi= xi/1000.0   ! in 1/fs

      z= abs(hxi)*dt   ! dimensionless

* expansions from Paterlini's article (refer to the Appendix,
* some of the numbers there are mistaken)
* Chem. Phys. vol. 236, 243-252, 1998

* refer to expressions \label{eqmopafl.9} of my
* Langevin chapter for the parameters, c0, c1, etc

* -----------------------
* c0= e^{-z}  ;  c1= (1-c0)/z  ;  c2= (1-c1)/z  ;  c3= (1/2-c2)/z
* -----------------------

* the coefficients are dimensionless

      c0= 1.    -z     +z**2/2.  -z**3/6.  +z**4/24. -z**5/120.
      c1= 1.    -z/2.  +z**2/6.  -z**3/24. +z**4/120.
      c2= 1./2. -z/6.  +z**2/24. -z**3/120.
      c3= 1./6. -z/24. +z**2/120.

      write(77,*)
      write(77,*) "dynamics & distrib function factors"
      write(77,*) "(for internal use)"
      write(77,13) "c0=", c0
      write(77,13) "c1=", c1
      write(77,13) "c2=", c2
      write(77,13) "c3=", c3
      write(77,*)

* facKT *(Angst/fs)^2= (Joule/Kelvin)*(Kelvin/amu)

      facKT= 1.0d17/1.6605402

* -----------------------
* refer to expressions \ref{intteqom.3} of my
* Langevin chapter for the parameters aa, bb, cc
* -----------------------

* we include the factor (kB*Tref/mass)
* after the call to this module

      aa= 1. -z +2.*z**2/3. -z**3/3. +2.*z**4/15. -2.*z**5/45.
c     aa= 2.*kB*T/mass *z *aa *facKT   <------------ in (Angst/fs)^2
      aa= 2.*z*aa *facKT

      bb= 1.-z + 7.*z**2/12. -z**3/4. +31.*z**4/360. -z**5/40.
c     bb= kB*T/mass *z*dt *bb *facKT  <-------- in Angst^2/fs
      bb= z*dt*bb *facKT

      cc= 1.-3.*z/4.+7.*z**2/20.-z**3/8.+31.*z**4/840.-3.*z**5/320.
c     cc= 2.*kB*T/(3.*mass) *z*dt**2 *cc *facKT  <-------- in Angst^2
      cc= 2./3. *z*dt**2 *cc *facKT

* print without the factor facKT

      write(77,13) "a=", aa/facKT 
      write(77,13) "b=", bb/facKT
      write(77,13) "c=", cc/facKT
      write(77,*)

      write(77,14) "ac=", aa*cc/facKT**2
      write(77,14) "ac-b^2=", (aa*cc-bb**2)/facKT**2
      write(77,14) "sqrt((ac-b^2)/c)=", sqrt((aa*cc-bb**2)/(cc*facKT))
      write(77,14) "b/c=", bb/cc
      write(77,*)

   13 format (1x,a,1x,f10.4)
   14 format (1x,a,d13.4)


      return
      end
* *****************************
* print the gaussian curve for the random speeds
* *****************************
      subroutine printGauss (nintervals,npoints,knt,xl,wid)
      implicit double precision (a-h,o-z)
      !include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension knt(m1)

        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida      
        CHARACTER (LEN = 400) RndGaussOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
     
        RndGaussOut = trim(RutaArchivosSalida)
     & // '/RndGaussOut_' // trim(NombreDinamica)
     & // '.out'


      open (unit = 90, file = trim(RndGaussOut), status = 'unknown')
 

c     open (16,file="rndgauss", status="unknown")

* numeric integration of the function: base * height
* vnorm= wid*n1 + wid*n2 + wid*n3= dx*(n1+n2+...)= wid*npoints

      vnorm= wid*npoints  ! normalization factor
      j= 0
      sum= 0.d0

* print number of points at each sub-interval

      do ki= 1, nintervals ! <--
                                !
      j= j+1                    ! j^th interval
                                !
      sum= sum + wid            !
      v1=  sum + wid            ! get central point
      vc= (sum + v1)/ 2.        ! of the interval
                                !
c     write (16,"(2f12.3)") vc, dble(knt(j))/vnorm
      write (90,"(2f12.3)") vc, dble(knt(j))/vnorm
                                !
      end do ! <----------------

c     close (16)
      CLOSE (90)

      return
      end
* *******************************
* evolve speeds using vVerlet integrators
* *******************************
      subroutine speeds (time, npar,ncage, atomass,
     & x,y,z, vx,vy,vz, Gradx,Grady,Gradz,
     & xnew,ynew,znew, vxnew,vynew,vznew, Gradxnew,Gradynew,Gradznew,
     & Xn,Yn,Zn, Vxn,Vyn,Vzn, sK, x0,y0,z0,
     & tempEquil, c0,c1,c2, aa,bb,cc, dt,ftest)

      implicit double precision (a-h,o-z)
      !include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      external gauss

      dimension atomass(nmax), ftest(nmax)

      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)
      dimension Gradx(nmax), Grady(nmax), Gradz(nmax)

      dimension  xnew(nmax),  ynew(nmax),  znew(nmax)
      dimension vxnew(nmax), vynew(nmax), vznew(nmax)
      dimension Gradxnew(nmax), Gradynew(nmax), Gradznew(nmax)

      dimension  Xn(nmax),  Yn(nmax),  Zn(nmax)
      dimension Vxn(nmax), Vyn(nmax), Vzn(nmax)

      dimension x0(nmax), y0(nmax), z0(nmax), sK(nmax)

      real*8 kB


      
* Boltzmann constant

      kB= 1.3806504d-23    ! in Joule/Kelvin= Kg (m/s)^2 /Kelvin

* fvel *(Angst/fs)= (1/amu) *(Hartree/Bohr) *fs

      fvel= 0.49616184     ! factors to change units

* fspring *(amu/fs^2) = Newton/cm

      fspring= 1.d0/ 16.605402

* facKT *(Angst/fs)^2= (Joule/Kelvin)*(Kelvin/amu)

      facKT= 1.0d17/1.6605402

* ------------------
* evolve speeds of cage atoms
* ------------------

      do i= 1, ncage ! <------------
                                    !
      amas= atomass(i)              !
                                    !
                        ! ---------------------
                        ! old systematic forces
                                    !
      a1ax= -fvel *Gradx(i) /amas   ! ab-initio
      a1ay= -fvel *Grady(i) /amas   ! in Angst/fs^2
      a1az= -fvel *Gradz(i) /amas   !
                                    !
                                    ! F= -k*dX
      a1bx= -fspring* sK(i)* (x(i)-x0(i)) /amas
      a1by= -fspring* sK(i)* (y(i)-y0(i)) /amas ! in Angst/fs^2
      a1bz= -fspring* sK(i)* (z(i)-z0(i)) /amas
                                    !
                                    !
      a1x= a1ax + a1bx              ! total
      a1y= a1ay + a1by              !
      a1z= a1az + a1bz              !
                                    !
                        ! ---------------------
                        ! new systematic forces
                                    !
      a2ax= -fvel *Gradxnew(i) /amas! ab-initio
      a2ay= -fvel *Gradynew(i) /amas! in Angst/fs^2
      a2az= -fvel *Gradznew(i) /amas!
                                    !
                                    ! F= -k*dX
      a2bx= -fspring* sK(i)* (xnew(i)-x0(i)) /amas
      a2by= -fspring* sK(i)* (ynew(i)-y0(i)) /amas ! in Angst/fs^2
      a2bz= -fspring* sK(i)* (znew(i)-z0(i)) /amas
                                    !
                                    !
      a2x= a2ax + a2bx              ! total
      a2y= a2ay + a2by              !
      a2z= a2az + a2bz              !
                                    !
                          ! -----------------
                          ! stochastic speeds
                                    !
      pKT= kB*tempEquil/amas        !
      pa= pKT *aa                   ! (Angst/fs)^2
      pb= pKT *bb                   !  Angst^2/fs
      pc= pKT *cc                   !  Angst^2
                                    !
      etax= gauss()                 ! random numbers
      etay= gauss()                 ! Eqs 1.13.15 of chapter Langevin
      etaz= gauss()                 !
                                    ! in Angst/fs
      Vxn(i)= etax *sqrt(pa-pb**2/pc) +pb/pc *Xn(i)
      Vyn(i)= etay *sqrt(pa-pb**2/pc) +pb/pc *Yn(i)
      Vzn(i)= etaz *sqrt(pa-pb**2/pc) +pb/pc *Zn(i)
                                    !
                                    !
                                    !
                 ! -------------------------------------
                 ! reproduce a Maxwellian distrib funct.
                                    !
* http://hyperphysics.phy-astr.gsu.edu/hbase/kinetic/maxspe.html#c3
* run in test form, use the interval below for speeds and plot the
* analytic and rndgauss data        !
                                    !
*     xl=  0.03     ! interval to scan speeds = [-xl,xl]
*     wid= xl/1000. ! resolution of sub-intervals to scan speeds
                                    !
c     Vxn(i)= etax *sqrt(pKT*facKT) ! to reproduce
c     Vyn(i)= etay *sqrt(pKT*facKT) ! the Maxwell
c     Vzn(i)= etaz *sqrt(pKT*facKT) ! distrib. funct.
                                    !
                                    !
                                    ! test function just for rnd speeds
      ftest(i)= sqrt (Vxn(i)**2 + Vyn(i)**2 + Vzn(i)**2)
                                    !
                                    !
                            ! ---------------
                            ! speed variables
                                    !
* v(t)= c0*v(t-dt) + (c1-c2)*a(t-dt)*dt + c2*a(t)*dt + Vn(dt)
                                    !
                                    !
                                    ! in Angst/fs
      vxnew(i)= c0*vx(i) + (c1-c2)*a1x*dt + c2*a2x*dt +Vxn(i)
      vynew(i)= c0*vy(i) + (c1-c2)*a1y*dt + c2*a2y*dt +Vyn(i)
      vznew(i)= c0*vz(i) + (c1-c2)*a1z*dt + c2*a2z*dt +Vzn(i)
                                    !
                                    !
                           ! ----------------
                           ! debugging speeds
                                    !
c     if ((time.lt.9.65).and.(i.eq.1)) then
c     if (i.le.3) then              !
c                                   !
c     write(77,51) "aki vx   ",vx(i),vy(i),vz(i)
c     write(77,51) "aki a1ax ",a1ax,a1ay,a1az
c     write(77,51) "aki a1bx ",a1bx,a1by,a1bz
c     write(77,51) "aki a2ax ",a2ax,a2ay,a2az
c     write(77,51) "aki a2bx ",a2bx,a2by,a2bz
c     write(77,51) "aki Vxn  ",Vxn(i),Vyn(i),Vzn(i)
c     write(77,51) "aki vxnew",vxnew(i),vynew(i),vznew(i)
c     write(77,51) "aki ", sqrt(pa-pb**2/pc) ,pb* Xn(i)/ pc
c     write(77,*)                        !
c                                   !
c  51 format (a,3f10.4)             !
c     end if                        !
                                    !
      end do   ! <------------------

* ------------------
* evolve speeds of confined atoms
* ------------------

      do i= ncage+1, npar ! <-------
                                    !
      amas= atomass(i)              !
                           ! ----------------
                           ! systematic forces
                                    !
      a1x= -fvel *Gradx(i) /amas    ! ab-initio
      a1y= -fvel *Grady(i) /amas    !
      a1z= -fvel *Gradz(i) /amas    !
                                    !
      a2x= -fvel *Gradxnew(i) /amas ! ab-initio
      a2y= -fvel *Gradynew(i) /amas !
      a2z= -fvel *Gradznew(i) /amas !
                                    !
                            ! ---------------
                            ! speed variables
                                    !
* v(t)= v(t-dt) + [F(t)/m + F(t-dt)/m]/2 * dt
                                    !
      vxnew(i)= vx(i)+ (a1x + a2x)/ 2.d0*dt ! in Angst/fs
      vynew(i)= vy(i)+ (a1y + a2y)/ 2.d0*dt
      vznew(i)= vz(i)+ (a1z + a2z)/ 2.d0*dt
                                    !
      end do   ! <------------------


      return
      end
* *******************************
* make the distribution of rnd numbers
* to satisfy a gaussian function with
* variance 1 and zero mean
* *******************************
      function gauss()
      implicit double precision (a-h,o-z)

      real*8 zbqlu01

* zbqlu01 is in [0,1]

   30 v1= 2.*zbqlu01() -1.
      v2= 2.*zbqlu01() -1.

      r= v1*v1 + v2*v2

      if (r.ge.1) goto 30

      fac= sqrt (-2. *log(r)/r)
      gauss= v2*fac

      return
      end
* *****************************
* ***********************************
* average the temp over the later ifreq steps
* ***********************************
      subroutine averageT (k,temp,ifreq,averTemp)
      implicit double precision (a-h,o-z)
!      include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension temp(maxt)

      icount= 0
      sum= 0.d0

      do i= k-ifreq, k
       icount= icount +1
       sum= sum + temp(i)
      end do

      averTemp= sum/dble(icount)

      return
      end
* ***********************************
* *****************************
* break a bond in stochastic form
* more info in:
* http://www.epa.gov/radiation/understand/ionize_nonionize.html
* *****************************
      subroutine breakBond (npar,ncage, atsymbol, x,y,z)
      implicit double precision (a-h,o-z)
!      include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      CHARACTER atsymbol(nmax)*5
      dimension x(nmax), y(nmax), z(nmax)
      integer itag(nmax)

      real*8 lambda, zbqlu01

        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
         NombreOllaOut = trim(RutaArchivosSalida)
     &//'/outgrl_'//trim(NombreDinamica)
     &// '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')

* initialize labels

      icycle= 1

      do i= 1, npar
      itag(i)= 0
      end do

* ------------------------
* stochastically chose an atom
* ------------------------

   10 eta= zbqlu01()       ! zbqlu01 is in [0,1]
      nconf= npar-ncage    ! num of confined atoms
      n1= int(eta*nconf)   ! n1 is in [0,nconf]
      if (n1.eq.0) goto 10

      k= ncage + n1        ! this is the chosen atom
      if (itag(k).eq.1) goto 10

      write(77,"(a,f8.3)") "rnd number", eta
      itag(k)= 1

* ------------------------
* find the closest neighbor to atom k
* ------------------------

      small= 1000.

      do i= ncage+1, npar ! <----
                                 !
      if (i.ne.k) then           !
      v1= (x(k)-x(i))**2+ (y(k)-y(i))**2+ (z(k)-z(i))**2
      dist= sqrt(v1)             !
                                 !
      if (dist.lt.small) then    !
      small= dist                !
      ismall= i                  !
      end if                     !
      end if                     !
                                 !
      end do ! <-----------------

      j= ismall  ! this is the closest atom to k
      itag(j)= 1

      write(77,"(a,2(2x,a2,i4),a,f8.3)") "broken bond", atsymbol(k),
     &       k, atsymbol(j), j, " with distance", small

* ------------------------
* eq of the line joining atoms J and K
* http://www.tec-digital.itcr.ac.cr/revistamatematica/cur
* sos-linea/Algebra-Lineal/algebra-vectorial-geova-walter/node5.html
* ------------------------

* points to use for the line
*   P= x(j),y(j),z(j)    ;    Q= x(k),y(k),z(k)

* director vector
*   vector(PQ)= x(k)-x(j), y(k)-y(j), z(k)-z(j)

* eq of the line joining P and Q
*   x,y,z= P + lambda * vector(PQ)

* ------------------------
* open the molecular bond
* ------------------------

      dlam= 0.25 ! factor to increase the bond per atom
      lambda= 1.0 + dlam

      a= x(j) +lambda *(x(k)-x(j))
      b= y(j) +lambda *(y(k)-y(j))
      c= z(j) +lambda *(z(k)-z(j))

      d= x(k) -lambda *(x(k)-x(j))
      e= y(k) -lambda *(y(k)-y(j))
      f= z(k) -lambda *(z(k)-z(j))

      dist1= sqrt((a-d)**2 + (b-e)**2 + (c-f)**2)
      dist2= sqrt((x(k)-a)**2 + (y(k)-b)**2 + (z(k)-c)**2) ! increment for atom k
      dist3= sqrt((x(j)-d)**2 + (y(j)-e)**2 + (z(j)-f)**2) ! increment for atom j

* assign the new positions after opening the bond

      x(j)= d
      y(j)= e
      z(j)= f

      x(k)= a
      y(k)= b
      z(k)= c

* -----------------------------
* write info of the molecular bond
* -----------------------------

c     write(77,"(a,3f9.3)") atsymbol(j), x(j), y(j), z(j)
c     write(77,"(a,3f9.3)") atsymbol(k), x(k), y(k), z(k)

      write(77,"(a,f9.3)") "new distance between atoms", dist1
c     write(77,"(a,f9.3)") "increment from atom k", dist2
c     write(77,"(a,f9.3)") "increment from atom j", dist3
c     write(77,"(a,f9.3)") "total increment", dist2+dist3
      write(77,"(a,f9.3)") "increase factor", dist1/small

      icycle= icycle +1
c     if (icycle.le.1) goto 10    ! in case of breaking several bonds

      return

      end
* *****************************
* ***************************************
* get cartesian coordinates (a,b,c) from
* the polar coordinates (rho,teta,fi)
* ***************************************
      subroutine cartesian (rho,teta,fi,a,b,c)
      implicit double precision (a-h,o-z)
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      pi= 3.1415926535898
      rad= pi/180.0

*  0 <= rho
*  0 <= teta <= 2pi
*  0 <= fi <= pi

      a= rho* sin(fi)* cos(teta)
      b= rho* sin(fi)* sin(teta)
      c= rho* cos(fi)

      return
      end
* ****************************

* ***********************************
* citations of IMD, TeraChem, etc.
* ***********************************
      subroutine citation (imodule)
      implicit double precision (a-h,o-z)
!      include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      CHARACTER imodule*40


        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')
      


      write(77,*) "--------------------------------------------"
      write(77,*) "  CITATION:"
      write(77,*)
      write(77,*) "Interactive Molecular Dynamics Program"
      write(77,*) "Module: ", imodule
      write(77,*) "Ruben Santamaria"
      write(77,*) "rso@fisica.unam.mx"
      write(77,*) "IFUNAM, 2014"
      write(77,*)
c     write(77,*) "TeraChem"
c     write(77,*) "I.S. Ufimtsev and T.J. Martinez"
c     write(77,*) "Quantum Chemistry on Graphical Processing Units. 3."
c     write(77,*) "Analytical Energy Grad. & First Princ. Molec. Dyn."
c     write(77,*) "Jour. Chem. Theory & Comput. Vol. 5, p2619, 2009"
c     write(77,*)
c     write(77,*) "random-number-generator module from:"
c     write(77,*) "Richard Chandler & Paul Northrop"
c     write(77,*) "richard@stats.ucl.ac.uk, northrop@stats.ox.ac.uk"
c     write(77,*)
c     write(77,*) "I, Ruben Santamaria, have taken every precaution for"
c     write(77,*) "the writing of this program. Although I have made my"
c     write(77,*) "best effort, it is not guaranteed to be free of error."
c     write(77,*) "Therefore, it should not be relied on for solving"
c     write(77,*) "problems where an error could result in injury or loss."
c     write(77,*) "The author disclaims all liability in this regard."
c     write(77,*) "Ruben Santamaria retains the rights to make changes to"
c     write(77,*) "this program at any time, without previous notice. The"
c     write(77,*) "program can be used, copied, modified and distributed"
c     write(77,*) "under the condition of giving the acknowledgments to"
c     write(77,*) "to the author, as indicated in the citation section"
c     write(77,*) "of the program."
c     write(77,*)
      write(77,*) "wall-clock time"
      call system ("date")
      write(77,*)
      write(77,*) "                  END"
      write(77,*) "* ----------------------------------- *"
      write(77,*)


      return
      end
* ********************

* ************************************
* measure distances between particles
* ************************************
      subroutine distances (time, npar,x,y,z)
      implicit double precision (a-h,o-z)
!      include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension x(nmax), y(nmax), z(nmax)

* ------------------------------------
* find the smallest & largest distances
* ------------------------------------

      small= 1.d+12
      big=  -1.d+12

      do i= 1, npar ! <----------
      do j= 1, npar              !
                                 !
      if (i.lt.j) then ! <-------!------
                                 !
      v1= (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
      Rij= sqrt (v1)             !
                                 !
      if (Rij.lt.small) then     ! get the smallest
      small= Rij                 ! distance
      ismall= i                  !
      jsmall= j                  !
      end if                     !
                                 !
      if (Rij.gt.big) then       ! get the largest
      big= Rij                   ! distance
      ibig= i                    !
      jbig= j                    !
      end if                     !
                                 !
      end if ! <-----------------!------
                                 !
      end do                     !
      end do !<------------------

c     write(77,20) "small", time, ismall, jsmall, small
c     write(77,20) "big  ", time, ibig,   jbig,   big

   20 format (1x,a,f8.3,2i5,f8.3)

      return
      end
* ***********************************

* ***********************************
* get ab-initio energy & gradients
* ***********************************
      subroutine eneNwchem (time,npar, Grdx,Grdy,Grdz, Epot,refPE)
      implicit double precision (a-h,o-z)
!  include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      CHARACTER qline*(maxchr)
      dimension itok(maxtok), jtok(maxtok)

      dimension xcrds(10)
      dimension Grdx(nmax), Grdy(nmax), Grdz(nmax)

      external kparse
 
        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut
        CHARACTER (LEN = 400) NwchemOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'
        
        NwchemOut = trim(RutaArchivosSalida)
     & // '/NwchemOut_' // trim(NombreDinamica)
     & // '.out'
      
c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')

c     open (unit=16, file="nwchem.out", status='unknown')
      open (unit=16, file=trim(NwchemOut), status='unknown')

* initialize variables

      abene= 0.d0
      Epot=  0.d0

      do i= 1, 10
      xcrds(i)= 0.d0
      end do

      do k= 1, npar
       Grdx(k)= 0.d0
       Grdy(k)= 0.d0
       Grdz(k)= 0.d0
      end do

* ------------------------
* find the section of energy
* ------------------------

   10 read(16,20, end=50) qline ! <------
                                         !
      ntok= kparse(qline,maxchr,itok,jtok,maxtok)
                                         !
      if ((qline(itok(1):jtok(1)) .eq. 'Total') .and.
     &    (qline(itok(2):jtok(2)) .eq. 'DFT')   .and.
     &    (qline(itok(3):jtok(3)) .eq. 'energy')) then
                                         !
                          ! ---------------------------
                          ! ab-initio energy [Hartree]
                                         !
      read(qline(itok(5):jtok(5)),*) abene
                                         !
                                         ! take the energy of the first
      if (time.lt.0.d0) refPE= abene     ! geometry as the ref pot ene
                                         !
      Epot= (abene-refPE)                ! in Hartree
      goto 13                            !
      end if                             !
                                         !
      goto 10 ! <------------------------

   13 continue

* ----------------------------------
* find the section of energy gradients
* ----------------------------------

   30 read(16,20, end=50) qline ! <------
                                         !
      ntok= kparse(qline,maxchr,itok,jtok,maxtok)
                                         !
      if ((qline(itok(1):jtok(1)) .eq. 'DFT') .and.
     &    (qline(itok(2):jtok(2)) .eq. 'ENERGY') .and.
     &    (qline(itok(3):jtok(3)) .eq. 'GRADIENTS')) then
                                         !
      do k= 1, 3                         ! junk files
      read(16,20) qline                  !
      end do                             !
                                         !
      do k= 1, npar               ! <----!--
      read(16,20) qline                  !
      ntok= kparse(qline,maxchr,itok,jtok,maxtok)
                                         !
      do j= 6, 8                         !
      read(qline(itok(j):jtok(j)),*) xcrds(j)
      end do                             !
                                         !
      Grdx(k)= xcrds(6)                  ! in Hartree/Bohr
      Grdy(k)= xcrds(7)                  !
      Grdz(k)= xcrds(8)                  !
      end do                      ! <----!--
                                         !
      goto 50                            !
      end if                             !
                                         !
      goto 30 ! <------------------------

   50 continue

      close (16)
      close (91)

c     do i= 1, npar
c      Grdx(i)= -Grdx(i)
c      Grdy(i)= -Grdy(i)
c     end do

* ---------------------------------
* obtain the largest forces in x,y,z
* ---------------------------------

      big=  -100.d0
      xbig= -100.d0
      ybig= -100.d0
      zbig= -100.d0

      do k= 1, npar ! <---------
                                !
       fx= abs (Grdx(k))        !
       fy= abs (Grdy(k))        !
       fz= abs (Grdz(k))        !
                                !
       ff= sqrt (fx**2 + fy**2 + fz**2)
                                !
       if (fx.gt.xbig) then     ! x-component
       ixbig= k                 !
       xbig= fx                 !
       end if                   !
                                !
       if (fy.gt.ybig) then     ! y-component
       iybig= k                 !
       ybig= fy                 !
       end if                   !
                                !
       if (fz.gt.zbig) then     ! z-component
       izbig= k                 !
       zbig= fz                 !
       end if                   !
                                !
       if (ff.gt.big) then      ! magnitude
       ibig= k                  !
       big= ff                  !
       end if                   !
                                !
      end do ! <----------------

      write(77,*)  "largest (absolute value) forces in x,y,z"
      write(77,12) "for atoms:", ixbig, iybig, izbig
      write(77,14) "forces [Hartree/Bohr]", xbig, ybig, zbig
      write(77,12) "largest force magnitude for atom", ibig
      write(77,14) "force magnitude [Hartree/Bohr]", big
      write(77,*)

   12 format (1x,a,3(1x,i4))
   14 format (1x,a,3(1x,f11.5))

* ------------------------------
* broadcast problems of terachem
* ------------------------------

      if (abene.eq.0.0) then
        write(77,*) "NWCHEM: got NO energy"
        write(77,*) "possible lack of energy convergence"
        write(77,*)
c       stop
      end if

   20 format (a)


      return
      end
* *********************************
* ***********************************
* get the potential energy of the springs
* ***********************************
      subroutine eneSpring (time,npar,ncage,x,y,z,x0,y0,z0,sK,Espring)
      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension x0(nmax), y0(nmax), z0(nmax)
      dimension x(nmax),  y(nmax),  z(nmax)
      dimension sK(nmax)

        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')
* ------------------------
* potential energy of springs
* ------------------------

      Espring= 0.d0
      sum= 0.d0

      do i= 1, ncage
      v1= (x(i)-x0(i))**2 +(y(i)-y0(i))**2 +(z(i)-z0(i))**2
      sum= sum +v1
      end do

* all spring constants are the same

      if (time.lt.0.d0) then
      write(77,*) "all spring constants are considered the same"
      write(77,*)
      end if

      Espring= sK(2)*sum/2.  ! in (N/cm)*Angst^2

* Angst^2/cm= 1x10^{-20} m^2/ 1x10^{-02} m= 1x10^{-18} m

* Espring= Espring* 1.0d-18  ! in N*m=Joules

* 1 Hartree= 4.35988 x10^{-18} J

      Espring= Espring/ 4.35988   ! in Hartree
      return
      end
* *********************************
* ***********************************
* estimate CPU time from the first computations
* ***********************************
      subroutine estimateCPU (value)
      implicit double precision (a-h,o-z)
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      CHARACTER qline*(maxchr)
      dimension itok(maxtok), jtok(maxtok)

      external kparse

      CHARACTER (LEN = 400) NombreDinamica
      CHARACTER (LEN = 400) RutaArchivosEntrada
      CHARACTER (LEN = 400) RutaArchivosSalida
      CHARACTER (LEN = 400) TerachemOut

      call DatosGrlsDinamica(NombreDinamica,
     & RutaArchivosEntrada, RutaArchivosSalida)


      TerachemOut = trim(RutaArchivosSalida) 
     & // '/TerachemOut_' // trim(NombreDinamica)
     & // '.out'


c     open (unit=16, file="terachem.out", status='unknown')
      open (unit=16, file=trim(TerachemOut), status='unknown')

* --------------------------
* find the section of interest
* --------------------------

   10 read(16,20, end=50) qline ! <------
                                         !
      ntok= kparse(qline,maxchr,itok,jtok,maxtok)
                                         !
      if ((qline(itok(1):jtok(1)) .eq. 'Total') .and.
     &    (qline(itok(2):jtok(2)) .eq. 'processing') .and.
     &    (qline(itok(3):jtok(3)) .eq. 'time:')) then
                                         !
c     write(77,*) "qline= ", qline           !
                                         ! read the time
      read(qline(itok(4):jtok(4)),*) value
                                         !
c     write(77,*) "valor encontrado", value  !
                                         !
      end if                             !
      goto 10 ! <------------------------

   50 continue
   20 format (a)

      close (16)

      return
      end
* *********************************
* *****************************
* gather info of the user, machine, etc
* *****************************
      subroutine infoUser

      implicit double precision (a-h,o-z)
             
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)

      CHARACTER pwd*60

        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'
        
c      open (unit = 77, file = trim(NombreOllaOut), status = 'old')
      
      
* ------------------------------------------
* identify initial cpu processes & directories
* ------------------------------------------

      write(77,*)
      call system ("whoami")
      call system ("hostname")

      call system ("pwd > t1")
      call workDir (pwd)
      write(77,"(3a)") "working directory ", pwd

c     write(77,*) "wall-clock time"
      call system ("date")
c     call system ("date -R | cut -c1-25")

      write(77,*) "fortran process ID   ", getpid()
      write(77,*) "ab-initio process ID ", getpid()+3
      write(77,*) "user ID  ", getuid()
      write(77,*) "group ID ", getgid()

      call system ("rm t1 t2")
c     call system ("rm -rf permanent scratch")
c     call system ("mkdir permanent")
c     call system ("mkdir scratch")
      return
      end
* *****************************

* ***********************************
* impose initial conditions on positions and speeds
* ***********************************
      subroutine initConditOlla (npar,ncage, atomass, x,y,z, vx,vy,vz,
     &    x0,y0,z0, cmx,cmy,cmz, tempEquil,
     &    nsteps,ibTime1,ibTime2,ibTime3)

      implicit double precision (a-h,o-z)
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      external gauss

      dimension atomass(nmax)
      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)

      dimension x0(nmax), y0(nmax), z0(nmax), sK(nmax)

      integer iBreakTime(3)
      real*8 zbqlu01, kB, iniTemp


        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')

* ------------------
* define constants
* ------------------

* Boltzmann constant

      kB= 1.3806504d-23  ! in Joule/Kelvin= Kg*(m/s)^2/Kelvin

* facKT *(Angst/fs)^2= (Joule/Kelvin)*(Kelvin/amu)

      facKT= 1.0d17/1.6605402

      iniTemp= tempEquil

* ------------------
* loop over the particles
* ------------------

      write(77,*) "moving the origin to the CM of the cage"
      write(77,*)

      write(77,12) "required temp for the simulation", iniTemp
      write(77,*) "giving maxwellian velocities to atoms"

      do i= 1, npar ! <----
                           !
                 ! -----------------
                 ! initial positions
                           !
         x(i)= x(i)-cmx    ! move the origin
         y(i)= y(i)-cmy    ! to the cage CM
         z(i)= z(i)-cmz    !
                           !
         x0(i)= x(i)       ! store initial
         y0(i)= y(i)       ! positions
         z0(i)= z(i)       !
                           !
                   ! --------------
                   ! initial speeds
                           !
      pKT= kB*iniTemp/ atomass(i)
                           !
      etax= gauss()        ! random numbers
      etay= gauss()        !
      etaz= gauss()        ! Maxwellian speeds
                           !
      vx(i)= etax *sqrt(pKT*facKT) ! in Angst/fs
      vy(i)= etay *sqrt(pKT*facKT)
      vz(i)= etaz *sqrt(pKT*facKT)
                           !
      end do ! <-----------

* ----------------------------
* get the speed of the mass center
* ----------------------------

      sumx= 0.d0
      sumy= 0.d0
      sumz= 0.d0
      sumass= 0.d0

      do i= ncage+1, npar   ! for the confined atoms
       sumass= sumass + atomass(i)
       sumx= sumx + atomass(i) *vx(i)
       sumy= sumy + atomass(i) *vy(i)
       sumz= sumz + atomass(i) *vz(i)
      end do

      v1= sumx/ sumass
      v2= sumy/ sumass
      v3= sumz/ sumass

* ---------------------------------
* refer the speeds to the speed of the CM
* ---------------------------------

      do i= ncage+1, npar ! <--
                               !
        vx(i)= vx(i)-v1        ! for the
        vy(i)= vy(i)-v2        ! confined
        vz(i)= vz(i)-v3        ! atoms
                               !
      end do ! <---------------

      write(77,*) "making zero the speed of the CM of confined atoms"

* --------------------
* kinetic energy of confined atoms
* --------------------

      Ekin= 0.d0
      nn= 0    ! momentum degrees of freedom

      do i= ncage+1, npar ! <-------
                                    ! in amu*(Angst/fs)^2
       nn= nn + 1                   !
       v2= vx(i)**2 + vy(i)**2 + vz(i)**2
       Ekin= Ekin + atomass(i)*v2/ 2.d0
                                    !
      end do   ! <------------------

* facEkin * Joule= amu *(Angst/fs)^2

      facEkin= 1.6605402d-17

      EkinConf= Ekin*facEkin    ! in Joules

* compute the temperature from:  Ekin=3NkT/2 ==> T=2Ekin/(3Nk)

      tempConf= 2.*EkinConf/ (3.*dble(nn)*kB)  ! in Kelvin

* --------------------
* kinetic energy of cage atoms
* --------------------

      Ekin= 0.d0
      nn= 0    ! momentum degrees of freedom

      do i= 1, ncage ! <---------
                                 ! in amu*(Angst/fs)^2
       nn= nn + 1                !
       v2= vx(i)**2 + vy(i)**2 + vz(i)**2
       Ekin= Ekin + atomass(i)*v2/ 2.d0
                                  !
      end do   ! <----------------

      EkinCage= Ekin*facEkin    ! in Joules

      tempCage= 2.*EkinCage/ (3.*dble(nn)*kB)  ! in Kelvin

      write(77,12) "initial temp of cage atoms [K]", tempCage
      write(77,12) "initial temp of confined atoms [K]", tempConf

   12 format (1x,a,f11.3)

* ---------------------------
* chose the times to break bonds in stochastic form
* ---------------------------

      do k= 1, 3 ! <---------
                             !
      iBreakTime(k)= 0       ! default value
                             !
   10 eta= zbqlu01()         ! zbqlu01 is in [0,1]
      n1= int(eta*nsteps)    ! n1 is in [0,nsteps]
                             !
      if ((n1.lt.100).or.(n1.gt.nsteps-100)) goto 10
      iBreakTime(k)= n1      ! this break time is taken
                             !
      end do ! <-------------

      n1= iBreakTime(1)
      n2= iBreakTime(2)
      n3= iBreakTime(3)

* find the order of the numbers n1,n2,n3

      if ((n1.lt.n2).and.(n1.lt.n3)) then
      ibTime1= n1
         if (n2.lt.n3) then
         ibTime2= n2
         ibTime3= n3   ! the order is n1<n2<n3
         else
         ibTime2= n3
         ibTime3= n2   ! the order is n1<n3<n2
         end if
      end if

      if ((n2.lt.n1).and.(n2.lt.n3)) then
      ibTime1= n2
         if (n1.lt.n3) then
         ibTime2= n1
         ibTime3= n3   ! the order is n2<n1<n3
         else
         ibTime2= n3   
         ibTime3= n1   ! the order is n2<n3<n1
         end if
      end if

      if ((n3.lt.n1).and.(n3.lt.n2)) then
      ibTime1= n3
         if (n1.lt.n2) then
         ibTime2= n1
         ibTime3= n2   ! the order is n3<n1<n2
         else
         ibTime2= n2
         ibTime3= n1   ! the order is n3<n2<n1
         end if
      end if

* a minimum difference of 10 is required between
* consecutive numbers to give chance to the charge uptake

      if (ibTime2-ibTime1.le.10) then
      ibTime2= ibTime2+13   ! increase the value of ibTime2
      end if

      if (ibTime3-ibTime2.le.10) then
      ibTime3= ibTime3+13   ! increase the value of ibTime3
      end if

      write(77,"(/,1x,a,3i5)") "rnd times to break bonds",
     &   ibTime1, ibTime2, ibTime3


      return

      end
* *******************************
* ***********************************
* create the nwchem input file
* ***********************************
      subroutine inptNwchem (ilabel,npar,atsymbol,
     & x,y,z)
     
      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      integer getcwd, status
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      CHARACTER atsymbol(nmax)*5
      dimension x(nmax), y(nmax), z(nmax)

        
        CHARACTER (LEN = 400) RutaScratch
        CHARACTER (LEN = 400) RutaPermanent
        CHARACTER (LEN = 400) BaseGrlNwchem

        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        CHARACTER (LEN = 400) NwchemInp

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')

        NwchemInp = trim(RutaArchivosEntrada)
     & // '/NwchemInp_' // trim(NombreDinamica)
     & //'.inp'
        
        call RutinaDatosNwchem (RutaScratch,
     & RutaPermanent, BaseGrlNwchem)


* ilabel= 1    generate the wavefunction (WF) from scratch
* ilabel= 2    reuse the WF from a previous computation
* ilabel= 3    print the nwchem input to the output file
* ilabel= 4    change the charge of the molecule & generate the WF
* ilabel= 5    change the charge of the molecule & reuse the WF

* -------------------------
* chose the appropriate run
* -------------------------

      if ((ilabel.eq.1).or.(ilabel.eq.4)) then
c       open (14,file="nwchem.inp",status="unknown")
        open (14,file=trim(NwchemInp),status="unknown")
        write(77,"(/,a,i5)") "generating the wave function from 
     & scratch", ilabel
        ij= 14
      end if

      if ((ilabel.eq.2).or.(ilabel.eq.5)) then
c       open (14,file="nwchem.inp",status="unknown")
        open (14,file=trim(NwchemInp),status="unknown")
        write(77,"(/,a,i5)") "reusing the wave function", ilabel
        ij= 14
      end if

      if (ilabel.eq.3) then
        write(77,*) "the ab-initio input file may be"
        write(77,*) "changed in module inptNwchem"
        ij= 6
      end if

* ---------------------
* write basic parameters
* ---------------------
      

      write(ij,*)
      write(ij,*) "title ", char(34), "mol computation", char(34)
      write(ij,*) "start molec"
c     write(ij,*) "echo"

      write(ij,*)

      write(ij,*) "memory total 4000 stack 1000 heap 1000 global 2000
     & mb"
      write(ij,*) trim(RutaPermanent)
      write(ij,*) trim(RutaScratch)

      write(ij,*)

      if ((ilabel.eq.4).or.(ilabel.eq.5)) then
      if (ilabel.eq.4) write(77,*) "charge changed in the input file"
       write(ij,*) "charge 0"   ! change charge
      else
       write(ij,*) "charge 0"  ! starting charge
      end if

* --------------
* write mol coords
* --------------

      write(ij,*)

      write(ij,*) "geometry units angstroms noautoz noautosym nocenter"

      if (ilabel.ne.3) then
      do i= 1, npar
      write(ij,"(2x,a,3f9.3)") atsymbol(i), x(i),y(i),z(i)
      end do
      end if

      write(ij,*) "end"

* --------------
* write basis sets
* --------------

      write(ij,*)
      write(ij,*) "basis"
      write(ij,*) "  * library "//trim(BaseGrlNwchem)
      write(ij,*) "end"

c    write(ij,*)
c     write(ij,*) "basis"
c     write(ij,*) "  H library 6-31g"
c     write(ij,*) "  C library 6-31g"
c     write(ij,*) "  N library 6-31g"
c     write(ij,*) "  O library 6-31g"
c     write(ij,*) "  S library 6-31g"
c     write(ij,*) " He library 6-31g"
c     write(ij,*) " Pd library ", char(34), "LANL2DZ ECP", char(34)
c     write(ij,*) "end"

c     write(ij,*)
c     write(ij,*) "set basis ", char(34), "cd basis", char(34)
c     write(ij,*) "basis ", char(34), "cd basis", char(34)
c     write(ij,*) "  H library 6-31g"
c     write(ij,*) "  C library 6-31g"
c     write(ij,*) "  N library 6-31g"
c     write(ij,*) "  O library 6-31g"
c     write(ij,*) "  S library 6-31g"
c     write(ij,*) " He library 6-31g"
c     write(ij,*) " Pd library ", char(34), "LANL2DZ ECP", char(34)
c     write(ij,*) "end"

      write(ij,*)
      write(ij,*) "ecp"
      write(ij,*) " Pd library ", char(34), "LANL2DZ ECP", char(34)
      write(ij,*) "end"

* --------------------
* write level of theory
* --------------------

      write(ij,*)

      write(ij,*) "dft"
      write(ij,*) " xc becke88 lyp"

      if ((ilabel.eq.4).or.(ilabel.eq.5)) then
      if (ilabel.eq.4) write(77,*) "multiplicity changed in the input"
     & //"file"
       write(ij,*) " mult 1"    ! change multiplicity
      else
       write(ij,*) " mult 1"    ! starting multiplicity
      end if

      write(ij,*) " iterations 1500"
      write(ij,*) " convergence energy  1e-5"
      write(ij,*) " convergence density 1e-5"
      write(ij,*) " grid medium"
      write(ij,*) " mulliken"

      if ((ilabel.eq.2).or.(ilabel.eq.5)) then
      write(ij,*)
     &  " vectors input "//RutaPermanent//"/molec.movecs"
      end if

      write(ij,*) "end"
      write(ij,*)

* ------------------
* write type of task
* ------------------

      write(ij,*) "task dft gradient"
      write(ij,*)

      if (ij.eq.14) close (ij)
      return
      end
* ********************
* **********************************
* compute kinetic energies & temperature
* **********************************
      subroutine kinEne (npar,ncage, atomass, vx,vy,vz,
     &     EkinConf,EkinCage, tempConf,tempCage,Vol,dynPress)

      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension atomass(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)

      real*8 kB, mdof

* Boltzmann constant

      kB= 1.3806504d-23  ! in Joule/Kelvin= Kg*(m/s)^2/Kelvin

* --------------------
* kinetic energy of confined atoms
* --------------------

      Ekin= 0.d0
      mdof= 0.   ! momentum degrees of freedom

      do i= ncage+1, npar ! <-------
                                    ! in amu*(Angst/fs)^2
       mdof= mdof + 1.              !
       v2= vx(i)**2 + vy(i)**2 + vz(i)**2
       Ekin= Ekin + atomass(i)*v2/ 2.d0
                                    !
      end do   ! <------------------

* facEkin * Joule= amu *(Angst/fs)^2

      facEkin= 1.6605402d-17

      EkinConf= Ekin*facEkin    ! in Joules

* compute the temperature from:  Ekin=3NkT/2 ==> T=2Ekin/(3Nk)

      tempConf= 2.*EkinConf/ (3.*mdof*kB)  ! in Kelvin

* --------------------
* kinetic energy of cage atoms
* --------------------

      Ekin= 0.d0
      mdof= 0.   ! momentum degrees of freedom

      do i= 1, ncage ! <---------
                                 ! in amu*(Angst/fs)^2
       mdof= mdof + 1.           !
       v2= vx(i)**2 + vy(i)**2 + vz(i)**2
       Ekin= Ekin + atomass(i)*v2/ 2.d0
                                  !
      end do   ! <----------------

      EkinCage= Ekin*facEkin    ! in Joules

      tempCage= 2.*EkinCage/ (3.*mdof*kB)  ! in Kelvin

* Francis W. Sears \& Gerhard L. Salinger,
* {\it Termodin\'amica, teor\'ia cin\'etica, y termodin\'amica estad\'istica},
* 2nd ed., Revert\'e S. A. M\'exico, 1978, page 300.
* http://www.seykota.com/rm/gas/gas.htm

* --------------------
* dynamic pressure from the ideal gas:
* PV= NkT= Nk [2Ekin/(3Nk)]= 2Ekin/3
* --------------------

* from above T=2Ekin/(3Nk)  ==> PV=NkT=Nk(2Ekin/(3Nk))=2Ekin/3

* facPres2 * Pa= Joule/Angst^3

      facPres2= 1.0d30

      dynPress= 2./3.* EkinConf/Vol *facPres2 ! in Pascals

* 1 GPa= 1x10^{9} Pa
c     dynPress= dynPress/1.d09  ! GPa

* 1 pascal= (1/101325.0) atmosphere

      dynPress= dynPress/ 101325.0  ! atmospheres

* 1 Hartree= 4.35988 x 10^{-18} J

      EkinConf= EkinConf/4.35988d-18   ! in Hartree
      EkinCage= EkinCage/4.35988d-18   ! in Hartree

      return
      end
* *********************************
* ***********************************
* this part illustrates the call
* to the function kparse
* ***********************************
c     subroutine test
c     implicit double precision (a-h,o-z)
c     parameter (maxchr=128, maxtok=50, m1=200)

c     CHARACTER qline*(maxchr)
c     dimension itok(maxtok),jtok(maxtok)
c     external kparse
c     external lenstr

* start debugging lines

c  10 read(11,900, end=50) qline
c     ntok= kparse(qline,maxchr,itok,jtok,maxtok)

c     if ((qline(itok(1):jtok(1)) .eq. 'DFT') .and.       ! 11111
c    &    (qline(itok(2):jtok(2)) .eq. 'ENERGY') .and.
c    &    (qline(itok(3):jtok(3)) .eq. 'GRADIENTS')) then

c     palabra= qline(itok(3):jtok(3))
c     nn= lenstr (palabra, maxchr)
c     write(77,*) "longitud de la palabra ", palabra, "es", nn
c     write(77,*) "longitud de la palabra ", jtok(3)-itok(3)+1

c  50 write(77,*) ' '
c 900 format(a)

* ************************************
* parte en palabras a la linea "qline"
* regresando los indices itok y jtok.
* ************************************
      function kparse(qline,lenq,itok,jtok,maxtok)
      implicit double precision (a-h,o-z)
      dimension itok(*),jtok(*)
      CHARACTER qline*(*)
      CHARACTER qspace*1,qtab*1
 
* qline(itok(1):jtok(1) <-- 1st word
* itok(1)   <-- first letter of 1st word
* jtok(1)   <-- last letter of 1st word

* qline(itok(2):jtok(2) <-- 2nd word
* itok(2)   <-- first letter of 2nd word
* jtok(2)   <-- last letter of 2nd word

      kparse= 0
      kount= 0
 
      qspace= ' '
      qtab= char(9) ! refers to a tab
 
* --------------------------------------------
* look for the first (non-blank or non-tab) token
* --------------------------------------------

   10 kount= kount + 1 ! <---------------------------
                                                     !
      if (kount.gt.lenq) goto 40                     !
                                                     !
* we are not interested on spaces or tabs            !
                                                     !
      if (qline(kount:kount) .eq. qspace .or.        !
     &    qline(kount:kount) .eq. qtab) goto 10 ! <--
 
* we've found a new token, check for maximum token count
 
      if (kparse .eq. maxtok) goto 40
 
* mark the beginning
 
      kparse= kparse + 1  ! label for the word
      itok(kparse)= kount ! first letter of the word

* -------------------------------
* look for the end of this token
* -------------------------------
 
   20 kount= kount + 1 ! <----------------------------
                                                      !
      if (kount .gt. lenq) goto 30                    !
                                                      !
      if (qline(kount:kount) .ne. qspace .and.        !
     &    qline(kount:kount) .ne. qtab) goto 20 ! <---
 
* end of the token, mark it
 
   30 jtok(kparse)= kount-1 ! last letter of the word

* ------------------------------- 
* continue searching for words
* ------------------------------- 
 
      if (kount.lt.lenq) goto 10
 
   40 return
      end
* ********************************
* computes the length of a string
* ********************************
      function lenq(qstring,max)
      implicit double precision (a-h,o-z)
      CHARACTER qstring*(*)
 
      lenq= max

* find a blank space

   10 if (qstring(lenq:lenq).ne.' ') return

      lenq = lenq - 1
      if (lenq.le.0) return

      goto 10
 
      end
* ********************************
* computes the length of a string
* ********************************
      integer function lenstr(string, max)
      implicit none
      CHARACTER*(*) string
      integer max

      lenstr= max

* find a blank space

   10 if (string(lenstr:lenstr).ne.' ') return
      write(77,*)"aqui stoy ", string, lenstr
      write(77,*)"aqui stoy 123456789012345678901234567890"

      lenstr= lenstr - 1

* if the lenght is 0 return

      if (lenstr.le.0) return

      goto 10

      return
      end
* ********************************
* ************************************
* get the effective radius of the confined particles
* ************************************
      subroutine largest (t,npar,ncage, x,y,z, radConf,radCage, Vol)
      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension x(nmax), y(nmax), z(nmax)
      dimension dist(nmax)

      integer itag(nmax), jtag(nmax)

        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')

      pi= 3.1415926535898

* -----------------
* initialize tags
* -----------------

      knt= 0

      do j= ncage+1, npar ! <----
       knt= knt + 1              !
                                 !
       itag(knt)= 0              !
       jtag(knt)= 0              ! radial distances
                                 !
       dist(knt)= sqrt(x(j)**2 +y(j)**2 +z(j)**2)
                                 !
      end do ! <-----------------

      idim= knt

* ---------------------------------
* find the largest numbers of the array
* ---------------------------------

      do i= 1, idim ! <---------
                                !
      max= 0                    !
      dmax= 0.d0                !
                                !
* get the largest distance without a label
                                !
      do j= 1, idim             !
        if ((dist(j).ge.dmax).and.(itag(j).eq.0)) then
        dmax= dist(j)           !
        max= j                  !
        end if                  !
      end do                    !
                                !
      itag(max)= 1 ! switch to avoid examination of this num
      jtag(i)= max ! store the index of the largest number
                                !
      end do ! <----------------

* jtag(1)    ! <---- index of the largest num
* jtag(2)
* ......
* jtag(idim) ! <---- index of the smallest num

* ---------------------------
* radius of the confined particles
* ---------------------------

* 30% of atoms determine the confinement radius

c     nn= (npar-ncage)/10 *3
      nn= 3

* sum radial distances of most external atoms

      knt= 0
      sum= 0.d0

      do i= 1, nn ! <-----
                          !
       knt= knt + 1       !
       sum= sum + dist(jtag(i))
                          !
c      write(77,*) i, jtag(i), dist(jtag(i))
                          !
      end do ! <----------

      radPart= sum /dble(knt)    ! average radius

      if (t.eq.-1.0) then
      write(77,14) "num of atoms defining the confinement radius",knt
      end if

* ------------------------
* get the effective confinement radius
* ------------------------

* radius of the cage

      sum= 0.d0

      do i= 1, ncage
      sum= sum +sqrt (x(i)**2 +y(i)**2 +z(i)**2)
      end do

      radCage= sum /dble(ncage) ! cage radius

* effective confinement radius

      radConf= radPart +(radCage-radPart)/ 2.

      area=   4.*pi* radConf**2    ! in Angst^2
      Vol= 4./3.*pi *radConf**3    ! in Angst^3

   14 format (1x,a,1x,i5)

      return
      end
* ************************************
* ***********************************
* make a movie from the atom coords
* ***********************************
      subroutine makeMovie (iframe,npar,atsymbol,x,y,z)
      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      CHARACTER atsymbol(nmax)*5
      dimension x(nmax), y(nmax), z(nmax)
       
  
* -------------------------
* write coords for a movie
* -------------------------

cc     write(12,*) npar + 1
c     write(12,12) npar
      write(92,12) npar
c     write(12,13) "frame", iframe
      write(92,13) "frame", iframe
      

      do i= 1, npar
c     write(12,15) atsymbol(i), x(i),y(i),z(i)
      write(92,15) atsymbol(i), x(i),y(i),z(i)
      end do

* writing an additional atom

c     write(12,16) iframe, 3.0, 0.000, 0.000

   12 format (1x,i5)
   13 format (1x,a5,2x,i5)
   15 format (1x,a5,3(2x,f10.3))
   16 format (1x,i5,"X", 3(1x,f8.4))

      return
      end
* **************************************
* ******************************************
* get the center of mass from atom positions
* ******************************************
      subroutine massCenter (nk1,nk2,atmass,x,y,z,r1,r2,r3)
      implicit double precision (a-h,o-z) 
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension atmass(nmax)
      dimension x(nmax),y(nmax),z(nmax)

        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')

* initialize summations

      sumx= 0.d0
      sumy= 0.d0
      sumz= 0.d0
      sumass= 0.d0

* get the mass center of the first nk1 particles

      do i= 1, nk1
       sumass= sumass + atmass(i)
       sumx= sumx + atmass(i) *x(i)
       sumy= sumy + atmass(i) *y(i)
       sumz= sumz + atmass(i) *z(i)
      end do

      r1= sumx/ sumass
      r2= sumy/ sumass
      r3= sumz/ sumass

* print the mass center position

      write(77,20) r1,r2,r3

   20 format(1x,"position of the CM",3(1x,f8.4))

      return
      end
* *************************************************
* **************************************
* collect information for the dynamics
* **************************************
      subroutine molInfo (npar,ncage, atsymbol,atomass, x,y,z,
     &   dt,Time, rfac,vfac, tempEquil, sK,xi,
     &   radiusf, ifreqShr,ifreqStor,irestDyn,iwavFunct,nsteps,
     &itest)

      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      CHARACTER atsymbol(nmax)*5, atname*12

      dimension atomass(nmax)
      dimension x(nmax), y(nmax), z(nmax)
      dimension sK(nmax)

      CHARACTER name, titleName*100
      integer atnum

        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut
        CHARACTER (LEN = 400) NombreOllaInp

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)

        NombreOllaInp = trim(RutaArchivosEntrada) //
     & '/olla_'//trim(NombreDinamica)//'.inp'

        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'
      
        print*, "NOMBREOLLAINP", NombreOllaInp
        pi= 3.1415926535898
      speedC= 299792458  ! m/s


* open the file with parameters for the dynamics

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')
        open (unit = 13, file = trim(NombreOllaInp))
* ------------------------------
* read parameters of the dynamics
* ------------------------------

      read(13,*) titleName       ! title of the dynamics
      read(13,*) name, dt        ! time step of the dynamics
      read(13,*) name, Time      ! total time of the dynamics
      read(13,*) name, rfac      ! scale factor of atom positions
      read(13,*) name, vfac      ! scale factor of atom velocities
      read(13,*) name, npar      ! num of particles
      read(13,*) name, ncage     ! number of cage atoms  [1,ncage]
      read(13,*) name, tempEquil ! required equilibrium temperature
      read(13,*) name, springK   ! spring constant K
      read(13,*) name, xi        ! friction constant of the solvent
      read(13,*) name, radiusf   ! factor to shrink the cage
      read(13,*) name, ifreqShr  ! freq to shrink the cage
      read(13,*) name, ifreqStor ! freq to store data in case of a crash
      read(13,*) name, irestDyn  ! restart the dynamics
      read(13,*) name, iwavFunct ! generate or reuse the wavefunction
      read(13,*) name, itest     ! ful prueba o dft (1 o 2)
      read(13,*) name

      write(77,*) "-----------------------------"
      write(77,*) titleName
      write(77,*) "initial parameters of the dynamics"
      write(77,16) "time step [fs]  ", dt
      write(77,16) "total simulation time [fs]", Time
      write(77,16) "rfac      ", rfac
      write(77,16) "vfac      ", vfac
      write(77,14) "number of atoms", npar
      write(77,14) "ncage     ",      ncage
      write(77,16) "tempEquil [K]",   tempEquil
      write(77,16) "spring constant [N/cm]",    springK
      write(77,16) "friction constant [ps^-1]", xi
      write(77,16) "factor to shrink the cage", radiusf
      write(77,14) "freq to shrink the cage",   ifreqShr
      write(77,14) "freq to store data",        ifreqStor

      write(77,*)
      if (irestDyn.eq.2) write(77,14) "restart the dynamics   yes"
      if (irestDyn.eq.1) write(77,14) "restart the dynamics   no"

      if (iwavFunct.eq.1) write(77,14) "wave function from scratch"
      if (iwavFunct.eq.2) write(77,14) "reusing the wave function"

      if ((iwavFunct.lt.1).or.(iwavFunct.gt.2)) then
        write(77,*)
        write(77,*) "BAD value for starting the wavefunction"
        write(77,*)
        stop
      end if

      write(77,*) "-----------------------------"
      write(77,*) "reading:"
      write(77,*) "coords [Angst], speeds [Angst/fs] & masses [amu]"

c     write(77,*)
c     write(77,*) "num symbol  atnum  atname        atmass        x
c    &y         z"

* ------------------------------
* identify atoms & their features
* ------------------------------

      isum= 0

      do i= 1, npar ! <------
                             !
      isum= isum +1          ! define defaults
                             !
      atsymbol(i)= "NADA"    ! must NOT be an atom symbol
      atnum= 500000          ! must NOT be an atom number
      atname= "NOSE"         ! must NOT be an atom name
      atmass= 500000.0d0     ! must NOT be an atom mass
                             !
      read(13,*) atsymbol(i), x(i), y(i), z(i)
                             !
      call periodicTable (i,atsymbol,atnum,atname,atmass)
                             !
      atomass(i)= atmass     !
                             !
c     write(77,13) i, atsymbol(i), atnum, atname, atmass,
c    &    x(i), y(i), z(i)   !
                             !
      end do ! <-------------

      close (13)   ! file olla.inp

      write(77,14) "number of cage atoms", ncage
      write(77,14) "number of confined atoms", npar-ncage
      write(77,14) "total num of atoms", isum

* -------------------------------
* get initial radius of the cage
* -------------------------------

      big= -1.d0
      ibig= 0
      sum= 0.d0

      do i= 1, ncage ! <-----
                             !
        dist= sqrt (x(i)**2 + y(i)**2 + z(i)**2)
        sum= sum + dist      !
                             !
       if (dist.gt.big) then ! get the largest
        big= dist            ! distance deviation
        ibig= i              !
       end if                !
                             ! all spring constants
        sK(i)= springK       ! are the same
                             !
      end do ! <-------------

      cradius= sum/ dble(ncage) ! central radius

* reduction of the cage after a num of steps

      xnum= int(Time/ real(ifreqShr))
      rfin= cradius *radiusf**xnum

      write(77,*)
      write(77,16) "average cage radius [Angst]", cradius
      write(77,16) "largest distance deviation", big
      write(77,14) "of cage atom", ibig

      if (xnum.lt.0.0001) then
       write(77,16) "the cage shall not be reduced"
      else
       write(77,16) "the cage shall be reduced down to [Angst]", rfin
       write(77,16) "after a number of reductions", xnum
      end if

* ------------------------------
* vibration freq of cage atoms
* ------------------------------

* 1 fs= 1.x10^{-15} s
* sqrt(k/m)= sqrt[(N/cm)/amu]= facFreq * s^{-1}

      facFreq= 1./ sqrt(1.6605402d-29)

* omega= 2pi*freq= sqrt(k/m)  ==> freq= sqrt(k/m)/ 2pi

      freq1= facFreq * sqrt(springK/atomass(2)) /(2.*pi)  ! in 1/s
      freq2= freq1/ (100.*speedC)    ! in cm^-1
      period= 1./(freq1 *1.0d-15)    ! in fs

      write(77,*)
      write(77,*)  "vibrations of cage atoms (void space)"
      write(77,17) "vib freq. [1/s] ", freq1
      write(77,16) "vib freq. [1/cm]", freq2
      write(77,16) "vib period [fs]", period
      write(77,*)

      v1= Time/dt
      nsteps= int(v1)  ! number of time steps to compute
      write(77,*) "energy calculations to perform", nsteps

* -----------------
* non-used variables
* -----------------

      rfac= 1.d0    ! non-used variables
      vfac= 1.d0    ! are given numbers

      write(77,*)
      write(77,*) "non-used variables have values:"
      write(77,16) "rfac", rfac
      write(77,16) "vfac", vfac

      write(77,*) "* -------------------------------"
      write(77,*) "   wall fluctuations in progress"
      write(77,*) "* -------------------------------"

      write(77,*) "CASI AL FINAL DE LA SUBRUTINA MOLINFO"

   13 format (1x,i3,1x,a,1x,i3,1x,a,4(1x,f9.4))
   14 format (1x,a,1x,i7)
   16 format (1x,a,1x,f11.3)
   17 format (1x,a,1x,e10.4)

      return
      end
* *******************************
* **************************************************
* given atom symbol or atom number, this module finds
* (atom symbol, atom number, atom name, atom mass)
* using the periodic table of the elements
* **************************************************
      subroutine periodicTable (i,atsymbol,atnum,atname,atmass)
      implicit double precision (a-h,o-z)
      parameter (m1=1500)

      CHARACTER atsymbol(m1)*5, atname*12
      integer atnum

* atom masses in au
* 1 atomic mass unit (amu) = 1.6605402 x 10^{-27} Kg

* WARNING: be careful in other programs with
* implicit real*8 (a-h,o-z)

      if ((atsymbol(i).eq.'H').or.(atnum.eq.1)) then
      atnum= 1
      atname= "Hydrogen"
      atsymbol(i)= 'H'
      atmass= 1.00794
      end if

      if ((atsymbol(i).eq.'He').or.(atnum.eq.2)) then
      atnum= 2
      atname= "Helium"
      atsymbol(i)= 'He'
      atmass= 4.002602
      end if

      if ((atsymbol(i).eq.'Li').or.(atnum.eq.3)) then
      atnum= 3
      atname= "Lithium"
      atsymbol(i)= 'Li'
      atmass= 6.941
      end if

      if ((atsymbol(i).eq.'Be'.or.(atnum.eq.4))) then
      atnum= 4
      atname= "Beryllium"
      atsymbol(i)= 'Be'
      atmass= 9.012182
      end if

      if ((atsymbol(i).eq.'B').or.(atnum.eq.5)) then
      atnum= 5
      atname= "Boron"
      atsymbol(i)= 'B'
      atmass= 10.811
      end if

      if ((atsymbol(i).eq.'C').or.(atnum.eq.6)) then
      atnum= 6
      atname= "Carbon"
      atsymbol(i)= 'C'
      atmass= 12.0107
      end if

      if ((atsymbol(i).eq.'N').or.(atnum.eq.7)) then
      atnum= 7
      atname= "Nitrogen"
      atsymbol(i)= 'N'
      atmass= 14.003070
      end if

      if ((atsymbol(i).eq.'O').or.(atnum.eq.8)) then
      atnum= 8
      atname= "Oxygen"
      atsymbol(i)= 'O'
      atmass= 15.99491
      end if

      if ((atsymbol(i).eq.'F').or.(atnum.eq.9)) then
      atnum= 9
      atname= "Fluorine"
      atsymbol(i)= 'F'
      atmass= 18.9984032
      end if

      if ((atsymbol(i).eq.'Ne').or.(atnum.eq.10)) then
      atnum= 10
      atname= "Neon"
      atsymbol(i)= 'Ne'
      atmass= 20.1797
      end if

      if ((atsymbol(i).eq.'Na').or.(atnum.eq.11)) then
      atnum= 11
      atname= "Sodium"
      atsymbol(i)= 'Na'
      atmass= 22.989770
      end if

      if ((atsymbol(i).eq.'Mg').or.(atnum.eq.12)) then
      atnum= 12
      atname= "Magnesium"
      atsymbol(i)= 'Mg'
      atmass= 23.985040
      end if

      if ((atsymbol(i).eq.'Al').or.(atnum.eq.13)) then
      atnum= 13
      atname= "Aluminum"
      atsymbol(i)= 'Al'
      atmass= 26.981540
      end if

      if ((atsymbol(i).eq.'Si').or.(atnum.eq.14)) then
      atnum= 14
      atname= "Silicon"
      atsymbol(i)= 'Si'
      atmass= 27.976930
      end if

      if ((atsymbol(i).eq.'P').or.(atnum.eq.15)) then
      atnum= 15
      atname= "Phosphorous"
      atsymbol(i)= 'P'
      atmass= 30.973760
      end if

      if ((atsymbol(i).eq.'S').or.(atnum.eq.16)) then
      atnum= 16
      atname= "Sulfur"
      atsymbol(i)= 'S'
      atmass= 32.066
      end if

      if ((atsymbol(i).eq.'Cl').or.(atnum.eq.17)) then
      atnum= 17
      atname= "Chlorine"
      atsymbol(i)= 'Cl'
      atmass= 35.4527
      end if

      if ((atsymbol(i).eq.'Ar').or.(atnum.eq.18)) then
      atnum= 18
      atname= "Argon"
      atsymbol(i)= 'Ar'
      atmass= 39.948
      end if

      if ((atsymbol(i).eq.'K').or.(atnum.eq.19)) then
      atnum= 19
      atname= "Potassium"
      atsymbol(i)= 'K'
      atmass= 39.0983
      end if

      if ((atsymbol(i).eq.'Ca').or.(atnum.eq.20)) then
      atnum= 20
      atname= "Calcium"
      atsymbol(i)= 'Ca'
      atmass= 40.078
      end if

      if ((atsymbol(i).eq.'Sc').or.(atnum.eq.21)) then
      atnum= 21
      atname= "Scandium"
      atsymbol(i)= 'Sc'
      atmass= 44.955910
      end if

      if ((atsymbol(i).eq.'Ti').or.(atnum.eq.22)) then
      atnum= 22
      atname= "Titanium"
      atsymbol(i)= 'Ti'
      atmass= 47.867
      end if

      if ((atsymbol(i).eq.'V').or.(atnum.eq.23)) then
      atnum= 23
      atname= "Vanadium"
      atsymbol(i)= 'V'
      atmass= 50.9415
      end if

      if ((atsymbol(i).eq.'Cr').or.(atnum.eq.24)) then
      atnum= 24
      atname= "Chromium"
      atsymbol(i)= 'Cr'
      atmass= 51.9961
      end if

      if ((atsymbol(i).eq.'Mn').or.(atnum.eq.25)) then
      atnum= 25
      atname= "Manganese"
      atsymbol(i)= 'Mn'
      atmass= 54.938049
      end if

      if ((atsymbol(i).eq.'Fe').or.(atnum.eq.26)) then
      atname= "Iron"
      atsymbol(i)= 'Fe'
      atmass= 55.845
      end if

      if ((atsymbol(i).eq.'Co').or.(atnum.eq.27)) then
      atnum= 27
      atname= "Cobalt"
      atsymbol(i)= 'Co'
      atmass= 58.933200
      end if

      if ((atsymbol(i).eq.'Ni').or.(atnum.eq.28)) then
      atnum= 28
      atname= "Nickel"
      atsymbol(i)= 'Ni'
      atmass= 58.6934
      end if

      if ((atsymbol(i).eq.'Cu').or.(atnum.eq.29)) then
      atnum= 29
      atname= "Copper"
      atsymbol(i)= 'Cu'
      atmass= 63.546
      end if

      if ((atsymbol(i).eq.'Zn').or.(atnum.eq.30)) then
      atnum= 30
      atname= "Zinc"
      atsymbol(i)= 'Zn'
      atmass= 65.39
      end if

      if ((atsymbol(i).eq.'Ga').or.(atnum.eq.31)) then
      atnum= 31
      atname= "Gallium"
      atsymbol(i)= 'Ga'
      atmass= 69.723
      end if

      if ((atsymbol(i).eq.'Ge').or.(atnum.eq.32)) then
      atnum= 32
      atname= "Germanium"
      atsymbol(i)= 'Ge'
      atmass= 72.61
      end if

      if ((atsymbol(i).eq.'As').or.(atnum.eq.33)) then
      atnum= 33
      atname= "Arsenic"
      atsymbol(i)= 'As'
      atmass= 74.92160
      end if

      if ((atsymbol(i).eq.'Se').or.(atnum.eq.34)) then
      atnum= 34
      atname= "Selenium"
      atsymbol(i)= 'Se'
      atmass= 78.96
      end if

      if ((atsymbol(i).eq.'Br').or.(atnum.eq.35)) then
      atnum= 35
      atname= "Bromine"
      atsymbol(i)= 'Br'
      atmass= 79.904
      end if

      if ((atsymbol(i).eq.'Kr').or.(atnum.eq.36)) then
      atnum= 36
      atname= "Kripton"
      atsymbol(i)= 'Kr'
      atmass= 83.80
      end if

      if ((atsymbol(i).eq.'Rb').or.(atnum.eq.37)) then
      atnum= 37
      atname= "Rubidium"
      atsymbol(i)= 'Rb'
      atmass= 85.4678
      end if

      if ((atsymbol(i).eq.'Sr').or.(atnum.eq.38)) then
      atnum= 38
      atname= "Strontium"
      atsymbol(i)= 'Sr'
      atmass= 87.62
      end if

      if ((atsymbol(i).eq.'Y').or.(atnum.eq.39)) then
      atnum= 39
      atname= "Yttrium"
      atsymbol(i)= 'Y'
      atmass= 88.90585
      end if

      if ((atsymbol(i).eq.'Zr').or.(atnum.eq.40)) then
      atnum= 40
      atname= "Zirconium"
      atsymbol(i)= 'Zr'
      atmass= 91.224
      end if

      if ((atsymbol(i).eq.'Nb').or.(atnum.eq.41)) then
      atnum= 41
      atname= "Niobium"
      atsymbol(i)= 'Nb'
      atmass= 92.90638
      end if

      if ((atsymbol(i).eq.'Mo').or.(atnum.eq.42)) then
      atnum= 42
      atname= "Molybdenum"
      atsymbol(i)= 'Mo'
      atmass= 95.94
      end if

      if ((atsymbol(i).eq.'Tc').or.(atnum.eq.43)) then
      atnum= 43
      atname= "Technetium"
      atsymbol(i)= 'Tc'
      atmass= 98
      end if

      if ((atsymbol(i).eq.'Ru').or.(atnum.eq.44)) then
      atnum= 44
      atname= "Ruthenium"
      atsymbol(i)= 'Ru'
      atmass= 101.07
      end if

      if ((atsymbol(i).eq.'Rh').or.(atnum.eq.45)) then
      atnum= 45
      atname= "Rhodium"
      atsymbol(i)= 'Rh'
      atmass= 102.90550
      end if

      if ((atsymbol(i).eq.'Pd').or.(atnum.eq.46)) then
      atnum= 46
      atname= "Palladium"
      atsymbol(i)= 'Pd'
      atmass= 106.42
      end if

      if ((atsymbol(i).eq.'Ag').or.(atnum.eq.47)) then
      atnum= 47
      atname= "Silver"
      atsymbol(i)= 'Ag'
      atmass= 107.8682
      end if

      if ((atsymbol(i).eq.'Cd').or.(atnum.eq.48)) then
      atnum= 48
      atname= "Cadmium"
      atsymbol(i)= 'Cd'
      atmass= 112.411
      end if

      if ((atsymbol(i).eq.'In').or.(atnum.eq.49)) then
      atnum= 49
      atname= "Indium"
      atsymbol(i)= 'In'
      atmass= 114.818
      end if

      if ((atsymbol(i).eq.'Sn').or.(atnum.eq.50)) then
      atnum= 50
      atname= "Tin"
      atsymbol(i)= 'Sn'
      atmass= 118.710
      end if

      if ((atsymbol(i).eq.'Sb').or.(atnum.eq.51)) then
      atnum= 51
      atname= "Antimony"
      atsymbol(i)= 'Sb'
      atmass= 121.76
      end if

      if ((atsymbol(i).eq.'Te').or.(atnum.eq.52)) then
      atnum= 52
      atname= "Tellurium"
      atsymbol(i)= 'Te'
      atmass= 127.60
      end if

      if ((atsymbol(i).eq.'I').or.(atnum.eq.53)) then
      atnum= 53
      atname= "Iodine"
      atsymbol(i)= 'I'
      atmass= 126.90447
      end if

      if ((atsymbol(i).eq.'Xe').or.(atnum.eq.54)) then
      atnum= 54
      atname= "Xenon"
      atsymbol(i)= 'Xe'
      atmass= 131.29
      end if

      
      return
      end
* ************************************************
* *************************
* perturb the confined atoms
* *************************
      subroutine perturb (npar,nstatic, x,y,z, vx,vy,vz,
     &    xnew,ynew,znew, vxnew,vynew,vznew, rfac,vfac)

      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)

      dimension  xnew(nmax),  ynew(nmax),  znew(nmax)
      dimension vxnew(nmax), vynew(nmax), vznew(nmax)

      do i= 1, npar ! <-------
                              !
      if (i.le.nstatic) then  !
       xnew(i)= x(i)          !
       ynew(i)= y(i)          ! in Angst
       znew(i)= z(i)          !
                              !
       vxnew(i)= vx(i)        !
       vynew(i)= vy(i)        ! in Angst/fs
       vznew(i)= vz(i)        !
      endif                   !
                              !
                              !
      if (i.gt.nstatic) then  !
       xnew(i)= x(i) *rfac    !
       ynew(i)= y(i) *rfac    ! in Angst
       znew(i)= z(i) *rfac    !
                              !
       vxnew(i)= vx(i) *vfac  !
       vynew(i)= vy(i) *vfac  ! in Angst/fs
       vznew(i)= vz(i) *vfac  !
      endif                   !
                              !
      end do ! <--------------

      return
      end
* ****************************
* **********************************
* get polar coordinates (rho,teta,fi)
* from the cartesian coordinates (a,b,c)
* **********************************
      subroutine polar (a,b,c,rho,teta,fi)
      implicit double precision (a-h,o-z)

        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')
      

      pi= 3.1415926535898
      rad= pi/180.0

*  0 <= rho
*  0 <= teta <= 2pi    teta is in the (x,y) plane
*  0 <= fi <= pi       goes from +z to -z

      rho= sqrt(a**2 + b**2 + c**2)
      fi= acos(c/rho)

* get teta in the 1st quadrant

      if ((a.ge.0.).and.(b.eq.0.)) then
      teta= 0.d0
      end if

      if ((a.gt.0.).and.(b.gt.0.)) then
      teta= atan(b/a)
      end if

      if ((a.eq.0.).and.(b.gt.0.)) then
      teta=  pi/2.d0
      end if

* get teta in the 2nd quadrant

      if ((a.lt.0.).and.(b.gt.0.)) then
      teta= pi - atan(abs(b/a))
      end if

      if ((a.lt.0.).and.(b.eq.0.)) then
      teta= pi
      end if

* get teta in the 3rd quadrant

      if ((a.lt.0.).and.(b.lt.0.)) then
      teta= pi + atan(abs(b/a))
      end if

      if ((a.eq.0.).and.(b.lt.0.)) then
      teta= 3./2.* pi
      end if

* get teta in the 4th quadrant

      if ((a.gt.0.).and.(b.lt.0.)) then
      teta= 2.*pi - atan(abs(b/a))
      end if

* the point has no teta angle

      if ((a.eq.0.).and.(b.eq.0.)) then
      write(77,*) "ERROR: x=0 and y=0"
      end if

      return
      end
* **********************************
* ***********************************
* store data of the dynamics for a restart
*
* WARNING:
* all arrays & summations in memory are lost after a crash,
* only information written to output files is preserved
* make a backup of all files before making a restart,
* ***********************************
      subroutine restartOlla (ilabel,t,npar, x,y,z, vx,vy,vz,
     &    Gradx,Grady,Gradz,x0,y0,z0, refPE)

      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      CHARACTER name

      dimension x0(nmax), y0(nmax), z0(nmax)
      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)
      dimension Gradx(nmax), Grady(nmax), Gradz(nmax)

      CHARACTER (LEN = 400) NombreDinamica
      CHARACTER (LEN = 400) RutaArchivosEntrada
      CHARACTER (LEN = 400) RutaArchivosSalida
      CHARACTER (LEN = 400) RestartDat


* ilabel= 1 ==> write information
* ilabel= 2 ==> read  information

* ------------------------
* open the file for storing data
* ------------------------

      call DatosGrlsDinamica(NombreDinamica,
     & RutaArchivosEntrada, RutaArchivosSalida)

      RestartDat = trim(RutaArchivosSalida)
     & // '/RestartDat_' // trim(NombreDinamica)
     & // '.dat'

c     open (unit=15, file= 'restart.dat',status= 'unknown')
      open (unit=15, file= trim(RestartDat), status= 'unknown')

* ------------------------
* write parameters & variables
* ------------------------

      if (ilabel.eq.1) then ! <---
                                  !
      write(15,*) "time", t       !
      write(15,*) "referencePE", refPE
      write(15,*) "positions"     !
                                  !
      do i= 1, npar               !
      write(15,*) i, x(i), y(i), z(i)
      end do                      !
                                  !
      write(15,*) "speeds"        !
                                  !
      do i= 1, npar               !
      write(15,*) i, vx(i), vy(i), vz(i)
      end do                      !
                                  !
      write(15,*) "gradients"     !
                                  !
      do i= 1, npar               !
      write(15,*) i, Gradx(i), Grady(i), Gradz(i)
      end do                      !
                                  !
      write(15,*) "starting positions"
                                  !
      do i= 1, npar               !
      write(15,*) i, x0(i), y0(i), z0(i)
      end do                      !
                                  !
      end if ! <------------------

* ------------------------
* read parameters & variables
* ------------------------

      if (ilabel.eq.2) then ! <---
                                  !
      read(15,*) name, t          !
      read(15,*) name, refPE      !
      read(15,*) name             !
                                  !
      do i= 1, npar               !
      read(15,*) k, x(i), y(i), z(i)
      end do                      !
                                  !
      read(15,*) name             !
                                  !
      do i= 1, npar               !
      read(15,*) k, vx(i), vy(i), vz(i)
      end do                      !
                                  !
      read(15,*) name             !
                                  !
      do i= 1, npar               !
      read(15,*) k, Gradx(i), Grady(i), Gradz(i)
      end do                      !
                                  !
      read(15,*) name             !
                                  !
      do i= 1, npar               !
      read(15,*) k, x0(i), y0(i), z0(i)
      end do                      !
                                  !
      end if ! <------------------

      close (15)

      return
      end
* *******************************
* ********************
* module to write the nwchem script
* the file created here shall be called from
* the main program to run nwchem
* ********************
      subroutine scriptNwchem
      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)

        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut
        CHARACTER (LEN = 400) NwchemOut
        CHARACTER (LEN = 400) NwchemInp
        CHARACTER (LEN = 400) NwchemScript

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

        NwchemOut = trim(RutaArchivosSalida)
     & // '/NwchemOut_' // trim(NombreDinamica)
     & //'.out'
       
        NwchemInp = trim(RutaArchivosEntrada)
     & // '/NwchemInp_' // trim(NombreDinamica)
     & //'.inp'
        
        NwchemScript = trim(RutaArchivosEntrada)
     & // '/NwchemScript_' // trim(NombreDinamica)

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')


* char(34)= "
* char(35)= #

      write(77,*) "creating the script to run nwchem"

c     open(14,file="script.nwchem",status="unknown")
      open(14,file=trim(NwchemScript),status="unknown")

      write(14,*)
      write(14,*) char(35), " script to submit nwchem"
      write(14,*) char(35), " commands may be changed in module"
      write(14,*) char(35), " scriptNwchem"
      write(14,*)
c     write(14,*) "  echo  "
      write(14,*) "  echo  ", char(34),
     &      "!!!!!! nwchem is working now !!!!!!", char(34)
      write(14,*)
      write(14,*) char(35), " submit the nwchem job"

c     write(14,*) "nwchem nwchem.inp > nwchem.out"
      write(14,*) "nwchem "//trim(NwchemInp)//">"//trim(NwchemOut) 
      write(14,*)

      write(14,*) char(35), " get the dft ene"
c     write(14,*) "  egrep  ", char(34), "DFT energy", char(34),
c    &   " nwchem.out"
      write(14,*) "  egrep  ", char(34), "DFT energy", char(34),
     &  trim(NwchemOut) 


      write(14,*)

      write(14,*) "  echo  ", char(34),
     &      "!!!!!! nwchem finished this step !!!!!!", char(34)
c     write(14,*)  "  tail -1  ", " nwchem.out"
      write(14,*)  "  tail -1  ", trim(NwchemOut) 
      write(14,*) "  echo  "
      write(14,*)

      close (14)

c     call system ("chmod a+rwx script.nwchem") ! make it executable
      call system ("chmod a+rwx " // trim(NwchemScript)) ! make it executable

      return
      end
* ********************
* ***********************************
* shrink the cage and atoms inside the cage
* ***********************************
      subroutine shrink (npar,ncage,x0,y0,z0,x,y,z,radiusf)
      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension x0(nmax), y0(nmax), z0(nmax)


        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')
      


 
* --------------------------------------------
* resize the reference radius of cage particles
* --------------------------------------------

      sum1= 0.d0
      sum2= 0.d0

      do i= 1, ncage ! <----
                            !
      a= x0(i)              !
      b= y0(i)              !
      c= z0(i)              !
                            !
      dold= sqrt (a**2 +b**2 +c**2)
                            !
                            ! transform to polar coords
      call polar (a,b,c,rho,teta,fi)
                            !
      rho= rho *radiusf     ! scale radial distances
                            !
                            ! back to cartesian coords
      call cartesian (rho,teta,fi,a,b,c)
                            !
      x0(i)= a              ! store reference
      y0(i)= b              ! positions
      z0(i)= c              !
                            !
      x(i)= a               ! present positions
      y(i)= b               !
      z(i)= c               !
                            !
      dnew= sqrt (a**2 +b**2 +c**2)
                            !
      sum1= sum1 + dold     ! sum distances
      sum2= sum2 + dnew     ! of cage atoms
                            !
      end do ! <------------

* average radius

      rold= sum1 /dble(ncage)
      rnew= sum2 /dble(ncage)

      write(77,"(1x,a,1x,f9.4)") "old reference rad [Angst]", rold
      write(77,"(1x,a,1x,f9.4)") "new reference rad [Angst]", rnew
      write(77,"(1x,a,1x,f9.4)") "shrinking factor", radiusf

* --------------------------------------------
* resize the present radii of confined particles
* --------------------------------------------

      do i= ncage+1, npar ! <--
                               !
      a= x(i)                  ! present positions
      b= y(i)                  !
      c= z(i)                  !
                               !
      call polar (a,b,c,rho,teta,fi)
                               !
      rho= rho *radiusf        !
                               !
      call cartesian (rho,teta,fi,a,b,c)
                               !
      x(i)= a                  ! scaled positions
      y(i)= b                  !
      z(i)= c                  !
                               !
      end do ! <---------------

      return
      end
* ***********************************
* ***********************************
* slow down speeds of confined atoms close to the cage
* to simulate stochastic loss of kinetic energy, thus
* simulating a cooling process
* ***********************************
      subroutine slowSpeeds (t,npar,ncage, x,y,z, vx,vy,vz,
     &     heatup,cooldown)

      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)

      CHARACTER*3 heatup, cooldown

      integer ilabel(nmax)


        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')    




* --------------------
* initialize parameters
* --------------------

      do i= 1, npar ! labels for atoms whose
      ilabel(i)= 0  ! speeds shall be reduced
      end do

      xlambda= 0.93 ! factor to reduce atom speeds
      dthresh= 3.50 ! threshold distance to the cage

      if (t.eq.1.) then
      write(77,12) "factor to reduce atom speeds", xlambda
      write(77,12) "distance from cage for cooling down", dthresh
      end if

   12 format (1x,a,2x,f6.3)

c     write(77,*)
c     write(77,*) "atoms with scaled speeds at time=", t
c     write(77,*) "atom number         old speed         new speed"

* --------------------
* classify atoms
* --------------------

      icont= 0

      do j= 1, ncage ! <-----  atoms of the cage
      do i= ncage+1, npar    ! atoms in the cage
                             !
       dx= (x(i)-x(j))**2    !
       dy= (y(i)-y(j))**2    !
       dz= (z(i)-z(j))**2    !
       dist= sqrt(dx+dy+dz)  ! find atoms close
                             ! to the cage
                             !
      if ((dist.lt.dthresh).and.(ilabel(i).eq.0)) then
       icont= icont +1       !
       Vold= sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
                             !
       vx(i)= vx(i) *xlambda !
       vy(i)= vy(i) *xlambda ! scale down
       vz(i)= vz(i) *xlambda ! these speeds
                             !
       ilabel(i)= i          ! tag this atom
       Vnew= xlambda* Vold   ! new speed
                             !
c     write(77,14) i, Vold, Vnew
      end if                 !
                             !
      end do                 !
      end do ! <-------------

      write(77,*) "num of cooled atoms", icont

* reset switches

      heatup=   "off"
      cooldown= "off"

   14 format (1x,i5,2(1x,f9.5))

      return
      end
* *****************************
* ***********************************
* evolve positions & speeds using the Taylor expansion
* ***********************************
      subroutine taylor (npar, nstatic, atomass,
     &         x,y,z, vx,vy,vz, Gradx,Grady,Gradz,
     &         xnew,ynew,znew, vxnew,vynew,vznew, dt)

      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      dimension  atomass(nmax)
      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)
      dimension Gradx(nmax), Grady(nmax), Gradz(nmax)

      dimension  xnew(nmax),  ynew(nmax),  znew(nmax)
      dimension vxnew(nmax), vynew(nmax), vznew(nmax)

      
        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')


      write(77,*) "evolving coords and speeds with Taylor"

* -------------------------
*  accel in terms of dE/dx
* -------------------------

* Fx= -dE/dx

*  Xnew= x + vx*dt +ax*dt^2/2
*      = x + vx*dt -(1/m)*Gradx*dt^2/2

* Vxnew= vx +ax*dt
*      = vx -(1/m)*Gradx*dt

* factor to change units

      faccl= 0.49616184

* ----------------------------------------
* progress positions & speeds of confined atoms
* ----------------------------------------

      do i= nstatic+1, npar ! <--
                                 ! in Angst/fs^2
      ax= -faccl *Gradx(i) /atomass(i)
      ay= -faccl *Grady(i) /atomass(i)
      az= -faccl *Gradz(i) /atomass(i)
                                 !
      vxnew(i)= vx(i) + ax*dt    ! in Angst/fs
      vynew(i)= vy(i) + ay*dt    !
      vznew(i)= vz(i) + az*dt    !
                                 ! in Angst
      xnew(i)= x(i) + vx(i)*dt + ax *dt**2 /2.d0
      ynew(i)= y(i) + vy(i)*dt + ay *dt**2 /2.d0
      znew(i)= z(i) + vz(i)*dt + az *dt**2 /2.d0
                                 !
      end do ! <-----------------

      return
      end
* *******************************

* **************************************
* informative module for the transformation of units
* **************************************

* -------------------------
* physical variables & their unis
* -------------------------

* [position]= Angts              [mass]= Kg
* [speeds]=   Angst/fs           [time]= fs
* [acceleration]= Angst/fs^2     [Force]= Newton
*                                [Energy]= Joules

* 1 Angst= 1x10^{-10} mt
* 1 mt=    1x10^{+10} Angst
* 1 Bohr= 0.529178 Angst= 0.529178 x10^{-10} mt
* 1 amu=   1.6605402 x10^{-27} Kg
* 1 Hartree= 4.35988 x10^{-18} J
* 1 Hartree= 27.211396132 eV
* 1 fs= 1x10^{-15} s
* 1 atmosphere= 101325.0 pascals
* 1 pascal= 1/101325.0 atmosphere

* Newton= Kg *mt/ s^2
* Joule= Newton * mt= Kg *(mt/s)^2

* kB= 1.3806504x10^{-23}    in Joule/Kelvin
* pi= 3.1415926535898
* c= 2.99792458 *10^8 m/s= 2.99792458*10^{10} cm/s      speed of light
* mole= 6.02214179 *10^{23}
* R= 8.3144621 J/(mol*Kelvin)   universal constant of gases
* Planck constant= hbar= 6.62606957 * 10^{-34} (kg*m^2)/s
* --------------------------------------
* units of positon= acceleration * dt^2
* --------------------------------------
* [atom mass]= atomass x amu
* [Force]= [dE/dx]= Hartree/Bohr
* [dt]= fs

* rx= ax * dt^2 /2
*  = (Fx/m) * dt^2 /2
*  = -(1/m) (dE/dx) *dt^2 /2
*  = -(1/m) (Gradx) *dt^2 /2  [(1/amu)*Hartree/Bohr *fs^2]

* [(1/amu) *Hartree/Bohr *fs^2]
* = 1/(1.6605402 x10^{-27} Kg) *
*     4.35988 x10^{-18} J/ 0.529178 x10^{-10} mt *1x10^{-30} s^2
* = 4.35988 /(1.6605402 * 0.529178) x10^{+27-18+10-30} J/(Kg*mt)*s^2
* = 4.9616184 x10^{-11} m

* transformation factor for positions derived from the force

* faccl= 0.49616184    dimensionless

* faccl*Angst= (1/amu)*(Hartree/Bohr)*(fs^2)

* ----------------------------------
* units of speed= acceleration * dt
* ----------------------------------
* from the result for positions we have:

* fvel= 0.49616184

* fvel *(Angst/fs)= (1/amu) *(Hartree/Bohr) *fs

* ------------------------
* units of kinetic energy
* ------------------------
* [Ek]= [mass *vx(i)^2]= amu *(Angst/fs)^2

* 1 amu = 1.6605402 x10^{-27} Kg
* 1 Angst= 1x10^{-10} mt
* 1 fsec= 1x10^{-15} sec

* [Ek]= [1.6605402x10^{-27} Kg] [1x10^{-10} mt]^2/ [1x10^{-15} sec]^2

* Joule= Newton*mt= Kg*(mt/s)^2

* transformation factor between units of kinetic energy

* facEkin= 1.6605402d-17

* facEkin * Joule= amu *(Angst/fs)^2

* -----------------------------
* units of pressure= force/area
* -----------------------------
* radial force= F . R/|R|= -(dE/dx,dE/dy,dE/dz) . R/|R| [Hartree/Bohr]

* [Hartree/Bohr]= 4.35988 x10^{-18} J / 0.529178 x10^{-10} mt
* = (4.35988/ 0.529178) x10^{-8} [J/m=Newton]
* = 8.2389668504738d-08 [N]

* area [Angst^2]= 1x10^{-20} mt^2

* pressure= radial force/area= 8.2389668504738d+12 [N/m^2=Pa]

* facPres= 8.2389668504738d+12

* facPres * Pa= (Hartree/Bohr)/(Angst^2)

* Joule/Angst^3= Joule/[1x10^{-30} mt^3]
*              = 1x10^{30} Newton/mt^2
*              = 1x10^{30} Pascal

* facPres2= 1.x10^{30}

* facPres2 * Pa= Joule/Angst^3

* -------------
* units of kT
* -------------
* 1 Joule= 2.29371248 x 10^{17}   Hartree

* kB= 1.3806504 x 10^{-23}       Joule/Kelvin
*   = 1.3806504x10^{-23} * 2.29371248x10^{17}   Hartree/Kelvin
*   = 3.16681505 x 10^{-6}     Hartree/Kelvin

* [Temp]= Kelvin
* [kB*Temp]= Joule

* 1 amu= 1.6605402 x10^{-27} Kg
* m/s= 1x10^{+10)/1x10^{+15) Angst/fs= 1x10^{-5} Angst/fs

* [kB*Temp/amu]= J/ (1.6605402x10^{-27} Kg)
*   = (m/s)^2/ (1.6605402x10^{-27})
*   = (Angst/fs)^2 * 1x10^{-10})/ (1.6605402x10^{-27})
*   = 1x10^{+17}/1.6605402 (Angst/fs)^2

* facKT= 1.0d17 /1.6605402

* facKT *(Angst/fs)^2= (Joule/Kelvin)*(Kelvin/amu)

* -------------------------
* units of spring constant K
* -------------------------
* [spring K]= Newton/cm

* 1 amu= 1.6605402 x10^{-27} Kg

* Newton/cm= Kg*mt/(s^2*cm)
*          = Kg*100cm/ (1x10^{30} fs^2 *cm)
*          = Kg/(1x10^{28} fs^2)
*          = amu/(1.6605402 x10^{-27} *1x10^{28} fs^2)
*          = (amu/fs^2) *1.0/16.605402

* fspring= 1.0/ 16.605402

* Newton/cm = fspring *(amu/fs^2)

* -------------------------
* units of thermostat mass
* -------------------------
* kB= 1.3806504 x 10^{-23}       Joule/Kelvin
* 1 Joule= 2.29371248 x 10^{17} Hartree
* 1 mole= 6.02214179 x 10^{23}

* Mt= (KJoule/Mole) x ps^2
*   = 1000 x (2.29371248 x 10^{17} Hartree) x (10^{-12} sec)^2 /6.02214179 x 10^{23}
*   = (2.29371248/ 6.02214179) x 10^(3+17-24-23)   Hartree * sec^2
*   = 3.808798530 x 10^{-27} Hartree * sec^2

* fmt= 3.808798530d-27
* (KJoule/Mole) x ps^2= fmt * (Hartree x sec^2)

* -------------------------
* units of frequency
* -------------------------
* (N/cm)/amu= (Kg mt/s^2)/ (cm * 1.6605402x10^{-27} Kg)
*           = (100 cm/s^2)/(cm * 1.6605402x10^{-27})
*           = (1/s^2)/1.6605402x10^{-29}

* [freq]= sqrt ([spring K/mass]
*       = sqrt [(N/cm)/amu]
*       = 1/sqrt(1.6605402x10^{-29}) s^{-1}

* facFreq= 1./sqrt(1.6605402d-29)

* sqrt[(N/cm)/amu]= facFreq * s^{-1}
* -------------------------------------------
* units for speed introducing the universal constant of gases
* -------------------------------------------
* universal constant of gases   R=8.3144621 J/(mol*Kelvin)
* Joule= newton*m= Kg*m/s^2 *m= Kg*m^2/s^2 
* mole= 6.02214179 *10^{23}

* 1 amu= 1.6605402 x10^{-27} Kg= 1.6605402*10^{-27} Kg *6.02214179*10^{23}/ mol 
*      = 1.6605402 *6.02214179 *10^{-4} Kg/mol 
*      = 1*10^{-3} Kg/mol 

* [sqrt(RT/mass)]= sqrt [J/(mol*Kelvin) * Kelvin/(1*10^{-3}*Kg/mol)]
*                = sqrt [10^3 (m/s)^2]
*                = sqrt [10^3] m/s

* m/s= 1x10^{+10)/1x10^{+15) Angst/fs= 1x10^{-5} Angst/fs
 
* fac= 10^{-5}* sqrt(1000.)

* [sqrt(RT/m)]= sqrt [J/(mol*Kelvin) * Kelvin/amu]= fac *Angst/fs
* -------------------------------------------
* convert energy to wavelength
* -------------------------------------------
* energy= hbar*freq
* freq= c/lambda

* energy= hbar* c/lambda
* lambda= c/freq= c* hbar/energy

* hbar= 6.62606957 * 10^{-34} (kg*m^2)/s
* c= 2.99792458 *10^8 m/s
* 1 Hartree= 4.35988 x10^{-18} J

* lambda= c* hbar/energy
*       = [2.99792458 *10^8 m/s] * [6.62606957 * 10^{-34} (kg*m^2)/s] / (alpha hartree)
*       = (2.99792458* 6.62606957/4.35988) *10^{-8} kg*m^3/s^2 / (alpha Joules)
*       = (4.55619348)/alfa *10^{-8} m

* calculator: frequency-wavelength-energy
* http://www2.chemistry.msu.edu/faculty/reusch/virttxtjml/cnvcalc.htm

* ionizing & non-ionizing radiation information
* http://www.epa.gov/radiation/understand/ionize_nonionize.html
* ---------------------------------
* sequence of events in breaking a bond (bb)
* ---------------------------------
* time intervals         wave function

* t is in [1,bb_1]           2    reuse the WF

* t= bb_1                    4   break a bond & lose charge
* t in [bb_1, bb_1+10]       5    reuse the WF
* t= bb_1+10                 1   recover the charge
* t in [bb_1+11, bb_2]       2    reuse the WF

* t= bb_2                    4   break a 2nd bond & lose charge
* t in [bb_2, bb_2+10]       5    reuse the WF
* t= bb_2+10                 1   recover the charge
* t in [bb_2+11, bb_3]       2    reuse the WF
* ************************************

* ***********************************
* get the path of the working directory
* ***********************************
      subroutine workDir (name)
      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      CHARACTER  qline*(maxchr), name*60
      dimension  itok(maxtok), jtok(maxtok)

      external kparse

        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) NombreOllaOut

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        NombreOllaOut = trim(RutaArchivosSalida)
     & // '/outgrl_' // trim(NombreDinamica)
     & // '.out'

c       open (unit = 77, file = trim(NombreOllaOut), status = 'old')

* --------------------------
* parse the present working directory
* --------------------------

      open (unit=16, file="t1", status='unknown')
      open (unit=17, file="t2", status='unknown')

      read (16,12) qline

      ntok= kparse (qline,maxchr,itok,jtok,maxtok)
      write (17,12) char(34), qline (itok(1):jtok(1)), char(34)

      close (16)
      close (17)

* ---------------------
* read the directory name
* ---------------------

      open (unit=16, file="t2", status='unknown')
      read (16,12) name
      close (16)

   12 format (3a)

      return
      end
* *********************************
* ***********************************
* molec to be computed by terachem
* ***********************************
      subroutine writeMol (npar,atsymbol,x,y,z)
      implicit double precision (a-h,o-z)
!     include 'RutaGrl/ArchivosEntrada/parametros'
      
      parameter (nmax=3000)
      parameter (maxt=50000)
      parameter (m1=1000000)
      parameter (maxchr=100, maxtok=10)


      CHARACTER atsymbol(nmax)*5
      dimension  x(nmax), y(nmax), z(nmax)

      open (unit=14, file= 'mol.xyz', status= 'unknown')

      write(14,*) npar
      write(14,*) "molec"

      do m= 1, npar
      write(14,35) atsymbol(m), x(m),y(m),z(m)
c     write(77,35)  atsymbol(m), x(m),y(m),z(m)
      end do

   35 format (1x,a5, 3(2x,f10.3))
c  35 format (1x,a, 3(1x,f9.4))

      close (14)

      return
      end
* **********************************

* *************************
* file: randgen.f
* richard chandler
* (richard@stats.ucl.ac.uk)
* paul northrop
* (northrop@stats.ox.ac.uk)
* last modified: 26/8/03
* see file randgen.txt for details
* *************************
      block data zbqlbd01
 
* initializes seed array etc. for random number generator.
* the values below have themselves been generated using the
* nag generator.
 
      common /zbql0001/ zbqlix,b,c
      double precision zbqlix(43),b,c
      integer i

      data (zbqlix(i),i=1,43) /8.001441d7,5.5321801d8,
     +1.69570999d8,2.88589940d8,2.91581871d8,1.03842493d8,
     +7.9952507d7,3.81202335d8,3.11575334d8,4.02878631d8,
     +2.49757109d8,1.15192595d8,2.10629619d8,3.99952890d8,
     +4.12280521d8,1.33873288d8,7.1345525d7,2.23467704d8,
     +2.82934796d8,9.9756750d7,1.68564303d8,2.86817366d8,
     +1.14310713d8,3.47045253d8,9.3762426d7 ,1.09670477d8,
     +3.20029657d8,3.26369301d8,9.441177d6,3.53244738d8,
     +2.44771580d8,1.59804337d8,2.07319904d8,3.37342907d8,
     +3.75423178d8,7.0893571d7 ,4.26059785d8,3.95854390d8,
     +2.0081010d7,5.9250059d7,1.62176640d8,3.20429173d8,
     +2.63576576d8/

      data b / 4.294967291d9 /
      data c / 0.0d0 /

      end
* *************************
      subroutine zbqlini (seed)
* *************************
* to initialize the random number generator - either
* repeatably or nonrepeatably. need double precision
* variables because integer storage can't handle the
* numbers involved

* arguments
* =========
* seed     (integer, input). user-input number which generates
* elements of the array zbqlix, which is subsequently used
* in the random number generation algorithm. if seed=0,
* the array is seeded using the system clock if the
* fortran implementation allows it.

* parameters
* ==========
* lflno     (integer). number of lowest file handle to try when
* opening a temporary file to copy the system clock into.
* default is 80 to keep out of the way of any existing
* open files (although the program keeps searching till
* it finds an available handle). if this causes problems,
* (which will only happen if handles 80 through 99 are
* already in use), decrease the default value.

      integer lflno
      parameter (lflno=80)

* variables
* =========
* seed     see above
* zbqlix     seed array for the random number generator. defined
* in zbqlbd01
* b,c     used in congruential initialisation of zbqlix
* ss,mm,\}     system clock secs, mins, hours and days
* hh,dd \}
* filno     file handle used for temporary file
* init     indicates whether generator has already been initialised
     
      integer seed,ss,mm,hh,dd,filno,i
      integer init
      double precision zbqlix(43),b,c
      double precision tmpvar1,dss,dmm,dhh,ddd

      common /zbql0001/ zbqlix,b,c
      save init

     
* ensure we don't call this more than once in a program    
     
cc      if (init.ge.1) then
cc       if(init.eq.1) then
cc        write(77,1)
cc        init = 2
cc       endif
cc       return
cc      else
cc       init = 1
cc      endif
     
* if seed = 0, cat the contents of the clock into a file    
* and transform to obtain zqblix(1), then use a congr.    
* algorithm to set remaining elements. otherwise take    
* specified value of seed.    
     
*---------------------------------------------  
* nb for systems which do not support the             
* (non-standard) 'call system' command,               
* this will not work, and the first clause            
* of the following if block should be                 
* commented out.                             
*---------------------------------------------  
      if (seed.eq.0) then
*---------------------------------------------  
* comment out from here if you don't have             
* 'call system' capability ...                     
*---------------------------------------------  
       call system(' date +%s%m%h%j > zbql1234.tmp')
     
*       try all file numbers for lflno to 999     
     
       filno = lflno
 10    open(filno,file='zbql1234.tmp',err=11)
       goto 12
 11    filno = filno + 1
       if (filno.gt.999) then
        write(77,2)
        return
       endif
       goto 10
 12    read(filno,'(3(i2),i3)') ss,mm,hh,dd
       close(filno)
       call system('rm zbql1234.tmp')
       dss = dint((dble(ss)/6.0d1) * b)
       dmm = dint((dble(mm)/6.0d1) * b)
       dhh = dint((dble(hh)/2.4d1) * b)
       ddd = dint((dble(dd)/3.65d2) * b)
       tmpvar1 = dmod(dss+dmm+dhh+ddd,b)
*---------------------------------------------   
* ... to here (end of commenting out for                   
* tab users without 'call system' capability                  
*---------------------------------------------   
      else
       tmpvar1 = dmod(dble(seed),b)
      endif
      zbqlix(1) = tmpvar1
      do 100 i = 2,43
       tmpvar1 = zbqlix(i-1)*3.0269d4
       tmpvar1 = dmod(tmpvar1,b)      
       zbqlix(i) = tmpvar1
 100  continue

 1    format(//5x,'****warning**** you have called routine zbqlini ',
     +'more than',/5x,'once. i''m ignoring any subsequent calls.',//)
 2    format(//5x,'**** error **** in routine zbqlini, i couldn''t',
     +' find an',/5x,
     +'available file number. to rectify the problem, decrease the ',
     +'value of',/5x,
     +'the parameter lflno at the start of this routine (in file ',
     +'randgen.f)',/5x,
     +'and recompile. any number less than 100 should work.')
      end
*---------------------------------------------
      function zbqlu01()
 
*       returns a uniform random number between 0 & 1, using
*       a marsaglia-zaman type subtract-with-borrow generator.
*       uses double precision, rather than integer, arithmetic
*       throughout because mz's integer constants overflow
*       32-bit integer storage (which goes from -2^31 to 2^31).
*       ideally, we would explicitly truncate all integer
*       quantities at each stage to ensure that the double
*       precision representations do not accumulate approximation
*       error; however, on some machines the use of dnint to
*       accomplish this is *seriously* slow (run-time increased
*       by a factor of about 3). this double precision version
*       has been tested against an integer implementation that
*       uses long integers (non-standard and, again, slow) -
*       the output was identical up to the 16th decimal place
*       after 10^10 calls, so we're probably ok ...
 
      double precision zbqlu01,b,c,zbqlix(43),x,b2,binv
      integer curpos,id22,id43

      common /zbql0001/ zbqlix,b,c
      save /zbql0001/
      save curpos,id22,id43
      data curpos,id22,id43 /1,22,43/

      b2 = b
      binv = 1.0d0/b
 5    x = zbqlix(id22) - zbqlix(id43) - c
      if (x.lt.0.0d0) then
       x = x + b
       c = 1.0d0
      else
       c = 0.0d0
      endif
      zbqlix(id43) = x
 
* update array pointers. do explicit check for bounds of each to
* avoid expense of modular arithmetic. if one of them is 0 the others
* won't be
 
      curpos = curpos - 1
      id22 = id22 - 1
      id43 = id43 - 1
      if (curpos.eq.0) then
       curpos=43
      elseif (id22.eq.0) then
       id22 = 43
      elseif (id43.eq.0) then
       id43 = 43
      endif
     
* the integer arithmetic there can yield x=0, which can cause     
* problems in subsequent routines (e.g. zbqlexp). the problem    
* is simply that x is discrete whereas u is supposed to     
* be continuous - hence if x is 0, go back and generate another    
* x and return x/b^2 (etc.), which will be uniform on (0,1/b).     
     
      if (x.lt.binv) then
       b2 = b2*b
       goto 5
      endif

      zbqlu01 = x/b2

      end
* ***********************************


      subroutine DatosGrlsDinamica (NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida) 

        integer getcwd, status
        CHARACTER (LEN = 400) dirname
        CHARACTER (LEN = 400) RutaCarpetaDinamica 
        CHARACTER (LEN = 400) RutaGrl
        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida

        status = getcwd( dirname ) 

        open (unit = 88, file = '.ckout', status = 'old')
      
        read(88,*) RutaCarpetaDinamica

        RutaGrl = trim(dirname)//'/'//RutaCarpetaDinamica
     
        OPEN(1, FILE=trim(RutaGrl)//"/.nm", status = 'old')
        READ(1, *) NombreDinamica 
        
        RutaArchivosSalida = trim(RutaGrl)//'/ArchivosSalida_'
     & //trim(NombreDinamica)        
        
        RutaArchivosEntrada = trim(RutaGrl)//'/ArchivosEntrada_'
     & //trim(NombreDinamica)
        
        CLOSE (88)
        CLOSE (1)
        end subroutine DatosGrlsDinamica


        subroutine RutinaDatosNwchem(RutaScratch, 
     & RutaPermanent, BaseGrlNwchem) 
    
        CHARACTER (LEN = 400) NombreDinamica
        CHARACTER (LEN = 400) RutaArchivosEntrada
        CHARACTER (LEN = 400) RutaArchivosSalida
        CHARACTER (LEN = 400) DatosNwchem

        CHARACTER (LEN = 400) RutaScratch
        CHARACTER (LEN = 400) RutaPermanent
        CHARACTER (LEN = 400) BaseGrlNwchem
        CHARACTER (LEN = 400) NombreCPU

        call DatosGrlsDinamica(NombreDinamica, 
     & RutaArchivosEntrada, RutaArchivosSalida)
        
        call system ("whoami > .quien")

        open(unit = 8, file = '.quien')      

        read(8, *) NombreCPU

        DatosNwchem = trim(RutaArchivosEntrada)
     & // '/DatosNwchem_' // trim(NombreDinamica)
     & // '.dat'
  
        OPEN(9, FILE = trim(DatosNwchem))

        READ(9, *) BaseGrlNwchem

        RutaScratch = '/home'//'/'//trim(NombreCPU)
     & //'/scratch'

        RutaPermanent = '/home'//'/'//trim(NombreCPU)
     & //'/permanent'
              
        CLOSE(8)          
        CLOSE(9)

        end subroutine RutinaDatosNwchem

      k
