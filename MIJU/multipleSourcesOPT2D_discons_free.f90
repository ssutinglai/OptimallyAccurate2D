program multipleSourcesOPT2D

  ! Computation of the synthetic seismograms in the time domain
  ! using the optimally accurate operators.
  ! 2D PSV heterogeneous medium
  ! CPML or Cerjan boundary conditions
  !
  !					originally from 1997.6  N. Takeuchi
  !                                                     2016.5. N. Fuji
  !                                  discon operators : 2016.5. O. Ovcharenko
  !                                         colorbars : 2016.5. K. Okubo
  !  

  implicit none
  character(120) :: filename
  integer :: tmpint
  !character(80), parameter :: modelname = 'homo'
  !character(80), parameter :: vpmodel = './2d_homo.vp'
  !character(80), parameter :: vsmodel = './2d_homo.vs'
  !character(80), parameter :: rhomodel = './2d_homo.rho'
  character(80) :: modelname,vpmodel,vsmodel,rhomodel
  !character(80), parameter :: modelname = 'layered'
  !character(80), parameter :: vpmodel = './2d_start.vp'
  !character(80), parameter :: vsmodel = './2d_start.vs'

  

  ! switch OPT / CONV
  !logical,parameter :: optimise = .false.
  logical :: optimise

  ! switch video
  logical, parameter :: videoornot = .true.

  ! writing strains
  logical, parameter :: writingStrain = .true.
  
  !integer :: nt, nx, nz
  !real :: dt, dx, dz
  ! integer :: isx, isz ! source position
  ! for Ricker wavelets
  !real :: f0, t0
  integer, parameter :: nReceiver = 7
  integer :: nrx(1:nReceiver), nrz(1:nReceiver)
  
  !integer, parameter :: iSourceStart = 2
  !integer, parameter :: iSourceInterval=2
  integer :: iSourceStart,iSourceInterval,nSource
  integer, parameter :: maxnSource = 1
  integer :: iisx(1:maxnSource), iisz(1:maxnSource)
  integer :: iSource, iReceiver



  integer, parameter :: maxnz = 600 
  integer, parameter :: maxnt = 3000
  double precision, parameter :: pi=3.1415926535897932d0 
  double precision, parameter :: ZERO = 0.d0
    
  !
  ! parameters for the gridding
  double precision dt,dx,dz
  ! parameters for the wavefield
  integer nt,nx,nz,it,ist,isx,isz,ix,iz,recl_size
  ! Attention ! nx and nz are modified with absorbing boundaries
  double precision ux(maxnz+1,maxnz+1),uz(maxnz+1,maxnz+1)
  double precision ux1(maxnz+1,maxnz+1),ux2(maxnz+1,maxnz+1)
  double precision uz1(maxnz+1,maxnz+1),uz2(maxnz+1,maxnz+1)
  double precision  e1(maxnz+1,maxnz+1), e2(maxnz+1,maxnz+1)
  double precision  e3(maxnz+1,maxnz+1), e4(maxnz+1,maxnz+1)
  double precision  e5(maxnz+1,maxnz+1), e6(maxnz+1,maxnz+1)
  double precision  e7(maxnz+1,maxnz+1), e8(maxnz+1,maxnz+1)
  double precision e13(maxnz+1,maxnz+1),e14(maxnz+1,maxnz+1)
  double precision e15(maxnz+1,maxnz+1),e16(maxnz+1,maxnz+1)
  double precision e17(maxnz+1,maxnz+1),e18(maxnz+1,maxnz+1)
  double precision e19(maxnz+1,maxnz+1),e20(maxnz+1,maxnz+1)
  double precision  f1(maxnz+1,maxnz+1), f2(maxnz+1,maxnz+1)
  double precision  f3(maxnz+1,maxnz+1), f4(maxnz+1,maxnz+1)
  double precision  f5(maxnz+1,maxnz+1), f6(maxnz+1,maxnz+1)
  double precision  f7(maxnz+1,maxnz+1), f8(maxnz+1,maxnz+1)
  double precision f13(maxnz+1,maxnz+1),f14(maxnz+1,maxnz+1)
  double precision f15(maxnz+1,maxnz+1),f16(maxnz+1,maxnz+1)
  double precision f17(maxnz+1,maxnz+1),f18(maxnz+1,maxnz+1)
  double precision f19(maxnz+1,maxnz+1),f20(maxnz+1,maxnz+1)
  double precision work(maxnz+1,32)

  ! for discontinuities

  double precision, dimension(maxnz+1,maxnz+1) :: ee12,ee34,ee56,ee65,ee78,ee87
  double precision, dimension(maxnz+1,maxnz+1) :: ff12,ff34,ff56,ff65,ff78,ff87

  
  ! parameter for the structure
  !double precision rrho(6),llam(6),mmu(6)
  character(80) :: vpfile, vsfile, rhofile   ! ,modelname
  double precision :: rho(maxnz+1,maxnz+1)
  double precision :: lam(maxnz+1,maxnz+1),mu(maxnz+1,maxnz+1)
  double precision :: fx(maxnz+1,maxnz+1),fz(maxnz+1,maxnz+1)
  double precision :: vs(maxnz+1,maxnz+1),vp(maxnz+1,maxnz+1)
  double precision :: cp ! maxvalue of vp
  
  double precision Courant_number
  ! parameter for the receiver
  !integer :: nReceiver ! number of receiver
  integer, parameter :: maxReceiver = nReceiver
  integer :: ir,j
  !integer :: nrx(1:maxReceiver),nrz(1:maxReceiver)
  real :: synx(0:maxnt,1:maxReceiver),synz(0:maxnt,1:maxReceiver),time(0:maxnt)
  real :: video(maxnz+1,maxnz+1)
  character(150) :: outfile
 
  
  ! parameter for the waveform
  double precision t
  !parameter for video
  real,dimension(maxnz+1,maxnz+1) :: snapux,snapuz
  integer, parameter :: IT_DISPLAY = 10
 
  integer(2) head(1:120)
  character(80) :: routine
  
  logical,parameter :: dummylog = .false.

  ! switch C-PML 
  logical, parameter :: USE_PML_XMIN = .true.
  logical, parameter :: USE_PML_XMAX = .true.
  logical, parameter :: USE_PML_YMIN = .true.
  logical, parameter :: USE_PML_YMAX = .true.
  ! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 100
  double precision, parameter :: CerjanRate = 0.0015
  double precision :: weightBC(maxnz+1,maxnz+1)
  ! Cerjan boundary condition
  integer :: lmargin(1:2),rmargin(1:2)
  

  ! Ricker wavelets source
  double precision f0,t0
  !double precision tp,ts


  ! for evolution of total energy in the medium
  double precision epsilon_xx,epsilon_yy,epsilon_xy
  double precision, dimension(maxnt+1) :: total_energy_kinetic,total_energy_potential
  
  ! power to compute d0 profile
  double precision, parameter :: NPOWER = 2.d0

  double precision, parameter :: K_MAX_PML = 1.d0 ! from Gedney page 8.11
  double precision :: ALPHA_MAX_PML

  ! for water
  
  integer :: liquidmarkers(maxnz+1,maxnz+1)

  ! for discontinuities

  integer :: markers(maxnz+1,maxnz+1)
  integer :: nDiscon ! number of discontinuities
  integer :: lengthDiscon ! with x,z coordinates
  double precision, allocatable :: dscr(:,:,:) ! discontinuity coordinates
  double precision :: tmpvaluex,tmpvaluez

  ! for free surface
  integer :: zerodisplacement(maxnz+1,maxnz+1)
  integer :: lengthFreeSurface ! with x,z coordinates
  double precision, allocatable :: free(:,:)
 
  

  ! for waveform inversion
  
  real(kind(0e0)) :: singleStrainDiagonal(maxnz+1,maxnz+1)
  real(kind(0e0)), allocatable:: tmpsingleStrain(:,:)

  
  character(140) :: commandline
  
  




  ! Reading Inf File
110 format(a80)
  read(5,110) modelname
  read(5,110) vpmodel
  read(5,110) vsmodel
  read(5,110) rhomodel
  read(5,'(L1)') optimise
  read(5,*) iSourceStart,iSourceInterval,nSource
 
  call system('mkdir ./inffile')
   
  commandline="mkdir synthetics"
  call system(commandline)
  commandline="mkdir snapshots"
  call system(commandline)
  commandline="mkdir videos"
  call system(commandline)
  commandline="mkdir synthetics/"//trim(modelname)
  call system(commandline)
  commandline="mkdir videos/"//trim(modelname)
  call system(commandline)
  commandline="mkdir strains"
  call system(commandline)
  commandline="mkdir strains/"//trim(modelname)
  call system(commandline)


  vpfile=vpmodel
  vsfile=vsmodel
  rhofile=rhomodel

  nt=2000
  nx=399
  nz=199
  dt=2.d-3
  dx=2.d-2
  dz=2.d-2
  
  f0=6.5d0
  t0=2.d-1

  allocate(tmpsingleStrain(1:nx+1,1:nz+1))




  ! Discontinuity configuration


  ! diagonal discontinuity

  if(1.eq.1) then
  
  nDiscon = 1
  lengthDiscon = 40*nx+1
  
  if(nDiscon.ne.0) then
     allocate(dscr(1:2,1:lengthDiscon,1:nDiscon))
     do ix =1,lengthDiscon
        dscr(1,ix,1) = dble(ix-1)*dx/40.d0
        dscr(2,ix,1) = 199.d0*dz-dscr(1,ix,1)*199.d0/399.d0
     enddo
  endif
  markers(1:maxnz,1:maxnz) = 0
  markers(1:nx+1,1:nz+1) = 1 ! for the moment NF will search for all the points (of course it is not good)

  markers(1:maxnz,1:maxnz) = 0
  do ix = 1,nx+1
     tmpint=nint(199.d0-dble(ix-1)*199.d0/399.d0)
     if(tmpint-3.ge.1) then
        markers(ix,tmpint-3:tmpint+3) = 1
     else if(tmpint-2.ge.1) then
        markers(ix,tmpint-2:tmpint+3) = 1
     else
        markers(ix,1:tmpint+3) = 1
     endif
  enddo
 
  endif

  ! Circle discontinuity
  
  if(0.eq.1) then
     
     nDiscon = 2
     lengthDiscon = 40*100+1
     
     if(nDiscon.ne.0) then
        allocate(dscr(1:2,1:lengthDiscon,1:nDiscon))
        do ix =1,lengthDiscon
           dscr(1,ix,1) = 99.d0*dx+dble(ix-1)*dx/40.d0
           dscr(1,ix,2) = dscr(1,ix,1)
           !dscr(2,ix,1) = 99.d0+sqrt(2500.d0*dx**2-(dble(ix-1)*dx/40.d0-50.d0)**2)
           dscr(2,ix,1) = 99.d0
           dscr(2,ix,2) = 99.d0-sqrt(2500.d0*dx**2-(dble(ix-1)*dx/40.d0-50.d0)**2)
        enddo
     endif
     markers(1:nx+1,1:nz+1)=1
  endif


  ! Oleg discontinuties

  if(0.eq.0) then

     open(1, file=

     stop


  endif
  
  ! Receiver position


  do iReceiver = 1, nReceiver
     nrx(iReceiver)=70+30*iReceiver
     nrz(iReceiver)=100
  enddo
  
  do iSource = 1, nSource
     iisx(iSource)=iSourceStart+iSourceInterval*(iSource-1)
     iisz(iSource)=10
     write(filename, '(I5,".",I5,".inf")') iisx(iSource),iisz(iSource)
     do j=1, 12
        if(filename(j:j).eq.' ') filename(j:j)='0'
     enddo
    filename= './inffile/'//trim(modelname)//'.'//filename
    
       

     open(1, file=filename, form='formatted')
     write(1,'(a)') 'c grids'
     write(1,*) nt, nx, nz
     write(1,*) dt, dx, dz
     write(1,'(a)') 'c modelname'
     write(1,'(a)') trim(modelname)
     write(1,'(a)') trim(vpmodel)
     write(1,'(a)') trim(vsmodel)
     write(1,'(a)') trim(rhomodel)

     write(1,'(a)') 'c source position (in grids) '
     write(1,*) iisx(iSource),iisz(iSource)
     write(1,'(a)') 'c source time function (Ricker wavelet)'
     write(1,*) f0,t0     
     write(1,'(a)') 'c receivers information'
     write(1,*) nReceiver
     do iReceiver = 1, nReceiver
        write(1,*) nrx(iReceiver), nrz(iReceiver)
     enddo
     

     !write(1,*) 'c'
     write(1,'(a)') 'end'     
  enddo


  !!!! for each source we calculate synthetics

  
  !computing the intermediate parameters
  
  call calstruct( maxnz,rhofile,dx,dz,nx,nz,rho )
  call calstruct( maxnz,vpfile,dx,dz,nx,nz,vp)
  call calstruct( maxnz,vsfile,dx,dz,nx,nz,vs )
    
    

  ! Free surface configuration

  lengthFreeSurface = 0
  if(lengthFreeSurface.ne.0) then
     allocate(free(1:2,1:lengthFreeSurface))
     do ix=1,lengthFreeSurface
        free(1,ix) = dble(ix-1)*dx/40.d0
        free(2,ix) = 5.d-1+3.d-1*sin(free(1,ix)/(dble(nx)*dx)*pi*4)
     enddo
     zerodisplacement(1:maxnz,1:maxnz)=0
     do ix=1,nx+1
        tmpint=nint((5.d-1+3.d-1*sin((dble(ix-1)*dx)/(dble(nx)*dx)*pi*4))/dx)
        zerodisplacement(ix,1:tmpint) = 1
        vp(ix,1:tmpint)=0.d0
        vs(ix,1:tmpint)=0.d0
     enddo

     
  endif
  

  

  lam = 0.d0
  mu = 0.d0
  !call datainit( maxnz,maxnz,lam )
  !call datainit( maxnz,maxnz,mu )
     
  ! Cerjan boundary
  
  lmargin(1)=NPOINTS_PML
  rmargin(1)=NPOINTS_PML
  lmargin(2)=NPOINTS_PML
  rmargin(2)=NPOINTS_PML
  
  liquidmarkers = 0


  call calstruct2(maxnz,nx,nz,rho,vp,vs,lam,mu,liquidmarkers)
  
  !print *, "OK ??"
  !print *, maxnz,nx,nz,markers,lmargin,rmargin
  !print *, rho
  
  write(12,*) rho
  write(13,*) vp
  write(14,*) vs
  !stop

  call calstructBC(maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)
  

  ! Smoothed version of CONV/OPT operators

  call cales( maxnz,nx,nz,rho,lam,mu,dt,dx,dz, &
       e1, e2, e3, e4, e5, e6, e7, e8, &
       e13,e14,e15,e16,e17,e18,e19,e20, &
       f1, f2, f3, f4, f5, f6, f7, f8, &
       f13,f14,f15,f16,f17,f18,f19,f20 )

  ! discontinuities
  
  ee12 = 0.d0
  ee34 = 0.d0
  ee56 = 0.d0
  ee65 = 0.d0
  ee78 = 0.d0
  ee87 = 0.d0

  ff12 = 0.d0
  ff34 = 0.d0
  ff56 = 0.d0
  ff65 = 0.d0
  ff78 = 0.d0
  ff87 = 0.d0
  
  if(nDiscon.ne.0) then
 
     ! changing dscr by putting lmargin(1) and (2)
     tmpvaluex=dble(lmargin(1))*dx
     tmpvaluez=dble(lmargin(2))*dz
     do ix=1,nDiscon
        do iz=1,lengthDiscon
           dscr(1,iz,ix)=dscr(1,iz,ix)+tmpvaluex
           dscr(2,iz,ix)=dscr(2,iz,ix)+tmpvaluez
        enddo
     enddo
          
     call cales_discon( maxnz,nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, & 
     ! hereafter are new variables for cales_discon
     ee12,ee34,ee56,ee65,ee78,ee87, &
     ff12,ff34,ff56,ff65,ff78,ff87, &
     markers,nDiscon,lengthDiscon,dscr)

  endif
  
  if(lengthFreeSurface.ne.0) then
          
     call cales_free( maxnz,nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, & 
     ! hereafter are new variables for cales_discon
     ee12,ee34,ee56,ee65,ee78,ee87, &
     ff12,ff34,ff56,ff65,ff78,ff87, &
     zerodisplacement,lengthFreeSurface,free)

     
  endif



  if(lengthFreeSurface.ne.0) then
     ! changing free by putting lmargin(1) and (2)
     tmpvaluex=dble(lmargin(1))*dx
     tmpvaluez=dble(lmargin(2))*dz
     
     do ix=1,lengthFreeSurface
           free(1,ix)=free(1,ix)+tmpvaluex
           free(2,ix)=free(2,ix)+tmpvaluez           
     enddo
  endif

  ! for Cerjan absorbing boundary


  weightBC=1.d0
     
  call compNRBCpre(weightBC(1:nx+1,1:nz+1),CerjanRate,lmargin,rmargin,nx+1,nz+1)
  

  do ir= 1, nReceiver
     nrx(ir)=nrx(ir)+lmargin(1)
     nrz(ir)=nrz(ir)+lmargin(2)
  enddo
  

  do iSource = 1, nSource
     !iisx(iSource)=iSourceStart+iSourceInterval*(iSource-1)
     !iisz(iSource)=1
     isx=iisx(iSource)
     isz=iisz(iSource)


     ! for video (without boundary)
     recl_size=(nx+1)*(nz+1)*kind(0e0)
    
     
     ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0) ! from Festa and Vilotte
     
     ! Initializing the data
     call datainit( maxnz,maxnz,ux )
     call datainit( maxnz,maxnz,uz )
     call datainit( maxnz,maxnz,ux1 )
     call datainit( maxnz,maxnz,uz1 )
     call datainit( maxnz,maxnz,ux2 )
     call datainit( maxnz,maxnz,uz2 )
     

     call datainit( maxnz,31,work )
     
     
 
     ! R. Courant et K. O. Friedrichs et H. Lewy (1928)
     cp=maxval(vp)
     Courant_number = cp * dt * sqrt(1.d0/dx**2 + 1.d0/dz**2)
     print *, 'Courant number is', Courant_number
     
     


     
     call datainit( maxnz,maxnz,fx)
     call datainit( maxnz,maxnz,fz)
     

     ! ist = dnint( 2 * tp / dt )
     ! isx = nx / 2 + 1
     ! isz = nz / 2 + 1
     
     
     ist=nt/4
     
     ! isx = 30s
     ! isz = 4
     
     
     isx=isx+lmargin(1)
     isz=isz+lmargin(2)
     
 
     
     t=0.d0
     time(0)=t
     do it=0,nt
        
        call calf2( maxnz,it,t,ist,isx,isz,dt,dx,dz,rho(isx,isz),f0,t0,fx,fz )
        t=t+dt
        !write(13,*) t, fx(isx,isz),fz(isx,isz)
        
     enddo
     !print *, maxnz,it,t,ist,isx,isz,dt,dx,dz,rho(isx,isz),f0,t0
     !stop

     
     t = 0.d0
     !write(14,*) real(t),real(ux(nrx,nrz)),real(uz(nrx,nrz))
     do ir = 1,nReceiver
        synx(0,ir)=ux(nrx(ir),nrz(ir))
        synz(0,ir)=uz(nrx(ir),nrz(ir))
     enddo

     

     do it=0,nt
        call calf2( maxnz,it,t,ist,isx,isz,dt,dx,dz,rho(isx,isz),f0,t0,fx,fz )
        ! evaluating the next step
        
        !if(nDiscon.eq.0) then
        !   call calstep( maxnz,nx,nz, &
        !        e1, e2, e3, e4, e5, e6, e7, e8, &
        !        e13,e14,e15,e16,e17,e18,e19,e20, &
        !        f1, f2, f3, f4, f5, f6, f7, f8, &
        !        f13,f14,f15,f16,f17,f18,f19,f20, &
        !        ux,uz,ux1,ux2,uz1,uz2,isx,isz,fx,fz, &
        !        work(1,1), work(1,5), work(1,9),work(1,13), &
        !        work(1,17),work(1,18),work(1,20),work(1,21), &
        !        work(1,23),work(1,24),work(1,28),work(1,29), optimise)
           
        !else
           call calstep_discon( maxnz,nx,nz, &
                e1, e2, e3, e4, e5, e6, e7, e8, &
                e13,e14,e15,e16,e17,e18,e19,e20, &
                f1, f2, f3, f4, f5, f6, f7, f8, &
                f13,f14,f15,f16,f17,f18,f19,f20, &
                ux,uz,ux1,ux2,uz1,uz2,isx,isz,fx,fz, &
                work(1,1), work(1,5), work(1,9),work(1,13), &
                work(1,17),work(1,18),work(1,20),work(1,21), &
                work(1,23),work(1,24),work(1,28),work(1,29), optimise, & 
                ! Hereafter are new variables for cales_discon
                ee12,ee34,ee56,ee65,ee78,ee87, &
                ff12,ff34,ff56,ff65,ff78,ff87)
           

        !endif


           if(lengthFreeSurface.ne.0) then
              do iz=1,nz+1
                 do ix=1,nx+1
                    if(zerodisplacement(ix,iz).eq.1) then
                       ux(ix,iz) = 0.d0
                       uz(ix,iz) = 0.d0
                    endif
                 enddo
              enddo
           endif


           ! increment of t
        t = t + dt
        time(it)=t
        !write(14,*) real(t),real(ux(nrx,nrz)),real(uz(nrx,nrz))
        do ir = 1,nReceiver
           synx(it,ir)=ux(nrx(ir),nrz(ir))
           synz(it,ir)=uz(nrx(ir),nrz(ir))
        enddo
        
     
        
        ! applying Cerjan boundary
        
        do iz=1,nz+1
           do ix=1,nx+1
              uz(ix,iz)=uz(ix,iz)*weightBC(ix,iz)
              ux1(ix,iz)=ux1(ix,iz)*weightBC(ix,iz)
              uz1(ix,iz)=uz1(ix,iz)*weightBC(ix,iz)
              ux2(ix,iz)=ux2(ix,iz)*weightBC(ix,iz)
              uz2(ix,iz)=uz2(ix,iz)*weightBC(ix,iz)
           enddo
        enddo
        
        
        
        ! calculating strains
        
        
        if(writingStrain) then
           singleStrainDiagonal=0.e0
           tmpsingleStrain=0.e0
           call calStrainDiagonal(maxnz,nx,nz,ux,uz,lmargin,rmargin,singleStrainDiagonal)
          


           if(optimise) then
              write(outfile,'("strain",I5,".",I5,".",I5,".OPT_dat") ') it,isx-lmargin(1),isz-lmargin(2)
           else
              write(outfile,'("strain",I5,".",I5,".",I5,".CON_dat") ') it,isx-lmargin(1),isz-lmargin(2)
           endif
           do j=1,24
              if(outfile(j:j).eq.' ') outfile(j:j)='0'
           enddo
           
           outfile = './strains/'//trim(modelname)//'/'//outfile
           open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)

           tmpsingleStrain(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2)) = &
                singleStrainDiagonal(lmargin(1)+1:nx+1-rmargin(1),lmargin(2)+1:nz+1-rmargin(2))
           write(1,rec=1)  tmpsingleStrain
           close(1,status='keep')
        endif


        
        !write(*,*) it, ' of ', nt
        if(mod(it,IT_DISPLAY) == 0)then
           !
           !head=0
           !head(58) = nx
           !head(59) = dz * 1E3
           !snapux=0.e0
           !snapuz=0.e0
           !snapux(1:nx,1:nz) = ux(1:nx,1:nz)
           !snapuz(1:nx,1:nz) = uz(1:nx,1:nz)
           !write(routine,'(a12,i5.5,a9)') './snapshots/',it,'snapUx.su'
           !open(21,file=routine,access='stream')
           !do j = 1,nx,1
           !   write(21) head,(real(snapux(k,j)),k=1,nz)
           !enddo
           !close(21)
           !write(routine,'(a12,i5.5,a9)') './snapshots/',it,'snapUz.su'
           !open(21,file=routine,access='stream')
           !do j = 1,nx,1
           !   write(21) head,(real(snapuz(k,j)),k=1,nz)
           !enddo
           !close(21)
           
           
           
           !call create_color_image(ux(1:nx+1,1:nz+1),nx+1,nz+1,it,isx,isz,ix_rec,iz_rec,1,0, &
           !     dummylog,dummylog,dummylog,dummylog,1)
           !call create_color_image(ux(1:nx+1,1:nz+1),nx+1,nz+1,it,isx,isz,ix_rec,iz_rec,1,&
           !    NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,1)
           if(videoornot) then
              call create_color_image(uz(1:nx+1,1:nz+1),nx+1,nz+1,it,isx,isz, &
                   nrx(1:nReceiver),nrz(1:nReceiver),nReceiver, &
                   NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,2)
              
              
           endif
           
           !if(optimise) then
           !   write(outfile,'("video",I5,".",I5,".",I5,".OPT_UX") ') it,isx-lmargin(1),isz-lmargin(2)
           !else
           !   write(outfile,'("video",I5,".",I5,".",I5,".CON_UX") ') it,isx-lmargin(1),isz-lmargin(2)
           !endif
           !do j=1,24
           !   if(outfile(j:j).eq.' ') outfile(j:j)='0'
           !enddo
           
           !outfile = './synthetics/'//trim(modelname)//'/'//outfile
           !open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
           !video(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2))= &
           !     ux(lmargin(1)+1:nx+1-rmargin(1),lmargin(2)+1:nz+1-rmargin(2))
           !write(1,rec=1)  video(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2))
           !close(1,status='keep')
           
           
           
           !if(optimise) then
           !   write(outfile,'("video",I5,".",I5,".",I5,".OPT_UX") ') it,isx-lmargin(1),isz-lmargin(2)
           !else
           !   write(outfile,'("video",I5,".",I5,".",I5,".CON_UX") ') it,isx-lmargin(1),isz-lmargin(2)
           !endif
           !do j=1,24
           !   if(outfile(j:j).eq.' ') outfile(j:j)='0'
           !enddo
           
           !outfile = './synthetics/'//trim(modelname)//'/'//outfile
           !open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
           !video(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2))= &
           !     ux(lmargin(1)+1:nx+1-rmargin(1),lmargin(2)+1:nz+1-rmargin(2))
           !write(1,rec=1)  video(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2))
           !close(1,status='keep')
           
        endif

        
        
        !call compNRBC2(ux(1:nx+1,1:nz+1),ux1(1:nx+1,1:nz+1),ux2(1:nx+1,1:nz+1), &
        !     uz(1:nx+1,1:nz+1),uz1(1:nx+1,1:nz+1),uz2(1:nx+1,1:nz+1), CerjanRate, lmargin, rmargin,nx+1,nz+1)
        
     enddo

     !write(18,*) singleStrainDiagonal(:,:)


     if(videoornot) then
        
        if(optimise) then
           write(outfile,'("video",".",I5,".",I5,".OPT.mp4") ') isx-lmargin(1),isz-lmargin(2)
        else
           write(outfile,'("video",".",I5,".",I5,".CON.mp4") ') isx-lmargin(1),isz-lmargin(2)
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './videos/'//trim(modelname)//'/'//outfile
        
        
        commandline="ffmpeg -framerate 5 -pattern_type glob -i 'snapshots/*.png' -c:v libx264 -pix_fmt yuv420p "//outfile
     
     endif
     
     
     call system(commandline)
  
  
     do ir = 1,nReceiver
        if(optimise) then
           write(outfile,'(I5,".",I5,".",I5,".",I5,".OPT_UX") ') nrx(ir)-lmargin(1),nrz(ir)-lmargin(2), &
                isx-lmargin(1),isz-lmargin(2)
        else
           write(outfile,'(I5,".",I5,".",I5,".",I5,".CON_UX") ') nrx(ir)-lmargin(1),nrz(ir)-lmargin(2), &
                isx-lmargin(1),isz-lmargin(2)
        endif
        
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './synthetics/'//trim(modelname)//'/'//outfile
        open(1, file=outfile,status='unknown',form='formatted')
        do it=0,nt
           write (1,*) time(it),synx(it,ir)
        enddo
        close(1)
        !open(1,file=outfile,form='unformatted',access='direct',recl=kind(0e0)*(nt+1))
        !write(1,rec=1) synx(0:nt,ir)
        !close(1,status='keep')

        
        if(optimise) then
           write(outfile,'(I5,".",I5,".",I5,".",I5,".OPT_UZ") ') nrx(ir)-lmargin(1),nrz(ir)-lmargin(2), &
                isx-lmargin(1),isz-lmargin(2)
        else
           write(outfile,'(I5,".",I5,".",I5,".",I5,".CON_UZ") ') nrx(ir)-lmargin(1),nrz(ir)-lmargin(2), &
                isx-lmargin(1),isz-lmargin(2)
        endif
        
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './synthetics/'//trim(modelname)//'/'//outfile
        open(1, file=outfile,status='unknown',form='formatted')
        do it=0,nt
           write (1,*) time(it), synz(it,ir)
        enddo
        close(1)

        !open(1,file=outfile,form='unformatted',access='direct',recl=kind(0e0)*(nt+1))
        !write(1,rec=1) synz(0:nt,ir)
        !close(1,status='keep')


        
     enddo
     
  enddo




end program multipleSourcesOPT2D




subroutine datainit( nx,nz,ux )
  
  integer nx,nz
  double precision ux(nx+1,*)
  integer i,j
  
  do j=1,nz+1
     do i=1,nx+1
        ux(i,j) = 0.d0
     enddo
  enddo
  
  return
end subroutine datainit


subroutine calstruct( maxnz,file2d,dx,dz,nx,nz,rho )
  implicit none
  integer maxnz,nx,nz
  double precision dx,dz,rho(1:maxnz+1,1:maxnz+1)
  real(kind(1.e0)),allocatable :: rrho(:,:)
  integer i,j,k,nox(6),noz(6)
  double precision x,z,xmax,zmax,trho,coef1,coef2
  integer recl_size
  character*80 file2d
  recl_size=kind(1.0)*(nx+1)*(nz+1)
  
  allocate(rrho(1:nx+1,1:nz+1))
  open (1,file=file2d,form='unformatted',access='direct',recl=recl_size)
  read(1,rec=1) rrho(1:nx+1,1:nz+1)
  close(1)
  rho(1:nx+1,1:nz+1)=1.d-3*rrho(1:nx+1,1:nz+1)
  
  deallocate(rrho)
  return
end subroutine calstruct

subroutine calstruct2(maxnz,nx,nz,rho,vp,vs,lam,mu,liquidmarkers)
  implicit none
  
  integer i,j,maxnz,nx,nz
  double precision rho(maxnz+1,maxnz+1),vp(maxnz+1,maxnz+1),vs(maxnz+1,maxnz+1)
  double precision lam(maxnz+1,maxnz+1),mu(maxnz+1,maxnz+1)
  integer liquidmarkers(maxnz+1,maxnz+1)

  do i=1,nx+1
     do j=1,nz+1
        if(vs(i,j).eq.0.d0) then
           liquidmarkers(i,j)=1
           !NF should take out this now
           !vs(i,j)=vp(i,j)/1.7d0
        endif
           
        mu(i,j)=rho(i,j)*vs(i,j)*vs(i,j)
        lam(i,j)=rho(i,j)*vp(i,j)*vp(i,j)-2*mu(i,j)

        
     enddo
  enddo
end subroutine calstruct2
  


subroutine calstructBC(maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)
  implicit none
  integer :: i,j,maxnz,nx,nz,nnx,nnz
  double precision ::lam(maxnz+1,maxnz+1),mu(maxnz+1,maxnz+1),rho(maxnz+1,maxnz+1)  
  integer :: rmargin(1:2), lmargin(1:2)
  integer :: markers(maxnz+1,maxnz+1) ! discontinuities
  integer :: zerodisplacement(maxnz+1,maxnz+1) ! above free surface
  integer :: liquidmarkers(maxnz+1,maxnz+1)
  ! real(kind(0d0)), dimension(maxnz+1,maxnz+1) ::mmu,rrho,llam
  double precision, allocatable :: mmu(:,:), rrho(:,:), llam(:,:)
  integer, allocatable :: mmarkers(:,:),lliquidmarkers(:,:)
  integer, allocatable :: zzerodisplacement(:,:)
  
  allocate(rrho(1:maxnz+1,1:maxnz+1))
  allocate(mmu(1:maxnz+1,1:maxnz+1))
  allocate(llam(1:maxnz+1,1:maxnz+1))
  allocate(mmarkers(1:maxnz+1,1:maxnz+1))
  allocate(lliquidmarkers(1:maxnz+1,1:maxnz+1))
  allocate(zzerodisplacement(1:maxnz+1,1:maxnz+1))

  mmu=0.d0
  rrho=0.d0
  mmarkers=0
  lliquidmarkers=0
  llam=0.d0
  zzerodisplacement=0


  llam(1+lmargin(1):nx+1+lmargin(1),1+lmargin(2):nz+1+lmargin(2))=lam(1:nx+1,1:nz+1)
  mmu(1+lmargin(1):nx+1+lmargin(1),1+lmargin(2):nz+1+lmargin(2))=mu(1:nx+1,1:nz+1)
  rrho(1+lmargin(1):nx+1+lmargin(1),1+lmargin(2):nz+1+lmargin(2))=rho(1:nx+1,1:nz+1)

  mmarkers(1+lmargin(1):nx+1+lmargin(1),1+lmargin(2):nz+1+lmargin(2))=markers(1:nx+1,1:nz+1)
  lliquidmarkers(1+lmargin(1):nx+1+lmargin(1),1+lmargin(2):nz+1+lmargin(2))= &
       liquidmarkers(1:nx+1,1:nz+1)
  zzerodisplacement(1+lmargin(1):nx+1+lmargin(1),1+lmargin(2):nz+1+lmargin(2))= &
       zerodisplacement(1:nx+1,1:nz+1)


  ! 4 corners

  llam(1:lmargin(1),1:lmargin(2))=lam(1,1)
  mmu(1:lmargin(1),1:lmargin(2))=mu(1,1)
  rrho(1:lmargin(1),1:lmargin(2))=rho(1,1)
  zzerodisplacement(1:lmargin(1),1:lmargin(2))=zerodisplacement(1,1)


  llam(1:lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2))=lam(1,nz+1)
  mmu(1:lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2))=mu(1,nz+1)
  rrho(1:lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2))=rho(1,nz+1)
  zzerodisplacement(1:lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2)) &
       = zerodisplacement(1,nz+1)
  !print *, llam(1,nz+lmargin(2)+5),mu(1,nz+1),rho(1,nz+1)

  llam(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1:lmargin(2))=lam(nx+1,1)
  mmu(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1:lmargin(2))=mu(nx+1,1)
  rrho(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1:lmargin(2))=rho(nx+1,1)
  zzerodisplacement(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1:lmargin(2)) &
       = zerodisplacement(nx+1,1)

  llam(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2)) &
       = lam(nx+1,nz+1)
  mmu(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2)) &
       = mu(nx+1,nz+1)
  rrho(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2)) &
       = rho(nx+1,nz+1)
  zzerodisplacement(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2)) &
       = zerodisplacement(nx+1,nz+1)
  

  ! 4 rectangles

  do i = 1,lmargin(1)
     llam(i,1+lmargin(2):nz+1+lmargin(2)) = lam(1,1:nz+1)
     mmu(i,1+lmargin(2):nz+1+lmargin(2)) = mu(1,1:nz+1)
     rrho(i,1+lmargin(2):nz+1+lmargin(2)) = rho(1,1:nz+1)
     zzerodisplacement(i,1+lmargin(2):nz+1+lmargin(2))  &
          = zerodisplacement(1,1:nz+1)
  enddo
  
  do i = 1+nx+1+lmargin(1),rmargin(1)+nx+1+lmargin(1)
     llam(i,1+lmargin(2):nz+1+lmargin(2)) = lam(nx+1,1:nz+1)
     mmu(i,1+lmargin(2):nz+1+lmargin(2)) = mu(nx+1,1:nz+1)
     rrho(i,1+lmargin(2):nz+1+lmargin(2)) = rho(nx+1,1:nz+1)
     zzerodisplacement(i,1+lmargin(2):nz+1+lmargin(2)) &
          = zerodisplacement(nx+1,1:nz+1)
  enddo
  
  do i = 1,lmargin(2)
     llam(1+lmargin(1):nx+1+lmargin(2),i)=lam(1:nx+1,1)
     mmu(1+lmargin(1):nx+1+lmargin(2),i)=mu(1:nx+1,1)
     rrho(1+lmargin(1):nx+1+lmargin(2),i)=rho(1:nx+1,1)
     zzerodisplacement(1+lmargin(1):nx+1+lmargin(2),i) &
          =zerodisplacement(1:nx+1,1)
  enddo

  do i = 1+nz+1+lmargin(2),rmargin(2)+nz+1+lmargin(2)
     llam(1+lmargin(1):nx+1+lmargin(2),i) = lam(1:nx+1,nz+1)
     mmu(1+lmargin(1):nx+1+lmargin(2),i) = mu(1:nx+1,nz+1)
     rrho(1+lmargin(1):nx+1+lmargin(2),i) = rho(1:nx+1,nz+1)
     zzerodisplacement(1+lmargin(1):nx+1+lmargin(2),i) &
          = zerodisplacement(1:nx+1,nz+1)

  enddo

  nnx=rmargin(1)+nx+lmargin(1)
  nnz=rmargin(2)+nz+lmargin(2)

  nx=nnx
  nz=nnz
  lam=0.d0
  rho=0.d0
  mu=0.d0
  liquidmarkers = 0
  zerodisplacement = 0
  markers = 0

  lam(1:nx+1,1:nz+1) = llam(1:nx+1,1:nz+1)
  rho(1:nx+1,1:nz+1) = rrho(1:nx+1,1:nz+1)
  mu(1:nx+1,1:nz+1) = mmu(1:nx+1,1:nz+1)
  markers(1:nx+1,1:nz+1)=mmarkers(1:nx+1,1:nz+1)
  zerodisplacement(1:nx+1,1:nz+1)=zzerodisplacement(1:nx+1,1:nz+1)
  liquidmarkers(1:nx+1,1:nz+1)=lliquidmarkers(1:nx+1,1:nz+1)
  !print *, nx,nz
  !write(12,*) rho(:,:)
  !write(13,*) lam(:,:)
  !write(14,*) mmu(1:nx+1,1:nz+1)
  !stop
  
  deallocate(llam)
  deallocate(rrho)
  deallocate(mmu)
  deallocate(mmarkers)
  deallocate(lliquidmarkers)
  deallocate(zzerodisplacement)
end subroutine calstructBC


 





subroutine cales( maxnz,nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20 )
  
  integer maxnz,nx,nz
  double precision rho(maxnz+1,*),lam(maxnz+1,*),mu(maxnz+1,*)
  double precision dt,dx,dz
  double precision  e1(maxnz+1,*), e2(maxnz+1,*), e3(maxnz+1,*)
  double precision  e4(maxnz+1,*), e5(maxnz+1,*), e6(maxnz+1,*)
  double precision  e7(maxnz+1,*), e8(maxnz+1,*)
  double precision e13(maxnz+1,*),e14(maxnz+1,*),e15(maxnz+1,*)
  double precision e16(maxnz+1,*),e17(maxnz+1,*),e18(maxnz+1,*)
  double precision e19(maxnz+1,*),e20(maxnz+1,*)
  double precision  f1(maxnz+1,*), f2(maxnz+1,*), f3(maxnz+1,*)
  double precision  f4(maxnz+1,*), f5(maxnz+1,*), f6(maxnz+1,*)
  double precision  f7(maxnz+1,*), f8(maxnz+1,*)
  double precision f13(maxnz+1,*),f14(maxnz+1,*),f15(maxnz+1,*)
  double precision f16(maxnz+1,*),f17(maxnz+1,*),f18(maxnz+1,*)
  double precision f19(maxnz+1,*),f20(maxnz+1,*)
  integer ix,iz
  double precision dt2,dx2,dz2,dxdz

  ! for smoothed part of the model :
  !  we use Zahradnik operators and optimally accurate operators

  
  dt2 = dt * dt
  dx2 = dx * dx
  dz2 = dz * dz
  dxdz = dx * dz
  
  do iz=2,nz
     do ix=2,nx
        e1(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix-1,iz) + lam(ix,iz) ) &
             + 2.d0 * ( mu(ix-1,iz) + mu(ix,iz) ) ) &
             / ( 2.d0 * dx2 )
        e2(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix,iz) + lam(ix+1,iz) ) &
             + 2.d0 * ( mu(ix,iz) + mu(ix+1,iz) ) ) &
             / ( 2.d0 * dx2 )
        e3(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix,iz-1) + mu(ix,iz) ) &
             / ( 2.d0 * dz2 )
        e4(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix,iz) + mu(ix,iz+1) ) &
             / ( 2.d0 * dz2 )
        e5(ix,iz) = dt2 / rho(ix,iz) * lam(ix-1,iz) &
             / ( 4.d0 * dxdz )
        e6(ix,iz) = dt2 / rho(ix,iz) * lam(ix+1,iz) &
             / ( 4.d0 * dxdz )
        e7(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz-1) &
             / ( 4.d0 * dxdz )
        e8(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+1) &
             / ( 4.d0 * dxdz )
        e13(ix,iz) = dt2 / rho(ix,iz) * lam(ix-1,iz) &
             * ( -5.d0 ) / ( 1728.d0 * dxdz )
        e14(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
             * ( -3.d0 ) / ( 1728.d0 * dxdz )
        e15(ix,iz) = dt2 / rho(ix,iz) * lam(ix+1,iz) &
           * ( +9.d0 ) / ( 1728.d0 * dxdz )
        if ( ix+2.le.nx+1 ) then
           e16(ix,iz) = dt2 / rho(ix,iz) * lam(ix+2,iz) &
                * ( -1.d0 ) / ( 1728.d0 * dxdz )
        else
           e16(ix,iz) = 0.d0
        endif
        e17(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz-1) &
             * ( -5.d0 ) / ( 1728.d0 * dxdz )
        e18(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
             * ( -3.d0 ) / ( 1728.d0 * dxdz )
        e19(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+1) &
             * ( +9.d0 ) / ( 1728.d0 * dxdz )
        if ( iz+2.le.nz+1) then
           e20(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+2) &
                * ( -1.d0 ) / ( 1728.d0 * dxdz )
        else
           e20(ix,iz) = 0.d0
        endif
        f1(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix-1,iz) + mu(ix,iz) ) &
             / ( 2.d0 * dx2 )
        f2(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix,iz) + mu(ix+1,iz) ) &
             / ( 2.d0 * dx2 )
        f3(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix,iz-1) + lam(ix,iz) ) &
             + 2.d0 * ( mu(ix,iz-1) + mu(ix,iz) ) ) &
             / ( 2.d0 * dz2 )
        f4(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix,iz) + lam(ix,iz+1) ) &
           + 2.d0 * ( mu(ix,iz) + mu(ix,iz+1) ) ) &
           / ( 2.d0 * dz2 )
        f5(ix,iz) = dt2 / rho(ix,iz) * mu(ix-1,iz) &
             / ( 4.d0 * dxdz )
        f6(ix,iz) = dt2 / rho(ix,iz) * mu(ix+1,iz) &
             / ( 4.d0 * dxdz )
        f7(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-1) &
             / ( 4.d0 * dxdz )
        f8(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz+1) &
             / ( 4.d0 * dxdz )
        if ( ix-2.ge.1 ) then
           f13(ix,iz) = dt2 / rho(ix,iz) * mu(ix-2,iz) &
                * (  1.d0 ) / ( 1728.d0 * dxdz )
        else
           f13(ix,iz) = 0.d0
        endif
        f14(ix,iz) = dt2 / rho(ix,iz) * mu(ix-1,iz) &
             * ( -9.d0 ) / ( 1728.d0 * dxdz )
        f15(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
             * (  3.d0 ) / ( 1728.d0 * dxdz )
        f16(ix,iz) = dt2 / rho(ix,iz) * mu(ix+1,iz) &
             * (  5.d0 ) / ( 1728.d0 * dxdz )
        if ( iz-2.ge.1 ) then
           f17(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-2) &
                * (  1.d0 ) / ( 1728.d0 * dxdz )
        else
           f17(ix,iz) = 0.d0
        endif
        f18(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-1) &
             * ( -9.d0 ) / ( 1728.d0 * dxdz )
        f19(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
             * (  3.d0 ) / ( 1728.d0 * dxdz )
        f20(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz+1) &
             * (  5.d0 ) / ( 1728.d0 * dxdz )
        
     enddo
  enddo
  
  return
end subroutine cales


subroutine calf( maxnz,it,t,ist,isx,isz,dt,dx,dz,rho,tp,ts,fx,fz )

  double precision pi
  parameter ( pi=3.1415926535897932d0 )
  integer maxnz,it,ist,isx,isz
  double precision t,dt,dx,dz,rho,tp,ts,fx(maxnz+1,*),fz(maxnz+1,*)
  double precision b
  
  if ( it.le.ist ) then
     b = pi * ( t - ts ) / tp
     fx(isx,isz) &
          = dsqrt(pi) / 2.d0 * (b*b-0.5d0) * dexp(-b*b) &
          / ( dx * dz )
     fx(isx,isz) = fx(isx,isz) * dt * dt / rho
     fz(isx,isz) = fx(isx,isz)
     if ( (it.eq.0).or.(it.eq.ist) ) then
        fx(isx,isz) = fx(isx,isz) / 2.d0
        fz(isx,isz) = fz(isx,isz) / 2.d0
     endif
  else
     fx(isx,isz) = 0.d0
     fz(isx,isz) = 0.d0
  endif

  ! NF for point source
  fx(isx,isx)=0.d0
  
  return
end subroutine calf

subroutine calf2( maxnz,it,t,ist,isx,isz,dt,dx,dz,rho,f0,t0,fx,fz )

  implicit none
  double precision pi
  parameter ( pi=3.1415926535897932d0 )
  integer maxnz,it,ist,isx,isz
  double precision t,dt,dx,dz,rho,f0,t0,fx(maxnz+1,*),fz(maxnz+1,*)
  double precision b,a,factor
  
  factor=1.d3
  
  if ( it.le.ist ) then

     a = pi*pi*f0*f0*(t-t0)*(t-t0)
     !print *,pi, t,t0,t-t0,a
     !print *, f0,t0,a
     ! Ricker source time function (second derivative of a Gaussian)
     fx(isx,isz) = factor * (1.d0 - 2.d0*a)*exp(-a);
     fx(isx,isz) = fx(isx,isz) * dt * dt / rho
     fz(isx,isz) = fx(isx,isz)
     if ( (it.eq.0).or.(it.eq.ist) ) then
        fx(isx,isz) = fx(isx,isz) / 2.d0
        fz(isx,isz) = fz(isx,isz) / 2.d0
     endif
  else
     fx(isx,isz) = 0.d0
     fz(isx,isz) = 0.d0
  endif

  !NF for point source
  fx(isx,isx)=0.d0
  
  return
end subroutine calf2


subroutine calStrainDiagonal(maxnz,nx,nz,ux,uz,lmargin,rmargin,singleStrainDiagonal)
  implicit none
  integer :: maxnz,nx,nz,ix,iz
  double precision :: ux(1:maxnz,1:maxnz),uz(1:maxnz,1:maxnz)
  integer :: lmargin(1:2), rmargin(1:2)
  real(kind(0e0)) :: singleStrainDiagonal(1:maxnz,1:maxnz)
  double precision :: straintmp
  double precision, parameter :: onetwelfth = 0.0833333333333333d0
  
  integer :: nxstart,nxend,nzstart,nzend
  ! we calculate only for the box of interest (without absorbing boundaries)
  ! and we suppose that we have lmargin != 0 and rmargin != 0
  ! i.e., we use the four-point first derivative operators with one extra point to the left

  nxstart=lmargin(1)+1
  nxend=nx+1-rmargin(1)

  nzstart=lmargin(2)+1
  nzend=nz+1-rmargin(2)

  
  do ix=nxstart,nxend
     do iz=nzstart,nzend
        straintmp=0.d0
        straintmp=(5.d0*ux(ix+1,iz)+3.d0*ux(ix,iz)-9.d0*ux(ix-1,iz)+ux(ix-2,iz))*onetwelfth
        straintmp=straintmp+(5.d0*uz(ix,iz+1)+3.d0*uz(ix,iz)-9.d0*uz(ix,iz-1)+uz(ix,iz-2))*onetwelfth
        singleStrainDiagonal(ix,iz)=straintmp
     enddo
  enddo
 
     
end subroutine calStrainDiagonal

subroutine calstep( maxnz,nx,nz, &
     e1, e2, e3, e4, e5, e6, e7, e8, &
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, &
     ux,uz,ux1,ux2,uz1,uz2,isx,isz,fx,fz, &
     work1,work2,work3,work4, &
     work5,work6,work7,work8, &
     work9,work10,work11,work12,optimise )

  integer maxnz,nx,nz,isx,isz
  double precision ux(maxnz+1,*),ux1(maxnz+1,*),ux2(maxnz+1,*)
  double precision uz(maxnz+1,*),uz1(maxnz+1,*),uz2(maxnz+1,*)
  double precision  e1(maxnz+1,*), e2(maxnz+1,*), e3(maxnz+1,*)
  double precision  e4(maxnz+1,*), e5(maxnz+1,*), e6(maxnz+1,*)
  double precision  e7(maxnz+1,*), e8(maxnz+1,*)
  double precision e13(maxnz+1,*),e14(maxnz+1,*),e15(maxnz+1,*)
  double precision e16(maxnz+1,*),e17(maxnz+1,*),e18(maxnz+1,*)
  double precision e19(maxnz+1,*),e20(maxnz+1,*)
  double precision  f1(maxnz+1,*), f2(maxnz+1,*), f3(maxnz+1,*)
  double precision  f4(maxnz+1,*), f5(maxnz+1,*), f6(maxnz+1,*)
  double precision  f7(maxnz+1,*), f8(maxnz+1,*)
  double precision f13(maxnz+1,*),f14(maxnz+1,*),f15(maxnz+1,*)
  double precision f16(maxnz+1,*),f17(maxnz+1,*),f18(maxnz+1,*)
  double precision f19(maxnz+1,*),f20(maxnz+1,*)
  double precision fx(maxnz+1,*),fz(maxnz+1,*)
  double precision work1(maxnz+1,-2:1),work2(maxnz+1,-1:2)
  double precision work3(maxnz+1,-2:1),work4(maxnz+1,-1:2)
  double precision work5(*),work6(maxnz+1,0:1)
  double precision work7(*),work8(maxnz+1,0:1)
  double precision work9(*),work10(maxnz+1,-2:1)
  double precision work11(*),work12(maxnz+1,-1:2)
  integer ix,iz,iz1,iz2,ix11,ix12,ix21,ix22
  logical optimise


  ! predicting the wavefield
 do iz=2,nz
    do ix=2,nx
       ux(ix,iz) = 2.d0 * ux1(ix,iz) - ux2(ix,iz) &
            + e1(ix,iz) * ( ux1(ix-1,iz) - ux1(ix,iz) ) &
            + e2(ix,iz) * ( ux1(ix+1,iz) - ux1(ix,iz) ) &
            + e3(ix,iz) * ( ux1(ix,iz-1) - ux1(ix,iz) ) &
            + e4(ix,iz) * ( ux1(ix,iz+1) - ux1(ix,iz) ) &
            - e5(ix,iz) * ( uz1(ix-1,iz+1) - uz1(ix-1,iz-1) ) &
            + e6(ix,iz) * ( uz1(ix+1,iz+1) - uz1(ix+1,iz-1) ) &
            - e7(ix,iz) * ( uz1(ix+1,iz-1) - uz1(ix-1,iz-1) ) &
            + e8(ix,iz) * ( uz1(ix+1,iz+1) - uz1(ix-1,iz+1) )
       uz(ix,iz) = 2.d0 * uz1(ix,iz) - uz2(ix,iz) &
            + f1(ix,iz) * ( uz1(ix-1,iz) - uz1(ix,iz) ) &
            + f2(ix,iz) * ( uz1(ix+1,iz) - uz1(ix,iz) ) &
            + f3(ix,iz) * ( uz1(ix,iz-1) - uz1(ix,iz) ) &
            + f4(ix,iz) * ( uz1(ix,iz+1) - uz1(ix,iz) ) &
            - f5(ix,iz) * ( ux1(ix-1,iz+1) - ux1(ix-1,iz-1) ) &
            + f6(ix,iz) * ( ux1(ix+1,iz+1) - ux1(ix+1,iz-1) ) &
            - f7(ix,iz) * ( ux1(ix+1,iz-1) - ux1(ix-1,iz-1) ) &
            + f8(ix,iz) * ( ux1(ix+1,iz+1) - ux1(ix-1,iz+1) )
    enddo
 enddo
 ux(isx,isz) = ux(isx,isz) + fx(isx,isz)
 uz(isx,isz) = uz(isx,isz) + fz(isx,isz)

 
 if(optimise) then
    ! correcting the wavefield
    !
    do ix=2,nx
       iz1 = 2
       iz2 = 3
       work1(ix,-2) = 0.d0
       work1(ix,-1) = 0.d0
       work1(ix,0) = 0.d0
       work1(ix,1) = ux(ix,iz1) - 2.d0 * ux1(ix,iz1) + ux2(ix,iz1)
       work2(ix,-1) = 0.d0
       work2(ix,0) = 0.d0
       work2(ix,1) = uz(ix,iz1) - 2.d0 * uz1(ix,iz1) + uz2(ix,iz1)
       work2(ix,2) = uz(ix,iz2) - 2.d0 * uz1(ix,iz2) + uz2(ix,iz2)
       work3(ix,-2) = 0.d0
       work3(ix,-1) = 0.d0
       work3(ix,0) = 0.d0
       work3(ix,1) = work1(ix,1) + 12.d0 * ux1(ix,iz1)
       work4(ix,-1) = 0.d0
       work4(ix,0) = 0.d0
       work4(ix,1) = work2(ix,1) + 12.d0 * uz1(ix,iz1)
       work4(ix,2) = work2(ix,2) + 12.d0 * uz1(ix,iz2)
    enddo
    
    do ix=1,nx+1
       ix11 = max0( ix-1,1 )
       ix12 = min0( ix+1,nx+1 )
       ix21 = max0( ix-2,1 )
       ix22 = min0( ix+2,nx+1 )
       work6(ix,0) = 0.d0
       work6(ix,1) = &
            (           ( -work3(ix11,1) ) &
            + 10.d0 * ( -work3(ix,  1) ) & 
            +         ( -work3(ix12,1) ) &
            )
       work8(ix,0) = 0.d0
       work8(ix,1) = &
            (           ( -work4(ix11,1) ) &
            + 10.d0 * ( -work4(  ix,1) ) &
            +         ( -work4(ix12,1) ) &
            )
       work10(ix,-2) = 0.d0
       work10(ix,-1) = 0.d0
       work10(ix,0) = 0.d0
       work10(ix,1)   = (          work3(ix21,1) - 9.d0 * work3(ix11,1) &
            + 3.d0 * work3(  ix,1) + 5.d0 * work3(ix12,1) )
       work12(ix,-1) = 0.d0
       work12(ix,0) = 0.d0
       work12(ix,1) = ( - 5.d0 * work4(ix11,1) - 3.d0 * work4(  ix,1) &
            + 9.d0 * work4(ix12,1) -        work4(ix22,1) )
       work12(ix,2) = ( - 5.d0 * work4(ix11,2) - 3.d0 * work4(  ix,2) &
            + 9.d0 * work4(ix12,2) -        work4(ix22,2) )
       
    enddo
    
    do iz=2,nz
       iz1 = iz + 1
       iz2 = min0( iz+2, nz+1 )
       do ix=2,nx
          work1(ix,-2) = work1(ix,-1)
          work1(ix,-1) = work1(ix,0)
          work1(ix,0) = work1(ix,1)
          work1(ix,1) = ux(ix,iz1) - 2.d0 * ux1(ix,iz1) + ux2(ix,iz1)
          work2(ix,-1) = work2(ix,0)
          work2(ix,0) = work2(ix,1)
          work2(ix,1) = work2(ix,2)
          work2(ix,2) = uz(ix,iz2) - 2.d0 * uz1(ix,iz2) + uz2(ix,iz2)
          work3(ix,-2) = work3(ix,-1)
          work3(ix,-1) = work3(ix,0)
          work3(ix,0) = work3(ix,1)
          work3(ix,1) = work1(ix,1) + 12.d0 * ux1(ix,iz1)
          work4(ix,-1) = work4(ix,0)
          work4(ix,0) = work4(ix,1)
          work4(ix,1) = work4(ix,2)
          work4(ix,2) = work2(ix,2) + 12.d0 * uz1(ix,iz2)
       enddo
       do ix=1,nx+1
          ix11 = max0( ix-1,1 )
          ix12 = min0( ix+1,nx+1 )
          ix21 = max0( ix-2,1 )
          ix22 = min0( ix+2,nx+1 )
          work5(ix) =   (           ( work3(ix11,-1)-work3(ix,-1) ) &
               + 10.d0 * ( work3(ix11, 0)-work3(ix, 0) ) &
               +         ( work3(ix11, 1)-work3(ix, 1) ) )
          work6(ix,0) = work6(ix,1)
          work6(ix,1) =  (  ( work3(ix11,0)-work3(ix11,1) ) &
               + 10.d0 * ( work3(  ix,0)-work3(ix,  1) ) &
               +         ( work3(ix12,0)-work3(ix12,1) ) )
          
          work7(ix) =( ( work4(ix11,-1)-work4(ix,-1) ) &
               + 10.d0 * ( work4(ix11, 0)-work4(ix, 0) ) &
               +         ( work4(ix11, 1)-work4(ix, 1) ))
          
          work8(ix,0) = work8(ix,1)
          work8(ix,1) =  (           ( work4(ix11,0)-work4(ix11,1) ) &
               + 10.d0 * ( work4(  ix,0)-work4(  ix,1) ) &
               +         ( work4(ix12,0)-work4(ix12,1) ))
          
          work9(ix) = (          work3(ix,-2) - 9.d0 * work3(ix,-1) &
               + 3.d0 * work3(ix,0)  + 5.d0 * work3(ix,1))
          
          work10(ix,-2) = work10(ix,-1)
          work10(ix,-1) = work10(ix,0)
          work10(ix,0) = work10(ix,1)
          work10(ix,1) = ( work3(ix21,1) - 9.d0 * work3(ix11,1) &
               + 3.d0 * work3(  ix,1) + 5.d0 * work3(ix12,1) )
          
          work11(ix) = ( - 5.d0 * work4(ix,-1)  - 3.d0 * work4(ix,0) &
               + 9.d0 * work4(ix, 1)  -        work4(ix,2) )
          
          work12(ix,-1) = work12(ix,0)
          work12(ix,0) = work12(ix,1)
          work12(ix,1) = work12(ix,2)
          work12(ix,2) = ( - 5.d0 * work4(ix11,2) - 3.d0 * work4(  ix,2) &
               + 9.d0 * work4(ix12,2) -        work4(ix22,2))
          
       enddo
       
       do ix=2,nx
          ix21 = max0( ix-2,1 )
          ix22 = min0( ix+2,nx+1 )
          ux(ix,iz) = ux(ix,iz) &
               + ( &
               - (           (   work1(ix-1,-1) + work1(ix-1,1) &
               + work1(ix+1,-1) + work1(ix+1,1) ) &
               + 10.d0 * (   work1(ix-1, 0) + work1(  ix,-1) &
               + work1(  ix, 1) + work1(ix+1, 0) ) &
               + 100.d0 * work1(ix,0) ) &
               + e1(ix,iz) * work5(ix) &
               - e2(ix,iz) * work5(ix+1) &
               + e3(ix,iz) * work6(ix,0) &
               - e4(ix,iz) * work6(ix,1) &
               ) / 144.d0 &
               + e13(ix,iz) * work11(ix-1) &
               + e14(ix,iz) * work11(ix) &
               + e15(ix,iz) * work11(ix+1) &
               + e16(ix,iz) * work11(ix22) &
               + e17(ix,iz) * work12(ix,-1) &
               + e18(ix,iz) * work12(ix,0) &
               + e19(ix,iz) * work12(ix,1) &
               + e20(ix,iz) * work12(ix,2)
          uz(ix,iz) = uz(ix,iz) &
               + ( &
               - (           (   work2(ix-1,-1) + work2(ix-1,1) &
               + work2(ix+1,-1) + work2(ix+1,1) ) &
               + 10.d0 * (   work2(ix-1, 0) + work2(  ix,-1) &
               + work2(  ix, 1) + work2(ix+1, 0) ) &
               + 100.d0 * work2(ix,0) ) &
               + f1(ix,iz) * work7(ix) &
               - f2(ix,iz) * work7(ix+1) &
               + f3(ix,iz) * work8(ix,0) &
               - f4(ix,iz) * work8(ix,1) &
               ) / 144.d0 &
               + f13(ix,iz) * work9(ix21) &
               + f14(ix,iz) * work9(ix-1) &
               + f15(ix,iz) * work9(ix) &
               + f16(ix,iz) * work9(ix+1) &
               + f17(ix,iz) * work10(ix,-2) &
               + f18(ix,iz) * work10(ix,-1) &
               + f19(ix,iz) * work10(ix,0) &
               + f20(ix,iz) * work10(ix,1)
       enddo
    enddo
     ux(isx,isz) = ux(isx,isz) + fx(isx,isz)
     uz(isx,isz) = uz(isx,isz) + fz(isx,isz)
 endif

 ! swapping u1 & u2 
 do iz=2,nz
    do ix=2,nx
       ux2(ix,iz) = ux1(ix,iz)
       ux1(ix,iz) =  ux(ix,iz)
       uz2(ix,iz) = uz1(ix,iz)
       uz1(ix,iz) =  uz(ix,iz)
    enddo
 enddo


 
 return
end subroutine calstep
	




subroutine create_color_image(image_data_2D,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
     NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,field_number)
	
  implicit none
  
  !       non linear display to enhance small amplitudes for graphics
  double precision, parameter :: POWER_DISPLAY = 1.d0
  
  !       amplitude threshold above which we draw the color point
  double precision, parameter :: cutvect = 0.01d0
  
  !       use black or white background for points that are below the threshold
  logical, parameter :: WHITE_BACKGROUND = .true.
  
  !       size of cross and square in pixels drawn to represent the source and the receivers
  integer, parameter :: width_cross=5,thickness_cross=1
  integer, parameter :: size_square=3
  
  integer NX,NY,it,field_number,ISOURCE,JSOURCE,NPOINTS_PML,nrec
  logical USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX
  
  double precision, dimension(NX,NY) :: image_data_2D
  
  integer, dimension(nrec) :: ix_rec,iy_rec
  logical, parameter :: kurama =.true.
  integer :: ix,iy,irec
  
  character(len=100) :: file_name,system_command1,system_command2,system_command3
	
  integer :: R, G, B
  double precision :: R1, G1, B1
  
  double precision :: normalized_value,max_amplitude
  
  !       open image file and create system command to convert image to more convenient format
  !       use the "convert" command from ImageMagick http://www.imagemagick.org
  if(field_number == 1) then
     write(file_name,"('image',i6.6,'_Ux.pnm')") it
     write(system_command1, "('convert image',i6.6,'_Ux.pnm snapshots/imageUx',i6.6,'.png')") it,it
     write(system_command2, "('rm image',i6.6,'_Ux.pnm')") it
  else if(field_number == 2) then
     write(file_name,"('image',i6.6,'_Uz.pnm')") it
     write(system_command1,"('convert image',i6.6,'_Uz.pnm snapshots/imageUz',i6.6,'.png')") it,it
     write(system_command2,"('rm image',i6.6,'_Uz.pnm')") it
  endif
  
  open(unit=27, file=file_name, status='unknown')
  
  write(27,"('P3')")	! write image in PNM P3 format
  
  write(27,*) NX,NY	! write image size
  write(27,*) '255'	! maximum value of each pixel color
	
  !       compute maximum amplitude
  max_amplitude = maxval(abs(image_data_2D))
  
  !       image starts in upper-left corner in PNM format
  do iy=1,NY,1
     do ix=1,NX
        
        !       define data as vector component normalized to [-1:1] and rounded to nearest integer
        !       keeping in mind that amplitude can be negative
        normalized_value = image_data_2D(ix,iy) / max_amplitude *1.5d0
        
        !       suppress values that are outside [-1:+1] to avoid small edge effects
        if(normalized_value < -1.d0) normalized_value = -1.d0
        if(normalized_value > 1.d0) normalized_value = 1.d0
        
        !       draw an orange cross to represent the source
        if((ix >= ISOURCE - width_cross .and. ix <= ISOURCE + width_cross .and. &
             iy >= JSOURCE - thickness_cross .and. iy <= JSOURCE + thickness_cross) .or. &
             (ix >= ISOURCE - thickness_cross .and. ix <= ISOURCE + thickness_cross .and. &
             iy >= JSOURCE - width_cross .and. iy <= JSOURCE + width_cross)) then
           R = 255
           G = 157
           B = 0
           
           !       display two-pixel-thick black frame around the image
        else if(ix <= 2 .or. ix >= NX-1 .or. iy <= 2 .or. iy >= NY-1) then
           R = 0
           G = 0
           B = 0
           
           !       display edges of the PML layers
        else if((USE_PML_XMIN .and. ix == NPOINTS_PML) .or. &
             (USE_PML_XMAX .and. ix == NX - NPOINTS_PML) .or. &
             (USE_PML_YMIN .and. iy == NPOINTS_PML) .or. &
             (USE_PML_YMAX .and. iy == NY - NPOINTS_PML)) then
           R = 255
           G = 150
           B = 0
           
           !       suppress all the values that are below the threshold
        else if(abs(image_data_2D(ix,iy)) <= max_amplitude * cutvect) then
           
           !       use a black or white background for points that are below the threshold
           if(WHITE_BACKGROUND) then
              R = 255
              G = 255
              B = 255
           else
              R = 0
              G = 0
              B = 0
           endif
           
           !       represent regular image points using red if value is positive, blue if negative
        elseif(kurama) then
           normalized_value=normalized_value**POWER_DISPLAY
           call plotcolor(normalized_value,R1,G1,B1)
           R = nint(R1)
           G = nint(G1)
           B = nint(B1)
        else
           if(normalized_value >= 0.d0) then
              R = 255
              G = nint(255.d0-255.d0*normalized_value**POWER_DISPLAY)
              B = nint(255.d0-255.d0*normalized_value**POWER_DISPLAY)
           else
              R = nint(255.d0-255.d0*abs(normalized_value)**POWER_DISPLAY)
              G = nint(255.d0-255.d0*abs(normalized_value)**POWER_DISPLAY)
              B = 255
           endif
        endif
        
        !       draw a green square to represent the receivers
        do irec = 1,nrec
           if((ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
                iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square) .or. &
                (ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
                iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square)) then
              !       use dark green color
              R = 30
              G = 180
              B = 60
           endif
        enddo
        
        !       write color pixel
        write(27,"(i3,' ',i3,' ',i3)") R,G,B
        
     enddo
  enddo
  
  !       close file
close(27)

!       call the system to convert image to GIF (can be commented out if "call system" is missing in your compiler)
call system(system_command1)
call system(system_command2)
end subroutine create_color_image



subroutine  compNRBC2(ux,ux1,ux2,uz,uz1,uz2, rrate, lmargin, rmargin,nnx,nnz)

  ! Cerjan boundary conditions (2D)
  implicit none
  integer :: nnx, nnz

  real*8, intent(inout) :: ux2(nnx,nnz),uz2(nnx,nnz)
  real*8, intent(inout) :: ux1(nnx,nnz),uz1(nnx,nnz)
  real*8, intent(inout) :: ux(nnx,nnz),uz(nnx,nnz)
  real*8, intent(in) :: rrate
  integer, dimension(3), intent(in) :: lmargin, rmargin
  integer ix, iy, iz
  integer i, j, k
  real*8 r
  
  do iz = 1, nnz
     !do iy = 1, nny
     do ix = 1, nnx
        
        i = 0
        j = 0
        k = 0
           
        if (ix < lmargin(1) + 1) i = lmargin(1) + 1 - ix
        !   if (iy < lmargin(2) + 1) j = lmargin(2) + 1 - iy
        if (iz < lmargin(2) + 1) k = lmargin(2) + 1 - iz
        if (nnx - rmargin(1) < ix) i = ix - nnx + rmargin(1)
        !if (nny - rmargin(2) < iy) j = iy - nny + rmargin(2)
        if (nnz - rmargin(2) < iz) k = iz - nnz + rmargin(2)
           
        if (i == 0 .and. j == 0 .and. k == 0) cycle
        
        r = rrate * rrate * dble( i * i + j * j + k * k )
        r = exp( - r )
        

        ux2(:,:) = ux2(:,:) * r
        ux1(:,:) = ux1(:,:) * r
        ux(:,:) = ux(:,:) * r

        uz2(:,:) = uz2(:,:) * r
        uz1(:,:) = uz1(:,:) * r
        uz(:,:) = uz(:,:) * r
        
     enddo
     !enddo
  enddo
  
end subroutine compNRBC2



subroutine  compNRBCpre(r,rrate, lmargin, rmargin,nnx,nnz)

  implicit none
  integer :: nnx, nnz
  ! Cerjan boundary conditions (2D)
  double precision :: r(nnx,nnz)
  real*8, intent(in) :: rrate
  integer, dimension(3), intent(in) :: lmargin, rmargin
  integer ix, iy, iz
  integer i, j, k
  double precision :: rr
  
  
  do iz = 1, nnz
     !do iy = 1, nny
     do ix = 1, nnx
        
        i = 0
        j = 0
        k = 0
           
        if (ix < lmargin(1) + 1) i = lmargin(1) + 1 - ix
        !   if (iy < lmargin(2) + 1) j = lmargin(2) + 1 - iy
        if (iz < lmargin(2) + 1) k = lmargin(2) + 1 - iz
   
        if (nnx - rmargin(1) < ix) i = ix - nnx + rmargin(1)
        !if (nny - rmargin(2) < iy) j = iy - nny + rmargin(2)
        if (nnz - rmargin(2) < iz) k = iz - nnz + rmargin(2)
           
        if (i == 0 .and. j == 0 .and. k == 0) cycle
        
        rr = rrate * rrate * dble( i * i + j * j + k * k )
        r(ix,iz) = exp( - rr )
        
        !if(r(ix,iz).ne.1.d0) then
        !   print *, ix,iz,r(ix,iz)
        !endif
        
        !print *, ix,iy,r
        !ux2(:,:) = ux2(:,:) * r
        !ux1(:,:) = ux1(:,:) * r
        !ux(:,:) = ux(:,:) * r

        !uz2(:,:) = uz2(:,:) * r
        !uz1(:,:) = uz1(:,:) * r
        !uz(:,:) = uz(:,:) * r
        
     enddo

     
     
     !enddo
  enddo
  
end subroutine compNRBCpre






subroutine plotcolor(v,r1,g1,b1)
!======================================================================
! Interpolate r1 g1 b1
!======================================================================

  real(8) v,v1,nl
  double precision r1, g1, b1
  real(8),allocatable,dimension(:) :: x,r,g,b

  integer n,i

  !read colormap data
  open (17, file='./colormap/colormap.dat', status='old')
  read (17, *) nl

  n = int(nl)
  allocate( x(n) )
  allocate( r(n) )
  allocate( g(n) )
  allocate( b(n) )

  do i = 1, int(nl)
      read (17, *) x(i), r(i), g(i), b(i)
  end do
  close(17)
    !rescale normalized value from [-1 1] to [0 1]
    !Depending on colormap.dat

    !when colormap domain is [0 1]
    !v1 = (1.0d0+v)/2.0d0

    !when colormap domain is [-1 1]
    v1 = v

    call colormap(v1,n,x,r,g,b,r1,g1,b1)

end subroutine plotcolor

subroutine colormap(v,n,x,r,g,b,r1,g1,b1)
!======================================================================
  !v: normalized value [-1 1]
  !n: number of the data in colormap
  !x: ampritude in dataset
  !r: dataset of R
  !g: dataset of G
  !b: dataset of B
  !r1: interpolated R
  !g1 interpolated G
  !b1: interpolated B
!----------------------------------------------------------------------

  implicit none
  integer i,n
  real(8) v,x(n),r(n),g(n),b(n)
  double precision ispline, r1, g1, b1
  real(8) b_spline(n),c_spline(n),d_spline(n)

  !for R
  call spline (x, r, b_spline, c_spline, d_spline,n)
  r1 = ispline(v, x, r, b_spline, c_spline, d_spline, n)

  !for G
  call spline (x, g, b_spline, c_spline, d_spline,n)
  g1 = ispline(v, x, g, b_spline, c_spline, d_spline, n)

  !for B
  call spline (x, b, b_spline, c_spline, d_spline,n)
  b1 = ispline(v, x, b, b_spline, c_spline, d_spline, n)


end subroutine colormap

!======================================================================
!======================================================================
!======================================================================


subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer n
double precision x(n), y(n), b(n), c(n), d(n)
integer i, j, gap
double precision h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline

function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
double precision ispline
integer n
double precision  u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline
