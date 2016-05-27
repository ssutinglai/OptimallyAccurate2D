program multipleSourcesOPT2D

  ! Computation of the synthetic seismograms in the time domain
  ! using the optimally accurate operators.
  ! 2D PSV heterogeneous medium
  ! CPML or Cerjan boundary conditions
  !
  !					originally from 1997.6  N.Takeuchi
  !                                                     2016.5. N.Fuji
  

  implicit none
  character(120) :: filename
  character(80), parameter :: modelname = 'elips'
  character(80), parameter :: vpmodel = './2d_elipsoid.vp'
  character(80), parameter :: vsmodel = './2d_elipsoid.vs'
  !character(80), parameter :: modelname = 'layered'
  !character(80), parameter :: vpmodel = './2d_start.vp'
  !character(80), parameter :: vsmodel = './2d_start.vs'

  character(80), parameter :: rhomodel = './2d_start.rho'
  !integer :: nt, nx, nz
  !real :: dt, dx, dz
  ! integer :: isx, isz ! source position
  ! for Ricker wavelets
  !real :: f0, t0
  integer, parameter :: nReceiver = 199
  integer :: nrx(1:nReceiver), nrz(1:nReceiver)
  
  integer, parameter :: nSource = 199
  integer :: iisx(1:nSource), iisz(1:nSource)
  integer :: iSource, iReceiver



  integer, parameter :: maxnz = 600 
  integer, parameter :: maxnt = 2000
  double precision, parameter :: pi=3.1415926535897932d0 
  double precision, parameter :: ZERO = 0.d0
    
  !
  ! parameters for the gridding
  double precision dt,dx,dz
  ! parameters for the wavefield
  integer nt,nx,nz,it,ist,isx,isz,ix,iz,recl_size
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
  ! parameter for the structure
  !double precision rrho(6),llam(6),mmu(6)
  character(80) :: vpfile, vsfile, rhofile   ! ,modelname
  double precision :: rho(maxnz+1,maxnz+1)
  double precision :: lam(maxnz+1,maxnz+1),mu(maxnz+1,maxnz+1)
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
  character(120) :: outfile
 
  
  ! parameter for the waveform
  double precision t
  !parameter for video
  real,dimension(maxnz+1,maxnz+1) :: snapux,snapuz
  integer, parameter :: IT_DISPLAY = 10
 
  integer(2) head(1:120)
  character(80) :: routine
  
  logical,parameter :: dummylog = .false.
  ! switch OPT / CONV
  logical,parameter :: optimise = .true.
  ! switch C-PML 
  logical, parameter :: USE_PML_XMIN = .true.
  logical, parameter :: USE_PML_XMAX = .true.
  logical, parameter :: USE_PML_YMIN = .true.
  logical, parameter :: USE_PML_YMAX = .true.
  ! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 40
  double precision, parameter :: CerjanRate = 0.015
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

  call system('mkdir ./inffile')


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


  do iReceiver = 1, nReceiver
     nrx(iReceiver)=2*iReceiver
     nrz(iReceiver)=1
  enddo
  
  do iSource = 1, nSource
     iisx(iSource)=2*iSource
     iisz(iSource)=1
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



  

  


  
  
  
  
  

end program multipleSourcesOPT2D
