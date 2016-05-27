program writingInf
  implicit none
  character(120) :: filename
  character(80), parameter :: modelname = 'elips'
  character(80), parameter :: vpmodel = './2d_elipsoid.vp'
  character(80), parameter :: vsmodel = './2d_elipsoid.vs'
  !character(80), parameter :: modelname = 'layered'
  !character(80), parameter :: vpmodel = './2d_start.vp'
  !character(80), parameter :: vsmodel = './2d_start.vs'

  character(80), parameter :: rhomodel = './2d_start.rho'
  integer :: nt, nx, nz
  real :: dt, dx, dz
  integer :: isx, isz ! source position
  ! for Ricker wavelets
  real :: f0, t0
  integer, parameter :: nReceiver = 199
  integer :: nrx(1:nReceiver), nrz(1:nReceiver)
  
  integer, parameter :: nSource = 199
  
  integer :: iSource, j, iReceiver

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
     isx=2*iSource
     isz=1
     write(filename, '(I5,".",I5,".inf")') isx,isz
     do j=1, 12
        if(filename(j:j).eq.' ') filename(j:j)='0'
     enddo
    filename= './'//trim(modelname)//'.'//filename
    
       

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
     write(1,*) isx,isz
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
  
  
  
  

 end program writingInf 
