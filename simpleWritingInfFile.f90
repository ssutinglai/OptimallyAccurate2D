program writingInf
  implicit none
  character(120) :: filename
  character(80), parameter :: modelname = 'layered'
  character(80), parameter :: vpmodel = './2d_sdafsa.vp'
  !vs
  !rho
  integer :: nt, nx, nz
  double precision :: dt, dx, dz
  integer :: isx, isz ! source position
  ! for Ricker wavelets
  double precision :: f0, t0
  integer, parameter :: nReceiver = 199
  integer :: nrx(1:nReceiver), nrz(1:nReceiver)
  
  integer, parameter :: nSource = 199
  
  integer :: iSource, j, iReceiver

  do iReceiver = 1, nReceiver
     nrx(iReceiver)=2*iReceiver
     nrz(iReceiver)=1
  enddo
  
  do iSource = 1, nSource
     isx=2*iSource
     isz=1
     write(filename, '(I5,".".I5,".inf")') isx,isz
     do j=1, 12
        if(filename(j:j).eq.' ') filename(j:j:)='0'
     enddo
     filename= './'//trim(modelname)//'.'//filename
     
     

     open(1, file=filename, form='formatted')
     write(1,*) 'c grids'
     write(1,*) nt, nx, nz
     write(1,*) dt, dx, dz
     write(1,*) 'c modelname'
     write(1,*) modelname 

     ! etc
     
     write(1,*) 'c receivers information'
     write(1,*) nReceiver
     do iReceiver = 1, nReceiver
        write(1,*) nrx(iReceiver), nrz(iReceiver)
     enddo
     

     write(1,*) 'c'
     write(1,*) 'end'
  enddo
  
  
  
  

  
