subroutine initialise1

  ! Nobuaki Fuji, 2016, Institut de Physique du Globe de Paris

  use paramters

  implicit none

  integer i

  fmax = 40.d0
  
  vmin = 1500.d0
  vmax = 2500.d0

  NX = 210
  NZ = 420
  NSTEP = 365

  disp_z=vmin/fmax/5.d0
  disp_x=vmin/fmax/5.d0
  
  condstability = disp_x/vmax/sqrt(2.d0)* 8.d-1

  dt = 0.0013d0/2.d0
  dx = 6.d0/2.d0
  dz = 6.d0/2.d0

  if ((dx>disp_x).or.(dz>disp_z)) then
     print *, "violation of the dispersion condition"
     stop
  endif


  if (dt>condstability) then
     print *, "violation of the stability condition"
     stop
  endif


  allocate(rho(0:NZ,0:NX))
  allocate(mu(0:NZ,0:NX))
  allocate(lambda(0:NZ,0:NX))
  allocate(cs(0:NZ,0:NX))
  allocate(cp(0:NZ,0:NX))
  
  rho(:,:) = 2400.d0

  cp(0:NZ/2,:) = 1800.d0
  cp(NZ/2+1:NZ,:) = 3000.d0

  cs(0:NZ/2,:) = 1800.d0/1.732d0
  cs(NZ/2+1:NZ,:) = 3000.d0/1.732d0

  
  allocate(nx_vec(0:NX))
  allocate(nz_vec(0:NZ))
  allocate(nt_vec(0:NSTEP))

  do i = 0, NX
     nx_vec(i) = dx*dble(i)
  enddo
  
  do i = 0, NZ
     nz_vec(i) = dz*dble(i)
  enddo

  do i = 0, NSTEP
     nt_vec(i) = dt*dble(i)
  enddo


  

  
  
  
 
