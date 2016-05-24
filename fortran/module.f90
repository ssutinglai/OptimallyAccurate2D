module parameters

  implicit none

  ! Nobuaki Fuji, 2016, Institut de Physique du Globe de Paris
  
  ! some parameters to handle

  real(kind(0d0)) :: fmax ! in Hz
  real(kind(0d0)) :: vmin, vmax ! max and min of wave speeds in m/s
  integer :: NX, NZ ! number of points in x and z direction
  integer :: NSTEP


  real(kind(0d0)) :: disp_z, disp_x ! dispersion condition (5 points per wavelength)
  real(kind(0d0)) :: condstability ! stability condition

  real(kind(0d0)) :: dt,dx,dz ! in m and s
  

  real(kind(0d0)), allocatable :: rho(:,:),mu(:,:),lambda(:,:),cs(:,:),cp(:,:)

  real(kind(0d0)), allocatable :: nx_vec(:), nz_vec(:), nt_vec(:), source(:) ! for visualisation

  
  

  
  

  
