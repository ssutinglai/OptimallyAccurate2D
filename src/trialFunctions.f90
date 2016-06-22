program SincInterpolation( nx,nz,rho,lam,mu,dx,dz,dt )
  implicit none
  integer nx,nz
  integer ix,iz
  double precision rho(nx+1,nz+1),lam(nx+1,nz+1),mu(nx+1,nz+1)
  double precision dx,dz,dt
  double precision dx2,dz2,dxdz,dt2
  double precision phi(ix,iz)
  parameter (pi = 3.141592653589793238462643383)

  dt2 = dt * dt
  dx2 = dx * dx
  dz2 = dz * dz
  dxdz = dx * dz

subroutine trialfunc
do ix=2,nx
   do iz=2,nz
	phi(ix,iz)=sin(pi*ix)*sin(pi*iz)/(pi*pi*ix*iz)
   enddo
enddo

   
