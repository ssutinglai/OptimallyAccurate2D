program SincInterpolation !( nx,nz,rho,lam,mu,dx,dz,dt )
  implicit none
  integer nx,nz
  integer ix,iz
  integer jx,jz
  integer m,n
  double precision xm,zm
  integer ndis,ngrid,npTF
  logical, parameter :: sincfunction = .FALSE.
  double precision x,xx,z,zz
  !double precision rho(nx+1,nz+1),lam(nx+1,nz+1),mu(nx+1,nz+1)
  double precision, allocatable :: lam(:,:)
  double precision dx,dz,dt
  double precision dx2,dz2,dxdz,dt2
  double precision, allocatable :: phix(:,:),phiz(:,:)
  double precision, allocatable :: phixderiv(:,:),phizderiv(:,:)
  double precision, allocatable :: H1(:,:),H2(:,:)
  double precision, parameter :: pi = 3.141592653589793238462643383
 

  dt = 1.d0
  dx = 1.d0
  dz = 1.d0


  

  
  !npTF defines points in scheme (3, 5, 7)
  npTF = 3
  ngrid = (npTF-1)/2
  ndis = 100


  allocate(phix(-ngrid*ndis:ngrid*ndis,-ngrid*ndis:ngrid*ndis))
  allocate(phiz(-ngrid*ndis:ngrid*ndis,-ngrid*ndis:ngrid*ndis))
  allocate(phixderiv(-ngrid*ndis:ngrid*ndis,-ngrid*ndis:ngrid*ndis))
  allocate(phizderiv(-ngrid*ndis:ngrid*ndis,-ngrid*ndis:ngrid*ndis))
  allocate(lam(-ngrid*ndis:ngrid*ndis,-ngrid*ndis:ngrid*ndis))
  allocate(H1(1,1))
  allocate(H2(1,1))
  

  lam =1.d0

  dt2 = dt * dt
  dx2 = dx * dx
  dz2 = dz * dz
  dxdz = dx * dz
  
  xm = 0.d0
  zm = 0.d0

  m = 1
  n = 1
  
  !trialfunction decides on sinc(true) or linear(false) interpolation
  !sincfunction = .true.


  
  if (sincfunction) then
     
     do ix=-ngrid*ndis, ngrid*ndis
        do iz=-ngrid*ndis, ngrid*ndis
           x=xm+dble(ix/ndis)*dx
           z=zm+dble(iz/ndis)*dz
           xx=(x-xm)/dx
           zz=(z-zm)/dz
           !phi(ix,iz)=sin(pi*xx)*sin(pi*zz)/(pi*pi*xx*zz)
           phix(ix,iz)=sin(pi*xx)/pi*xx
           phiz(ix,iz)=sin(pi*zz)/pi*zz
           phixderiv(ix,iz)=pi*cos(pi*xx)/(pi*xx) - sin(pi*xx)/(pi*xx*xx)
           phizderiv(ix,iz)=pi*cos(pi*zz)/(pi*zz) - sin(pi*zz)/(pi*zz*zz)

           open(unit=8,file="phix.dat",form="formatted"&
           ,status="replace",action="write")
           
           write(8,*)'phix', phix(ix,iz)
           
           close(8)
           
        enddo
     enddo
     
     
     
  else
     
     do ix=-ngrid*ndis,0
        do iz=-ngrid*ndis,0
           x=xm+dble(ix/ndis)*dx
           z=zm+dble(iz/ndis)*dz
           !phi(ix,iz)=(x+dx)*(z+dz)/dx*dz
           phix(ix,iz)=(x+dx)/dx
           phiz(ix,iz)=(z+dz)/dz
           phixderiv=1/dx
           phizderiv=1/dz
        end do
     end do
     
     do ix=0,ngrid*ndis
        do iz=0,ngrid*ndis
           x=xm+dble(ix/ndis)*dx
           z=zm+(iz/ndis)*dz
           !phi(ix,iz)=(-x+dx)*(-z+dz)/dx*dz
           phix(ix,iz)=(-x+dx)/dx
           phiz(ix,iz)=(-z+dz)/dz
           phixderiv(ix,iz)=-1/dx
           phizderiv(ix,iz)=-1/dz
        end do
     end do
  endif

     do jx=-ngrid, ngrid
        do jz=-ngrid, ngrid
           
           H1(m,n)=(phix(jx,jz)*phix(jx,jz))+ &
                (phiz(jx,jz)*phiz(jx,jz))
           
           
           H2(m,n)=lam(jx,jz)*(phixderiv(jx,jz)+phizderiv(jx,jz))&
               *(phixderiv(jx,jz)+phizderiv(jx,jz))

           print *,'jx,jz',jx,jz,'H1', H1(m,n), 'H2', H2(m,n)
           
        end do
     end do
     






end program SincInterpolation
  
  
