subroutine invbyCG
  use paramFWI
  use parameters
  implicit none

  double complex :: a, b, tmp_r2
  integer :: ii,jj
  double complex :: x(1:(nx+1)*(nz+1)*2)
  double complex :: r(1:(nx+1)*(nz+1)*2)
  double complex :: w(1:(nx+1)*(nz+1)*2)
  double complex :: z(1:(nx+1)*(nz+1)*2)
  double complex :: x0(1:(nx+1)*(nz+1)*2)
 
  character(3) :: num
  double precision :: ND
  double precision :: AIC(0:(nx+1)*(nz+1)*2)
  logical :: doCG

  doCG=.true.
  
  AIC=0.d0
  
  x0 = 0.d0

  r = atd - matmul(ata,x0)
  w = -r
  z = matmul(ata,w)
  a = dot_product(r,w) / dot_product(w,z)
  x = x0 +a*w
  b = 0
  

  ii=0
  ND = dble(nFreq*nSource*nReceiver)/alphaAIC
  AIC(ii) = ND*log(2.d0*pi)+ND*log(dot_product(conjg(r),r))+ND+2.d0*dble(ii+1)
  
  

  do while (doCG)
     ii=ii+1

     

     r = r - a*z
     b = dot_product(r,z)/dot_product (w,z)
     w = -r + b*w
     z = matmul(ata,w)
     a = dot_product(r,w)/dot_product(w,z)
     x = x+a*w

     AIC(ii) = ND*log(2.d0*pi)+ND*log(dot_product(conjg(r),r))+ND+2.d0*dble(ii+1)
     if(AIC(ii)>AIC(ii-1)) then
        doCG=.false.
     else
        x0 = x ! x0 to be updated
     endif
  enddo


  atd=x0
  

  return
  
end subroutine invbyCG
  
