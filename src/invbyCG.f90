subroutine invbyCG
  use paramFWI
  use parameters
  implicit none

  double complex :: a, b, tmp_r2
  integer :: ii,jj,ixz

  double complex :: x(1:(boxnx+1)*(boxnz+1)*2)
  double complex :: r(1:(boxnx+1)*(boxnz+1)*2)
  double complex :: w(1:(boxnx+1)*(boxnz+1)*2)
  double complex :: z(1:(boxnx+1)*(boxnz+1)*2)
  double complex :: x0(1:(boxnx+1)*(boxnz+1)*2)
 
  character(3) :: num
  double precision :: ND
  double precision :: AIC(0:(boxnx+1)*(boxnz+1)*2)
  logical :: doCG

  print *, "inversion by CG"
  open(unit=1,form="unformatted",file='ataatd')
  read(1) ata,atd
  close(1)
  print *, "ata read"

  


  doCG=.true.
  
  AIC=cmplx(0.d0)
  
  x0 = 0.d0

  r = atd ! - matmul(ata,x0)
  w = -r
  z = matmul(ata,w)
  a = dot_product(r,w) / dot_product(w,z)
  x = x0 +a*w
  b = 0
  

  ii=0
  ND = dble(nnFreq*nSource*nReceiver)/alphaAIC
  AIC(ii) = ND*log(2.d0*pi)+ND*log(dot_product(conjg(r),r))+ND+2.d0*dble(ii+1)
  
  

  do while (doCG)

     print *, ii
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
  
  print *, ii, " CG vectors were used"

  do ixz=1,(boxnx+1)*(boxnz+1)
     iz=(ixz-1)/(boxnx+1)+1
     ix=mod(ixz-1,boxnx+1)+1
     kernelP(ix,iz)=dble(x0(2*(ixz-1)+1))
     kernelS(ix,iz)=dble(x0(2*(ixz-1)+2))
  enddo



  return
  
end subroutine invbyCG
  
