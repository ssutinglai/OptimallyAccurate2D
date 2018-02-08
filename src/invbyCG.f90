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
  double precision :: VAR(0:(boxnx+1)*(boxnz+1)*2)
  logical :: doCG

  print *, "inversion by CG"
  open(unit=1,form="unformatted",file='ataatd')
  read(1) ata,atd
  close(1)
  print *, "ata read"

  !@ata=real(ata)
  !atd=real(atd)

  print *, ata(1:81,200:200)
  
  !stop

  doCG=.true.
  
  AIC=0.d0
  VAR=0.d0
  

  x0 = 0.d0

  r = atd ! - matmul(ata,x0)
  w = -r
  !z = matmul(ata,w)
  call ataMatmul(w,z)
  a = dot_product(conjg(r),w) / dot_product(conjg(w),z)
  x = x0 +a*w
  b = 0
  

  ii=0
  ND = dble(nnFreq*nSource*nReceiver)/alphaAIC
  VAR(ii) = dot_product(conjg(r),r)
  AIC(ii) = ND*log(2.d0*pi)+ND*log(VAR(ii))+ND+2.d0*dble(ii+1)
  
  

  do while (doCG)
 
     ii=ii+1

     r = r - a*z
     b = dot_product(conjg(r),z)/dot_product (conjg(w),z)
     w = -r + b*w
     !z = matmul(ata,w)
     call ataMatmul(w,z)
     a = dot_product(conjg(r),w)/dot_product(conjg(w),z)
     x = x+a*w

     VAR(ii) = dot_product(conjg(r),r)
     AIC(ii) = ND*log(2.d0*pi)+ND*log(VAR(ii))+ND+2.d0*dble(ii+1)
     
     write(13,*) ii, VAR(ii), AIC(ii)
     
     if(AIC(ii)>AIC(ii-1)) then
        doCG=.false.
     else
        x0 = x ! x0 to be updated
     endif
  enddo
  
  print *, ii, " CG vectors were used"
  open(1,form='formatted',file='aicvar.dat')
  do jj=0,ii
     write(1,*) jj,VAR(jj),AIC(jj)
  enddo
  close(1)

  
  do ixz=1,(boxnx+1)*(boxnz+1)
     iz=(ixz-1)/(boxnx+1)+1
     ix=mod(ixz-1,boxnx+1)+1
   
     kernelP(ix,iz)=dble(x0(2*(ixz-1)+1))
     kernelS(ix,iz)=dble(x0(2*(ixz-1)+2))
  enddo

  

  return
  
end subroutine invbyCG
  


subroutine ataMatmul(w,z)
  ! ata is conjugate transpose
  ! w is input; z is output
  use paramFWI
  use parameters
  implicit none
  integer :: ixz
  integer :: jxz,jx,jz
  integer :: jxzlocal
  integer :: iTypeParam,jTypeParam ! 1 for Vp and 2 for Vs
  ! matrix multiplication with a truncated ata
  double complex :: w(1:(boxnx+1)*(boxnz+1)*2)
  double complex :: z(1:(boxnx+1)*(boxnz+1)*2)

  z=cmplx(0.d0)

  do ixz=1,(boxnx+1)*(boxnz+1)
     iz=(ixz-1)/(boxnx+1)+1
     ix=mod(ixz-1,boxnx+1)+1     
     do iTypeParam=1,2
        do jz=max(iz-nNeighbours/2,1),min(iz+nNeighbours/2,boxnz+1)
           do jx=max(ix-nNeighbours/2,1),min(ix+nNeighbours/2,boxnx+1)
              jxz=(jz-1)*(boxnx+1)+jx
              jxzlocal=(jz-iz+nNeighbours/2)*nNeighbours+(jx-ix+nNeighbours/2+1)
              do jTypeParam=1,2
                 z(2*(ixz-1)+iTypeParam)=&
                      z(2*(ixz-1)+iTypeParam)+&
                      ata(2*(jxzlocal-1)+jTypeParam,2*(ixz-1)+iTypeParam)*&
                      conjg(w(2*(jxz-1)+jTypeParam))
              enddo
           enddo
        enddo
     enddo
  enddo
                 
              
  return
end subroutine ataMatmul
