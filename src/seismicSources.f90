
subroutine calf2(nx,nz,it,t,ist,isx,isz,dt,dx,dz,rho,f0,t0,fx,fz )

  implicit none
  real*8 pi
  parameter ( pi=3.1415926535897932d0 )
  integer nx,nz,it,ist,isx,isz,i
  real*8 t,dt,dx,dz,rho,f0,t0,fx(nx+1,nz+1),fz(nx+1,nz+1)
  real*8 b,a,factor
  real*8 wt(100000),wx(100000)

  !%! Read wavelet file from Katayama

  open(1,file='Wavelet.txt')
  do i=1,100000
  read (1,*) wt(i), wx(i)
  enddo
  close(1)

   factor=1.d3
  !factor=5.d17
  !factor=1.d3
  !factor =1.d-7
  if ( it.le.10000 ) then
     !a = pi*pi*f0*f0*(t-t0)*(t-t0) (04.04.2018)
     !print *,pi, t,t0,t-t0,a
     !print *, f0,t0,a

     ! Ricker source time function (second derivative of a Gaussian)
     !%! Change for fitting with Specfem2D (08.03.2018)
     !fx(isx,isz) = -factor * (1.d0 - 2.d0*a)*exp(-a); !(04.04.2018)
     if(it.le.1)then
     fx(isx,isz) = 0.d0
     fz(isx,isz) = 0.d0
     else
     fx(isx,isz) = wx(it)
     !fx(isx,isz) = fx(isx,isz) * dt * dt / rho !%! Closed for matching with Specfem2D (08.03.2018)
     fz(isx,isz) = fx(isx,isz)
     endif
     !%! Commented for matching for Specfem2D (08.03.2018)
     !if ( (it.eq.0).or.(it.eq.ist) ) then
     !   fx(isx,isz) = fx(isx,isz) / 2.d0
     !   fz(isx,isz) = fz(isx,isz) / 2.d0
     !endif
!write(*,*)t+t0
!open(13,file='waveletOpt',form='formatted',access='append',status='old')
!write(13,*) it, fz(isx,isz)
!close(13)

  else
     fx(isx,isz) = 0.d0
     fz(isx,isz) = 0.d0

! Output Wavelet

  endif
!open(13,file='waveletOpt',form='formatted',status='old')
!write(13,*) it, fz(isx,isz)
!!write(*,*) it,ist
!close(13)

  !NF for point source
  fx(isx,isz)=0.d0







  return
end subroutine calf2



subroutine calf( nx,nz,it,t,ist,isx,isz,dt,dx,dz,rho,tp,ts,fx,fz )
  implicit none
  real*8 pi
  parameter ( pi=3.1415926535897932d0 )
  integer nx,nz,it,ist,isx,isz
  real*8 t,dt,dx,dz,rho,tp,ts,fx(nx+1,nz+1),fz(nx+1,nz+1)
  real*8 b
  
  if ( it.le.ist ) then
     b = pi * ( t - ts ) / tp
     fx(isx,isz) &
          = dsqrt(pi) / 2.d0 * (b*b-0.5d0) * dexp(-b*b) &
          / ( dx * dz )
     fx(isx,isz) = fx(isx,isz) * dt * dt / rho
     fz(isx,isz) = fx(isx,isz)
     if ( (it.eq.0).or.(it.eq.ist) ) then
        fx(isx,isz) = fx(isx,isz) / 2.d0
        fz(isx,isz) = fz(isx,isz) / 2.d0
     endif
  else
     fx(isx,isz) = 0.d0
     fz(isx,isz) = 0.d0
  endif

  ! NF for point source
  !! NOT change fz and fx into different direction 28/02/2018
  fx(isx,isz)=0.d0
  fz(isx,isz)=fz(isx,isz)*1.d4   !Change for matching with specfem2d

!%! Plot source wavelet!


  return
end subroutine calf
