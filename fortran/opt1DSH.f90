program opt1d
  
  ! Computation of the synthetic seismograms in the time domain
  !using the normal operators.
  !
  !					originally from	1995.10  N.Takeuchi
  !                                                     2016. 5  N.Fuji
  
  implicit none

  integer, parameter :: maxnz = 5000
  double precision, parameter :: pi=3.1415926535897932d0 
  integer :: ir, iir
  !  parameters for the gridding
  integer nt,nz,it,ist,isz,i
  real*8 dt,dz
  ! parameter for the structure
  double precision rrho(4),kkappa(4)
  double precision rho(maxnz+1),kappa(maxnz+1)
  ! parameter for the wavefield
  double precision t
  double precision u(maxnz+1),u1(maxnz+1),u2(maxnz+1),f(maxnz+1)
  double precision e1(maxnz+1),e2(maxnz+1),e3(maxnz+1)
  double precision f1(maxnz+1),f2(maxnz+1),f3(maxnz+1)
  real utotal(0:50,0:2000)
  double precision uanalytic(1:50,0:2000),source(0:2000),xrec(1:50)
  double precision traveltime,coefamp
  ! parameter for the source
  double precision f0,t0,a
  double precision Courant_number
  ! parameter for the receiver
  integer nr
  ! switch OPT / CONV
  logical,parameter :: optimise = .true.
  
  rho=0.d0
  kappa=0.d0
  ! reading the parameter files
  call pinput( maxnz,nt,nz,dt,dz,rho,kappa,f0,t0,nr )
  
  !Initializing the data
  u=0.d0
  utotal=0.e0
  u1=0.d0
  u2=0.d0
  f=0.d0
  uanalytic=0.d0
  !call datainit( maxnz,u )
  !call datainit( maxnz,u1 )
  !call datainit( maxnz,u2 )
  !call datainit( maxnz,f )
  !call datainit( maxnz,rho )
  !call datainit( maxnz,kappa )
  ! computing the intermediate parameters
  !call calstruct( rrho,dz,nz,rho )
  !call calstruct( kkappa,dz,nz,kappa )
  call cales( nz,rho,kappa,dt,dz,e1,e2,e3,f1,f2,f3 )
  ist = nt/4
  isz = nz*3 / 4 + 1

    

  Courant_number = sqrt(kappa(isz)/rho(isz)) * dt /dz
  print *, 'Courant_number ', Courant_number
  do i=1,nz
     write(13,*) dz*dble(i),kappa(i),rho(i)
  enddo
  
  t=0.d0
  do it=0,nt

     
     call calf2(it,t,ist,isz,dt,dz,rho(isz),f0,t0,f)
     write(14,*) f(isz)
     source(it)=f(isz)
     t=t+dt
  enddo
  !stop

 
  do ir=1,50
     iir=nz/2+6*ir
     xrec(ir)=dble(iir)*dz-dble(isz)*dz
  enddo

  t = 0.d0
  write(15,*) real(t),real(u(nr))
  do it=0,nt
     call calf2( it,t,ist,isz,dt,dz,rho(isz),f0,t0,f )
     ! evaluating the next step
     call calstep( nz,e1,e2,e3,f1,f2,f3,u,u1,u2,f,optimise)
     ! increment of t
     !utotal(:,it)=u(:)
     t = t + dt
     !if ( mod(it,10).eq.9 ) write(6,*) real(t),real(u(nr))
     utotal(0,it)=real(t)
     do ir=1,50
        iir=nz/2+6*ir
        utotal(ir,it)=real(u(iir))
     enddo
     !write(15,*) real(t),real(u(nr))
  enddo
  
  ! analytic solution
  coefamp= 1.d0/2.d0*sqrt(kappa(isz)/rho(isz))*kappa(isz)
  a=pi*pi*f0*f0
  do ir=1,50
     iir=nz/2+6*ir
     traveltime=abs(xrec(ir))/(sqrt(kappa(isz)/rho(isz)))
     do it=0,nt
        if(dt*dble(it)>traveltime) then

           t=dble(it)*dt-traveltime-t0
           
           uanalytic(ir,it)= t*exp(-a*t**2)
           !print *,uanalytic(ir,it)
        endif
     enddo
  enddo

  open(2,file="1Danalytical_homo.dat",form="formatted")
  write(2,*) real(uanalytic(:,:))
  close(2)


  !if(1.eq.0) then
    ! OPEN(1,FILE="1Dsynthetic.dat",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=kind(0e0)*51*2001, &
    !      STATUS="REPLACE") ! gfortran defines RECL as bytes
     open(1,file="1Dsynthetic.dat",form="formatted")
    ! write(1,rec=nt) utotal(:,:)
     write(1,*)utotal(:,:)
     
     close(1)
  !endif


  !call create_color_image2(utotal,51,2001)


  
end program opt1d



subroutine pinput (maxnz,nt,nz,dt,dz,rho,kappa,f0,t0,nr )

  implicit none
  double precision f0,t0
  integer maxnz,nt,nx,nz,nrx,nrz,nr,nzz
  double precision dt,dx,dz !rho(*),lam(*),mu(*),tp,ts
  double precision rho(maxnz+1),kappa(maxnz+1),vp(maxnz+1),tmpvector(maxnz+1)
  double precision lam(maxnz+1),mu(maxnz+1),vs(maxnz+1)
  character*80 tmpfile,dummy
  character*80 vpfile,vsfile,rhofile
  tmpfile='tmpfileforwork'
  
  
  ! temporary file open
  open( unit=11, file=tmpfile, status='unknown' )
  !writing to the temporary file
100 continue
  read(5,110) dummy
110 format(a80)
  if ( dummy(1:1).eq.'c' ) goto 100
  if ( dummy(1:3).eq.'end' ) goto 120
  write(11,110) dummy
  goto 100
120 continue
  ! temporary file close
  close(11)
  
  ! temporary file open
  open( unit=11, file=tmpfile, status='unknown' )
  ! reading the parameter
  read(11,*) nt,nx,nz
  if ( nx.gt.maxnz ) pause 'nx is too large (pinput).'
  if ( nz.gt.maxnz ) pause 'nz is too large (pinput).'
  read(11,*) dt,dx,dz
111 format(a80)
  read(11,111) vpfile
  read(11,111) vsfile
  read(11,111) rhofile
  read(11,*) f0,t0
  read(11,*) nrx,nrz
  ! temporary file close
  close(11)
  !
 
  rho=0.d0
  vp=0.d0
  vs=0.d0
  lam=0.d0
  mu=0.d0


  call calstruct2D1D(maxnz,rhofile,dx,dz,nx,nz,rho)
  call calstruct2D1D(maxnz,vpfile,dx,dz,nx,nz,vp)
  call calstruct2D1D(maxnz,vsfile,dx,dz,nx,nz,vs)

  call calstruct2(maxnz,nz,rho,vp,vs,lam,mu)
  

  nzz=nz*2

  ! mirroring

  tmpvector(1:nz)=rho(1:nz)
  tmpvector(nz+1:nzz)=rho(nz)
  rho(:)=0.d0
  rho(maxnz/2-nzz+1:maxnz/2)=tmpvector(nzz:1:-1)
  rho(maxnz/2+1:maxnz/2+nzz)=tmpvector(1:nzz)
  
  tmpvector(1:nz)=lam(1:nz)
  tmpvector(nz+1:nzz)=lam(nz)
  lam(:)=0.d0
  lam(maxnz/2-nzz+1:maxnz/2)=tmpvector(nzz:1:-1)
  lam(maxnz/2+1:maxnz/2+nzz)=tmpvector(1:nzz)

  tmpvector(1:nz)=mu(1:nz)
  tmpvector(nz+1:nzz)=mu(nz)
  mu(:)=0.d0
  mu(maxnz/2-nzz+1:maxnz/2)=tmpvector(nzz:1:-1)
  mu(maxnz/2+1:maxnz/2+nzz)=tmpvector(1:nzz)

  nz=nzz

  tmpvector(1:2*nz)=rho(maxnz/2-nz+1:maxnz/2+nz)
  rho(:)=0.d0
  rho(1:2*nz)=tmpvector(1:2*nz)

  tmpvector(1:2*nz)=lam(maxnz/2-nz+1:maxnz/2+nz)
  lam(:)=0.d0
  lam(1:2*nz)=tmpvector(1:2*nz)
  
  tmpvector(1:2*nz)=mu(maxnz/2-nz+1:maxnz/2+nz)
  mu(:)=0.d0
  mu(1:2*nz)=tmpvector(1:2*nz)
    
  kappa(1:2*nz)=2.d0*mu(1:2*nz)+lam(1:2*nz)

  nz=2*nz-1
  nr=3*nz/4

  ! homogeneous case

  tmpvector(1:nz)=rho(1:nz)
  rho(1:nz)=tmpvector(nz/2)
  
  tmpvector(1:nz)=kappa(1:nz)
  kappa(1:nz)=tmpvector(nz/2)

  return
end subroutine pinput


subroutine calstruct2(maxnz,nz,rho,vp,vs,lam,mu)
  implicit none
  
  integer i,j,maxnz,nz
  double precision rho(maxnz+1),vp(maxnz+1),vs(maxnz+1)
  double precision lam(maxnz+1),mu(maxnz+1)

 
  do i=1,nz+1
     mu(i)=rho(i)*vs(i)*vs(i)
     lam(i)=rho(i)*vp(i)*vp(i)-2*mu(i)
  enddo

end subroutine calstruct2


subroutine calstruct2D1D( maxnz,file2d,dx,dz,nx,nz,rho )
  implicit none
  integer maxnz,nx,nz
  double precision dx,dz,rho(maxnz+1)
  real(kind(1.0)) rrho(maxnz+1,maxnz+1)
  integer i,j,k,nox(6),noz(6)
  double precision x,z,xmax,zmax,trho,coef1,coef2
  integer recl_size
  character*80 file2d
  recl_size=kind(1.0)*(nx+1)*(nz+1)

  open (1,file=file2d,form='unformatted',access='direct',recl=recl_size)
  read(1,rec=1) rrho(1:nx+1,1:nz+1)
  close(1)
  rho(1:nz+1)=1.d-3*rrho(1,1:nz+1)
  
  return
end subroutine calstruct2D1D


subroutine pinput_old( maxnz,nt,nz,dt,dz,rho,kappa,tp,ts,nr )
  
  integer maxnz,nt,nz,nr
  double precision dt,dz,rho(4),kappa(4),tp,ts
  character*80 tmpfile,dummy
  
  data tmpfile / '/tmp/work' /
  
  ! temporary file open
  open( unit=11, file=tmpfile, status='unknown' )
  !writing to the temporary file
100 continue
  read(5,110) dummy
110 format(a80)
  if ( dummy(1:1).eq.'c' ) goto 100
  if ( dummy(1:3).eq.'end' ) goto 120
  write(11,*) dummy
  goto 100
120 continue
  ! temporary file close
  close(11)
  ! 
  ! temporary file open
  open( unit=11, file=tmpfile, status='unknown' )
  ! reading the parameter
  read(11,*) nt,nz
  if ( nz.gt.maxnz ) pause 'nz is too large (pinput).'
  read(11,*) dt,dz
  read(11,*) rho(1),rho(2),rho(3),rho(4)
  read(11,*) kappa(1),kappa(2),kappa(3),kappa(4)
  read(11,*) tp,ts
  read(11,*) nr
  ! temporary file close
  close(11)
  
  return
end subroutine pinput_old


subroutine datainit( maxnz,u )
  integer maxnz
  double precision u(*)
  integer i
  
  do i=1,maxnz+1
     u(i) = 0.d0
  enddo
  return
end subroutine datainit


subroutine calstruct( rrho,dz,nz,rho )
  integer nz
  double precision rrho(4),dz,rho(*)
  integer i,j
  double precision r,rmax,coef,trho
  
  rmax = dz * nz
  do i=1,nz+1
     r = dble(i-1) * dz
     trho = 0.d0
     do j=1,4
        if ( j.eq.1 ) then
           coef = 1.d0
        else
           coef = coef * ( r / rmax )
        endif
        trho = trho + rrho(j) * coef
     enddo
     rho(i) = trho
  enddo
  return
end subroutine calstruct



subroutine cales( nz,rho,kappa,dt,dz,e1,e2,e3,f1,f2,f3 )
  integer nz
  double precision rho(*),kappa(*),dt,dz,e1(*),e2(*),e3(*),f1(*),f2(*),f3(*)
  integer iz
  double precision dt2,dz2
  
  dt2 = dt * dt
  dz2 = dz * dz
  
  e1(1) = 0.d0
  e2(1) = 2.d0 - ( kappa(1) + kappa(2) ) &
       / rho(1) * dt2 / dz2
  e3(1) = ( kappa(1) + kappa(2) ) / rho(1) * dt2 / dz2
  f1(1) = 0.d0
  f2(1) = e2(1) / 12.d0
  f3(1) = ( e3(1) - 2.d0 ) / 12.d0
  do iz=2,nz
     e1(iz) = ( kappa(iz-1) + kappa(iz) ) &
          / ( 2.d0 * rho(iz) ) * dt2 / dz2
     e2(iz) = 2.d0 &
          - ( kappa(iz-1) + 2.d0 * kappa(iz) + kappa(iz+1) ) &
          / ( 2.d0 * rho(iz) ) * dt2 / dz2
     e3(iz) = ( kappa(iz) + kappa(iz+1) ) &
          / ( 2.d0 * rho(iz) ) * dt2 / dz2
     f1(iz) = ( e1(iz) - 1.d0 ) / 12.d0
     f2(iz) = e2(iz) / 12.d0
     f3(iz) = ( e3(iz) - 1.d0 ) / 12.d0
  enddo
  e1(nz+1) = ( kappa(nz) + kappa(nz+1) ) &
       / rho(nz+1) * dt2 / dz2
  e2(nz+1) = 2.d0 - ( kappa(nz) + kappa(nz+1) ) &
       / rho(nz+1) * dt2 / dz2
  e3(nz+1) = 0.d0
  f1(nz+1) = ( e1(nz+1) - 2.d0 ) / 12.d0
  f2(nz+1) = e2(nz+1) / 12.d0
  f3(nz+1) = 0.d0
  
  return
end subroutine cales


subroutine calf( it,t,ist,isz,dt,dz,rho,tp,ts,f )
  double precision pi
  parameter ( pi=3.1415926535897932d0 )
  integer it,ist,isz
  double precision t,dt,dz,rho,tp,ts,f(*)
  double precision b
  
  if ( it.le.ist ) then
     b = pi * ( t - ts ) / tp
     f(isz) = dsqrt(pi) / 2.d0 * (b*b-0.5d0) * dexp(-b*b) / dz
     if ( (it.eq.0).or.(it.eq.ist) ) f(isz) = f(isz) / 2.d0
     f(isz) = f(isz) * dt * dt / rho
  else
     f(isz) = 0.d0
  endif

  return
end subroutine calf


subroutine calf2( it,t,ist,isz,dt,dz,rho,f0,t0,f )
  implicit none
  double precision pi
  parameter ( pi=3.1415926535897932d0 )
  integer it,ist,isz
  double precision t,dt,dz,rho,f0,t0,f(*)
  double precision b,a,factor

  factor=1.d3
  
  if ( it.le.ist ) then

     a = pi*pi*f0*f0*(t-t0)*(t-t0)
     !print *,pi, t,t0,t-t0,a
     !print *, f0,t0,a
     ! Ricker source time function (second derivative of a Gaussian)
     !print *,exp(-a)
     f(isz)=factor*(1.d0-2.d0*a)*exp(-a)
    ! f(isz) = factor * (1.d0 - 2.d0*a)*exp(-a);

     f(isz) = f(isz) * dt * dt / rho
    
     if ( (it.eq.0).or.(it.eq.ist) ) then
        f(isz) = f(isz) / 2.d0
        
     endif
  else
     f(isz) = 0.d0
  endif
  
  
end subroutine calf2


subroutine calstep( nz,e1,e2,e3,f1,f2,f3,u,u1,u2,f,optimise )
  integer nz
  double precision e1(*),e2(*),e3(*),f1(*),f2(*),f3(*)
  double precision u(*),u1(*),u2(*),f(*)
  integer iz
  double precision tmp1,tmp2,tmp3
  logical optimise
  
  ! evaluating the u using the unmodified operators
  u(1) = - u2(1) &
       + e2(1) * u1(1) &
       + e3(1) * u1(2) &
       + f(1)
  do iz=2,nz
     u(iz) = - u2(iz) &
          + e1(iz) * u1(iz-1) &
          + e2(iz) * u1(iz) &
          + e3(iz) * u1(iz+1) &
          + f(iz)
  enddo
  u(nz+1) = - u2(nz+1) &
       + e1(nz+1) * u1(nz) &
       + e2(nz+1) * u1(nz+1) &
       + f(nz+1)

  if(optimise) then
     
     ! computing du using 1-D Born approximation
     tmp2 = ( u(1) - u1(1) - u1(1) + u2(1) )
     do iz=1,nz
        tmp3 = ( u(iz+1) - u1(iz+1) - u1(iz+1) + u2(iz+1) )
        u(iz) = u(iz) &
             + f1(iz) * tmp1 &
             + f2(iz) * tmp2 &
             + f3(iz) * tmp3
        tmp1 = tmp2
        tmp2 = tmp3
     enddo
     u(nz+1) = u(nz+1) &
          + f1(nz+1) * tmp1 &
          + f2(nz+1) * tmp2
  endif
  ! swapping u1 & u2
  do iz=1,nz+1
     u2(iz) = u1(iz)
     u1(iz) = u(iz)
  enddo
  
  return
end subroutine calstep



subroutine create_color_image2(image_data_2D,NX,NY)
	
  implicit none
  
  !       non linear display to enhance small amplitudes for graphics
  double precision, parameter :: POWER_DISPLAY = 0.30d0
  
  !       amplitude threshold above which we draw the color point
  double precision, parameter :: cutvect = 0.01d0
  
  !       use black or white background for points that are below the threshold
  logical, parameter :: WHITE_BACKGROUND = .true.
  
  !       size of cross and square in pixels drawn to represent the source and the receivers
  integer, parameter :: width_cross=5,thickness_cross=1
  integer, parameter :: size_square=3
  
  integer NX,NY,it,field_number,ISOURCE,JSOURCE,NPOINTS_PML
  logical USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX
  integer, parameter :: nrec=1
  double precision, dimension(NX,NY) :: image_data_2D
  
  
  integer, dimension(nrec) :: ix_rec,iy_rec
  
  integer :: ix,iy,irec,ii
  
  character(len=100) :: file_name,system_command1,system_command2,system_command3
	
  integer :: R, G, B
  
  double precision :: normalized_value,max_amplitude

  it = 999
  field_number = 1
  !       open image file and create system command to convert image to more convenient format
  !       use the "convert" command from ImageMagick http://www.imagemagick.org
  if(field_number == 1) then
     write(file_name,"('image',i6.6,'_Ux.pnm')") it
     write(system_command1, "('convert image',i6.6,'_Ux.pnm synthe1D',i6.6,'.png')") it,it
     write(system_command2, "('rm image',i6.6,'_Ux.pnm')") it
  else if(field_number == 2) then
  !   write(file_name,"('image',i6.6,'_Uz.pnm')") it
  !   write(system_command1,"('convert image',i6.6,'_Uz.pnm snapshots/imageUz',i6.6,'.png')") it,it
  !   write(system_command2,"('rm image',i6.6,'_Uz.pnm')") it
  endif
  
  open(unit=27, file=file_name, status='unknown')
  
  write(27,"('P3')")	! write image in PNM P3 format
  
  write(27,*) NX*100,NY	! write image size
  write(27,*) '255'	! maximum value of each pixel color
	
  !       compute maximum amplitude
  max_amplitude = maxval(abs(image_data_2D))
  
  !       image starts in upper-left corner in PNM format
  do iy=1,NY,1
     do ix=1,NX
        
        !       define data as vector component normalized to [-1:1] and rounded to nearest integer
        !       keeping in mind that amplitude can be negative
        normalized_value = image_data_2D(ix,iy) / max_amplitude /2.d0
        
        !       suppress values that are outside [-1:+1] to avoid small edge effects
        if(normalized_value < -1.d0) normalized_value = -1.d0
        if(normalized_value > 1.d0) normalized_value = 1.d0
        
       
           !       represent regular image points using red if value is positive, blue if negative
        if(normalized_value >= 0.d0) then
           R = nint(255.d0*normalized_value**POWER_DISPLAY)
           G = 0
           B = 0
        else
           R = 0
           G = 0
           B = nint(255.d0*abs(normalized_value)**POWER_DISPLAY)
        endif
        
          
     
        !       write color pixel
        do ii=1,100
        write(27,"(i3,' ',i3,' ',i3)") R,G,B
        enddo
     enddo
  enddo
  
  !       close file
close(27)

!       call the system to convert image to GIF (can be commented out if "call system" is missing in your compiler)
call system(system_command1)
call system(system_command2)
end subroutine create_color_image2
