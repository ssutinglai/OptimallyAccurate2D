program opt22
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Computation of the synthetic seismograms in the time domain
  ! using the normal operators.
  ! --- heterogeneous medium
  !
  !						1997.6  N.Takeuchi
  !                                               2016.5. N.Fuji
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer maxnz
  parameter ( maxnz = 1000 )
  !
  ! parameters for the gridding
  real*8 dt,dx,dz
  ! parameters for the wavefield
  integer nt,nx,nz,it,ist,isx,isz
  real*8 ux(maxnz+1,maxnz+1),uz(maxnz+1,maxnz+1)
  real*8 ux1(maxnz+1,maxnz+1),ux2(maxnz+1,maxnz+1)
  real*8 uz1(maxnz+1,maxnz+1),uz2(maxnz+1,maxnz+1)
  real*8  e1(maxnz+1,maxnz+1), e2(maxnz+1,maxnz+1)
  real*8  e3(maxnz+1,maxnz+1), e4(maxnz+1,maxnz+1)
  real*8  e5(maxnz+1,maxnz+1), e6(maxnz+1,maxnz+1)
  real*8  e7(maxnz+1,maxnz+1), e8(maxnz+1,maxnz+1)
  real*8 e13(maxnz+1,maxnz+1),e14(maxnz+1,maxnz+1)
  real*8 e15(maxnz+1,maxnz+1),e16(maxnz+1,maxnz+1)
  real*8 e17(maxnz+1,maxnz+1),e18(maxnz+1,maxnz+1)
  real*8 e19(maxnz+1,maxnz+1),e20(maxnz+1,maxnz+1)
  real*8  f1(maxnz+1,maxnz+1), f2(maxnz+1,maxnz+1)
  real*8  f3(maxnz+1,maxnz+1), f4(maxnz+1,maxnz+1)
  real*8  f5(maxnz+1,maxnz+1), f6(maxnz+1,maxnz+1)
  real*8  f7(maxnz+1,maxnz+1), f8(maxnz+1,maxnz+1)
  real*8 f13(maxnz+1,maxnz+1),f14(maxnz+1,maxnz+1)
  real*8 f15(maxnz+1,maxnz+1),f16(maxnz+1,maxnz+1)
  real*8 f17(maxnz+1,maxnz+1),f18(maxnz+1,maxnz+1)
  real*8 f19(maxnz+1,maxnz+1),f20(maxnz+1,maxnz+1)
  real*8 work(maxnz+1,32)
  ! parameter for the structure
  !real*8 rrho(6),llam(6),mmu(6)
  character(80) :: vpfile, vsfile, rhofile
  real*8 rho(maxnz+1,maxnz+1)
  real*8 lam(maxnz+1,maxnz+1),mu(maxnz+1,maxnz+1)
  ! parameter for the source
  real*8 tp,ts
  ! parameter for the receiver
  integer nrx,nrz
  ! parameter for the waveform
  real*8 t
  !parameter for video
  real,dimension(maxnz+1,maxnz+1) :: snapux,snapuz
  integer, parameter :: IT_DISPLAY = 10
  integer :: ix_rec(1),iz_rec(1)
  integer(2) head(1:120)
  character(80) :: routine
  
  logical,parameter :: dummylog = .false.
  ! switch
  logical,parameter :: optimise = .true.
  
  ! reading the parameter files
  call pinput( maxnz,nt,nx,nz,dt,dx,dz,rrho,llam,mmu,tp,ts,nrx,nrz )

  ! Initializing the data
  call datainit( maxnz,maxnz,ux )
  call datainit( maxnz,maxnz,uz )
  call datainit( maxnz,maxnz,ux1 )
  call datainit( maxnz,maxnz,uz1 )
  call datainit( maxnz,maxnz,ux2 )
  call datainit( maxnz,maxnz,uz2 )
  call datainit( maxnz,maxnz,rho )
  call datainit( maxnz,maxnz,lam )
  call datainit( maxnz,maxnz,mu )
  call datainit( maxnz,31,work )
  !computing the intermediate parameters
  call calstruct( maxnz,rhofile,dx,dz,nx,nz,rho )
  call calstruct( maxnz,vpfile,dx,dz,nx,nz,lam )
  call calstruct( maxnz,vsfile,dx,dz,nx,nz,mu )
  call cales( maxnz,nx,nz,rho,lam,mu,dt,dx,dz, &
       e1, e2, e3, e4, e5, e6, e7, e8, &
       e13,e14,e15,e16,e17,e18,e19,e20, &
       f1, f2, f3, f4, f5, f6, f7, f8, &
       f13,f14,f15,f16,f17,f18,f19,f20 )
  
  call datainit( maxnz,maxnz,lam )
  call datainit( maxnz,maxnz,mu )
  ist = dnint( 2 * tp / dt )
  isx = nx / 2 + 1
  isz = nz / 2 + 1
  
  t=0.d0
  do it=0,nt
     call calf( maxnz,it,t,ist,isx,isz,dt,dx,dz,rho(isx,isz),tp,ts,lam,mu )
    
     write(13,*) t, lam(isx,isz),mu(isx,isz)
     t=t+dt
  enddo
  !stop
  
  
  t = 0.d0
  !write(6,*) real(t),real(ux(nrx,nrz)),real(uz(nrx,nrz))
  do it=0,nt
     call calf( maxnz,it,t,ist,isx,isz,dt,dx,dz,rho(isx,isz),tp,ts,lam,mu )
     ! evaluating the next step
     call calstep( maxnz,nx,nz, &
          e1, e2, e3, e4, e5, e6, e7, e8, &
          e13,e14,e15,e16,e17,e18,e19,e20, &
          f1, f2, f3, f4, f5, f6, f7, f8, &
          f13,f14,f15,f16,f17,f18,f19,f20, &
          ux,uz,ux1,ux2,uz1,uz2,isx,isz,lam,mu, &
          work(1,1), work(1,5), work(1,9),work(1,13), &
          work(1,17),work(1,18),work(1,20),work(1,21), &
          work(1,23),work(1,24),work(1,28),work(1,29), optimise)
     ! increment of t
     t = t + dt
     !write(6,*) real(t),real(ux(nrx,nrz)),real(uz(nrx,nrz))
     
     
     write(*,*) it, ' of ', nt
     if(mod(it,IT_DISPLAY) == 0)then
        !
        !head=0
        !head(58) = nx
        !head(59) = dz * 1E3
        !snapux=0.e0
        !snapuz=0.e0
        !snapux(1:nx,1:nz) = ux(1:nx,1:nz)
        !snapuz(1:nx,1:nz) = uz(1:nx,1:nz)
        !write(routine,'(a12,i5.5,a9)') './snapshots/',it,'snapUx.su'
        !open(21,file=routine,access='stream')
        !do j = 1,nx,1
        !   write(21) head,(real(snapux(k,j)),k=1,nz)
        !enddo
        !close(21)
        !write(routine,'(a12,i5.5,a9)') './snapshots/',it,'snapUz.su'
        !open(21,file=routine,access='stream')
        !do j = 1,nx,1
        !   write(21) head,(real(snapuz(k,j)),k=1,nz)
        !enddo
        !close(21)
        
        ix_rec(1)=nrx
        iz_rec(1)=nrz
        
        !call create_color_image(ux(1:nx,1:nz),nx,nz,it,isx,isz,ix_rec,iz_rec,1,0, &
        !     dummylog,dummylog,dummylog,dummylog,1)
        call create_color_image(uz(1:nx+1,1:nz+1),nx+1,nz+1,it,isx,isz,ix_rec,iz_rec,1,0, &
             dummylog,dummylog,dummylog,dummylog,2)
        
        
        
        
        
     endif
  enddo


  

  !
end program opt22


subroutine pinput( maxnz,nt,nx,nz,dt,dx,dz,vpfile,vsfile,rhofile,tp,ts,nrx,nrz )
  
  integer maxnz,nt,nx,nz,nrx,nrz
  real*8 dt,dx,dz !rho(*),lam(*),mu(*),tp,ts
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
  write(11,*) dummy
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
  read(11,*) vpfile
  read(11,*) vsfile
  read(11,*) rhofile
  read(11,*) tp,ts
  read(11,*) nrx,nrz
  ! temporary file close
  close(11)
  !
  return
end subroutine pinput


subroutine datainit( nx,nz,ux )
  
  integer nx,nz
  real*8 ux(nx+1,*)
  integer i,j
  
  do j=1,nz+1
     do i=1,nx+1
        ux(i,j) = 0.d0
     enddo
  enddo
  
  return
end subroutine datainit


subroutine calstruct( maxnz,file2d,dx,dz,nx,nz,rho )

  integer maxnz,nx,nz
  real*8 dx,dz,rho(maxnz+1,*)
  integer i,j,k,nox(6),noz(6)
  real*8 x,z,xmax,zmax,trho,coef1,coef2
  character*80 file2d


  open (1,file=file2d,form='unformatted',access='direct',recl=recl_size)
  read(1,rec=1) rho(1:nx+1,1:nz+1)
  close(1)
  
  return
end subroutine calstruct


subroutine cales( maxnz,nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20 )
  
  integer maxnz,nx,nz
  real*8 rho(maxnz+1,*),lam(maxnz+1,*),mu(maxnz+1,*)
  real*8 dt,dx,dz
  real*8  e1(maxnz+1,*), e2(maxnz+1,*), e3(maxnz+1,*)
  real*8  e4(maxnz+1,*), e5(maxnz+1,*), e6(maxnz+1,*)
  real*8  e7(maxnz+1,*), e8(maxnz+1,*)
  real*8 e13(maxnz+1,*),e14(maxnz+1,*),e15(maxnz+1,*)
  real*8 e16(maxnz+1,*),e17(maxnz+1,*),e18(maxnz+1,*)
  real*8 e19(maxnz+1,*),e20(maxnz+1,*)
  real*8  f1(maxnz+1,*), f2(maxnz+1,*), f3(maxnz+1,*)
  real*8  f4(maxnz+1,*), f5(maxnz+1,*), f6(maxnz+1,*)
  real*8  f7(maxnz+1,*), f8(maxnz+1,*)
  real*8 f13(maxnz+1,*),f14(maxnz+1,*),f15(maxnz+1,*)
  real*8 f16(maxnz+1,*),f17(maxnz+1,*),f18(maxnz+1,*)
  real*8 f19(maxnz+1,*),f20(maxnz+1,*)
  integer ix,iz
  real*8 dt2,dx2,dz2,dxdz
  
  dt2 = dt * dt
  dx2 = dx * dx
  dz2 = dz * dz
  dxdz = dx * dz
  
  do iz=2,nz
     do ix=2,nx
        e1(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix-1,iz) + lam(ix,iz) ) &
             + 2.d0 * ( mu(ix-1,iz) + mu(ix,iz) ) ) &
             / ( 2.d0 * dx2 )
        e2(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix,iz) + lam(ix+1,iz) ) &
             + 2.d0 * ( mu(ix,iz) + mu(ix+1,iz) ) ) &
             / ( 2.d0 * dx2 )
        e3(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix,iz-1) + mu(ix,iz) ) &
             / ( 2.d0 * dz2 )
        e4(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix,iz) + mu(ix,iz+1) ) &
             / ( 2.d0 * dz2 )
        e5(ix,iz) = dt2 / rho(ix,iz) * lam(ix-1,iz) &
             / ( 4.d0 * dxdz )
        e6(ix,iz) = dt2 / rho(ix,iz) * lam(ix+1,iz) &
             / ( 4.d0 * dxdz )
        e7(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz-1) &
             / ( 4.d0 * dxdz )
        e8(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+1) &
             / ( 4.d0 * dxdz )
        e13(ix,iz) = dt2 / rho(ix,iz) * lam(ix-1,iz) &
             * ( -5.d0 ) / ( 1728.d0 * dxdz )
        e14(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
             * ( -3.d0 ) / ( 1728.d0 * dxdz )
        e15(ix,iz) = dt2 / rho(ix,iz) * lam(ix+1,iz) &
           * ( +9.d0 ) / ( 1728.d0 * dxdz )
        if ( ix+2.le.nx+1 ) then
           e16(ix,iz) = dt2 / rho(ix,iz) * lam(ix+2,iz) &
                * ( -1.d0 ) / ( 1728.d0 * dxdz )
        else
           e16(ix,iz) = 0.d0
        endif
        e17(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz-1) &
             * ( -5.d0 ) / ( 1728.d0 * dxdz )
        e18(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
             * ( -3.d0 ) / ( 1728.d0 * dxdz )
        e19(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+1) &
             * ( +9.d0 ) / ( 1728.d0 * dxdz )
        if ( iz+2.le.nz+1) then
           e20(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+2) &
                * ( -1.d0 ) / ( 1728.d0 * dxdz )
        else
           e20(ix,iz) = 0.d0
        endif
        f1(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix-1,iz) + mu(ix,iz) ) &
             / ( 2.d0 * dx2 )
        f2(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix,iz) + mu(ix+1,iz) ) &
             / ( 2.d0 * dx2 )
        f3(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix,iz-1) + lam(ix,iz) ) &
             + 2.d0 * ( mu(ix,iz-1) + mu(ix,iz) ) ) &
             / ( 2.d0 * dz2 )
        f4(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix,iz) + lam(ix,iz+1) ) &
           + 2.d0 * ( mu(ix,iz) + mu(ix,iz+1) ) ) &
           / ( 2.d0 * dz2 )
        f5(ix,iz) = dt2 / rho(ix,iz) * mu(ix-1,iz) &
             / ( 4.d0 * dxdz )
        f6(ix,iz) = dt2 / rho(ix,iz) * mu(ix+1,iz) &
             / ( 4.d0 * dxdz )
        f7(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-1) &
             / ( 4.d0 * dxdz )
        f8(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz+1) &
             / ( 4.d0 * dxdz )
        if ( ix-2.ge.1 ) then
           f13(ix,iz) = dt2 / rho(ix,iz) * mu(ix-2,iz) &
                * (  1.d0 ) / ( 1728.d0 * dxdz )
        else
           f13(ix,iz) = 0.d0
        endif
        f14(ix,iz) = dt2 / rho(ix,iz) * mu(ix-1,iz) &
             * ( -9.d0 ) / ( 1728.d0 * dxdz )
        f15(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
             * (  3.d0 ) / ( 1728.d0 * dxdz )
        f16(ix,iz) = dt2 / rho(ix,iz) * mu(ix+1,iz) &
             * (  5.d0 ) / ( 1728.d0 * dxdz )
        if ( iz-2.ge.1 ) then
           f17(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-2) &
                * (  1.d0 ) / ( 1728.d0 * dxdz )
        else
           f17(ix,iz) = 0.d0
        endif
        f18(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-1) &
             * ( -9.d0 ) / ( 1728.d0 * dxdz )
        f19(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
             * (  3.d0 ) / ( 1728.d0 * dxdz )
        f20(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz+1) &
             * (  5.d0 ) / ( 1728.d0 * dxdz )
     enddo
  enddo
  
  return
end subroutine cales


subroutine calf( maxnz,it,t,ist,isx,isz,dt,dx,dz,rho,tp,ts,fx,fz )

  real*8 pi
  parameter ( pi=3.1415926535897932d0 )
  integer maxnz,it,ist,isx,isz
  real*8 t,dt,dx,dz,rho,tp,ts,fx(maxnz+1,*),fz(maxnz+1,*)
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
  fz(isx,isx)=0.d0
  
  return
end subroutine calf




subroutine calstep( maxnz,nx,nz, &
     e1, e2, e3, e4, e5, e6, e7, e8, &
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, &
     ux,uz,ux1,ux2,uz1,uz2,isx,isz,fx,fz, &
     work1,work2,work3,work4, &
     work5,work6,work7,work8, &
     work9,work10,work11,work12,optimise )

  integer maxnz,nx,nz,isx,isz
  real*8 ux(maxnz+1,*),ux1(maxnz+1,*),ux2(maxnz+1,*)
  real*8 uz(maxnz+1,*),uz1(maxnz+1,*),uz2(maxnz+1,*)
  real*8  e1(maxnz+1,*), e2(maxnz+1,*), e3(maxnz+1,*)
  real*8  e4(maxnz+1,*), e5(maxnz+1,*), e6(maxnz+1,*)
  real*8  e7(maxnz+1,*), e8(maxnz+1,*)
  real*8 e13(maxnz+1,*),e14(maxnz+1,*),e15(maxnz+1,*)
  real*8 e16(maxnz+1,*),e17(maxnz+1,*),e18(maxnz+1,*)
  real*8 e19(maxnz+1,*),e20(maxnz+1,*)
  real*8  f1(maxnz+1,*), f2(maxnz+1,*), f3(maxnz+1,*)
  real*8  f4(maxnz+1,*), f5(maxnz+1,*), f6(maxnz+1,*)
  real*8  f7(maxnz+1,*), f8(maxnz+1,*)
  real*8 f13(maxnz+1,*),f14(maxnz+1,*),f15(maxnz+1,*)
  real*8 f16(maxnz+1,*),f17(maxnz+1,*),f18(maxnz+1,*)
  real*8 f19(maxnz+1,*),f20(maxnz+1,*)
  real*8 fx(maxnz+1,*),fz(maxnz+1,*)
  real*8 work1(maxnz+1,-2:1),work2(maxnz+1,-1:2)
  real*8 work3(maxnz+1,-2:1),work4(maxnz+1,-1:2)
  real*8 work5(*),work6(maxnz+1,0:1)
  real*8 work7(*),work8(maxnz+1,0:1)
  real*8 work9(*),work10(maxnz+1,-2:1)
  real*8 work11(*),work12(maxnz+1,-1:2)
  integer ix,iz,iz1,iz2,ix11,ix12,ix21,ix22
  logical optimise


  ! predicting the wavefield
 do iz=2,nz
    do ix=2,nx
       ux(ix,iz) = 2.d0 * ux1(ix,iz) - ux2(ix,iz) &
            + e1(ix,iz) * ( ux1(ix-1,iz) - ux1(ix,iz) ) &
            + e2(ix,iz) * ( ux1(ix+1,iz) - ux1(ix,iz) ) &
            + e3(ix,iz) * ( ux1(ix,iz-1) - ux1(ix,iz) ) &
            + e4(ix,iz) * ( ux1(ix,iz+1) - ux1(ix,iz) ) &
            - e5(ix,iz) * ( uz1(ix-1,iz+1) - uz1(ix-1,iz-1) ) &
            + e6(ix,iz) * ( uz1(ix+1,iz+1) - uz1(ix+1,iz-1) ) &
            - e7(ix,iz) * ( uz1(ix+1,iz-1) - uz1(ix-1,iz-1) ) &
            + e8(ix,iz) * ( uz1(ix+1,iz+1) - uz1(ix-1,iz+1) )
       uz(ix,iz) = 2.d0 * uz1(ix,iz) - uz2(ix,iz) &
            + f1(ix,iz) * ( uz1(ix-1,iz) - uz1(ix,iz) ) &
            + f2(ix,iz) * ( uz1(ix+1,iz) - uz1(ix,iz) ) &
            + f3(ix,iz) * ( uz1(ix,iz-1) - uz1(ix,iz) ) &
            + f4(ix,iz) * ( uz1(ix,iz+1) - uz1(ix,iz) ) &
            - f5(ix,iz) * ( ux1(ix-1,iz+1) - ux1(ix-1,iz-1) ) &
            + f6(ix,iz) * ( ux1(ix+1,iz+1) - ux1(ix+1,iz-1) ) &
            - f7(ix,iz) * ( ux1(ix+1,iz-1) - ux1(ix-1,iz-1) ) &
            + f8(ix,iz) * ( ux1(ix+1,iz+1) - ux1(ix-1,iz+1) )
    enddo
 enddo
 ux(isx,isz) = ux(isx,isz) + fx(isx,isz)
 uz(isx,isz) = uz(isx,isz) + fz(isx,isz)

 
 if(optimise) then
    ! correcting the wavefield
    !
    do ix=2,nx
       iz1 = 2
       iz2 = 3
       work1(ix,-2) = 0.d0
       work1(ix,-1) = 0.d0
       work1(ix,0) = 0.d0
       work1(ix,1) = ux(ix,iz1) - 2.d0 * ux1(ix,iz1) + ux2(ix,iz1)
       work2(ix,-1) = 0.d0
       work2(ix,0) = 0.d0
       work2(ix,1) = uz(ix,iz1) - 2.d0 * uz1(ix,iz1) + uz2(ix,iz1)
       work2(ix,2) = uz(ix,iz2) - 2.d0 * uz1(ix,iz2) + uz2(ix,iz2)
       work3(ix,-2) = 0.d0
       work3(ix,-1) = 0.d0
       work3(ix,0) = 0.d0
       work3(ix,1) = work1(ix,1) + 12.d0 * ux1(ix,iz1)
       work4(ix,-1) = 0.d0
       work4(ix,0) = 0.d0
       work4(ix,1) = work2(ix,1) + 12.d0 * uz1(ix,iz1)
       work4(ix,2) = work2(ix,2) + 12.d0 * uz1(ix,iz2)
    enddo
    
    do ix=1,nx+1
       ix11 = max0( ix-1,1 )
       ix12 = min0( ix+1,nx+1 )
       ix21 = max0( ix-2,1 )
       ix22 = min0( ix+2,nx+1 )
       work6(ix,0) = 0.d0
       work6(ix,1) = &
            (           ( -work3(ix11,1) ) &
            + 10.d0 * ( -work3(ix,  1) ) & 
            +         ( -work3(ix12,1) ) &
            )
       work8(ix,0) = 0.d0
       work8(ix,1) = &
            (           ( -work4(ix11,1) ) &
            + 10.d0 * ( -work4(  ix,1) ) &
            +         ( -work4(ix12,1) ) &
            )
       work10(ix,-2) = 0.d0
       work10(ix,-1) = 0.d0
       work10(ix,0) = 0.d0
       work10(ix,1)   = (          work3(ix21,1) - 9.d0 * work3(ix11,1) &
            + 3.d0 * work3(  ix,1) + 5.d0 * work3(ix12,1) )
       work12(ix,-1) = 0.d0
       work12(ix,0) = 0.d0
       work12(ix,1) = ( - 5.d0 * work4(ix11,1) - 3.d0 * work4(  ix,1) &
            + 9.d0 * work4(ix12,1) -        work4(ix22,1) )
       work12(ix,2) = ( - 5.d0 * work4(ix11,2) - 3.d0 * work4(  ix,2) &
            + 9.d0 * work4(ix12,2) -        work4(ix22,2) )
       
    enddo
    
    do iz=2,nz
       iz1 = iz + 1
       iz2 = min0( iz+2, nz+1 )
       do ix=2,nx
          work1(ix,-2) = work1(ix,-1)
          work1(ix,-1) = work1(ix,0)
          work1(ix,0) = work1(ix,1)
          work1(ix,1) = ux(ix,iz1) - 2.d0 * ux1(ix,iz1) + ux2(ix,iz1)
          work2(ix,-1) = work2(ix,0)
          work2(ix,0) = work2(ix,1)
          work2(ix,1) = work2(ix,2)
          work2(ix,2) = uz(ix,iz2) - 2.d0 * uz1(ix,iz2) + uz2(ix,iz2)
          work3(ix,-2) = work3(ix,-1)
          work3(ix,-1) = work3(ix,0)
          work3(ix,0) = work3(ix,1)
          work3(ix,1) = work1(ix,1) + 12.d0 * ux1(ix,iz1)
          work4(ix,-1) = work4(ix,0)
          work4(ix,0) = work4(ix,1)
          work4(ix,1) = work4(ix,2)
          work4(ix,2) = work2(ix,2) + 12.d0 * uz1(ix,iz2)
       enddo
       do ix=1,nx+1
          ix11 = max0( ix-1,1 )
          ix12 = min0( ix+1,nx+1 )
          ix21 = max0( ix-2,1 )
          ix22 = min0( ix+2,nx+1 )
          work5(ix) =   (           ( work3(ix11,-1)-work3(ix,-1) ) &
               + 10.d0 * ( work3(ix11, 0)-work3(ix, 0) ) &
               +         ( work3(ix11, 1)-work3(ix, 1) ) )
          work6(ix,0) = work6(ix,1)
          work6(ix,1) =  (  ( work3(ix11,0)-work3(ix11,1) ) &
               + 10.d0 * ( work3(  ix,0)-work3(ix,  1) ) &
               +         ( work3(ix12,0)-work3(ix12,1) ) )
          
          work7(ix) =( ( work4(ix11,-1)-work4(ix,-1) ) &
               + 10.d0 * ( work4(ix11, 0)-work4(ix, 0) ) &
               +         ( work4(ix11, 1)-work4(ix, 1) ))
          
          work8(ix,0) = work8(ix,1)
          work8(ix,1) =  (           ( work4(ix11,0)-work4(ix11,1) ) &
               + 10.d0 * ( work4(  ix,0)-work4(  ix,1) ) &
               +         ( work4(ix12,0)-work4(ix12,1) ))
          
          work9(ix) = (          work3(ix,-2) - 9.d0 * work3(ix,-1) &
               + 3.d0 * work3(ix,0)  + 5.d0 * work3(ix,1))
          
          work10(ix,-2) = work10(ix,-1)
          work10(ix,-1) = work10(ix,0)
          work10(ix,0) = work10(ix,1)
          work10(ix,1) = ( work3(ix21,1) - 9.d0 * work3(ix11,1) &
               + 3.d0 * work3(  ix,1) + 5.d0 * work3(ix12,1) )
          
          work11(ix) = ( - 5.d0 * work4(ix,-1)  - 3.d0 * work4(ix,0) &
               + 9.d0 * work4(ix, 1)  -        work4(ix,2) )
          
          work12(ix,-1) = work12(ix,0)
          work12(ix,0) = work12(ix,1)
          work12(ix,1) = work12(ix,2)
          work12(ix,2) = ( - 5.d0 * work4(ix11,2) - 3.d0 * work4(  ix,2) &
               + 9.d0 * work4(ix12,2) -        work4(ix22,2))
          
       enddo
       
       do ix=2,nx
          ix21 = max0( ix-2,1 )
          ix22 = min0( ix+2,nx+1 )
          ux(ix,iz) = ux(ix,iz) &
               + ( &
               - (           (   work1(ix-1,-1) + work1(ix-1,1) &
               + work1(ix+1,-1) + work1(ix+1,1) ) &
               + 10.d0 * (   work1(ix-1, 0) + work1(  ix,-1) &
               + work1(  ix, 1) + work1(ix+1, 0) ) &
               + 100.d0 * work1(ix,0) ) &
               + e1(ix,iz) * work5(ix) &
               - e2(ix,iz) * work5(ix+1) &
               + e3(ix,iz) * work6(ix,0) &
               - e4(ix,iz) * work6(ix,1) &
               ) / 144.d0 &
               + e13(ix,iz) * work11(ix-1) &
               + e14(ix,iz) * work11(ix) &
               + e15(ix,iz) * work11(ix+1) &
               + e16(ix,iz) * work11(ix22) &
               + e17(ix,iz) * work12(ix,-1) &
               + e18(ix,iz) * work12(ix,0) &
               + e19(ix,iz) * work12(ix,1) &
               + e20(ix,iz) * work12(ix,2)
          uz(ix,iz) = uz(ix,iz) &
               + ( &
               - (           (   work2(ix-1,-1) + work2(ix-1,1) &
               + work2(ix+1,-1) + work2(ix+1,1) ) &
               + 10.d0 * (   work2(ix-1, 0) + work2(  ix,-1) &
               + work2(  ix, 1) + work2(ix+1, 0) ) &
               + 100.d0 * work2(ix,0) ) &
               + f1(ix,iz) * work7(ix) &
               - f2(ix,iz) * work7(ix+1) &
               + f3(ix,iz) * work8(ix,0) &
               - f4(ix,iz) * work8(ix,1) &
               ) / 144.d0 &
               + f13(ix,iz) * work9(ix21) &
               + f14(ix,iz) * work9(ix-1) &
               + f15(ix,iz) * work9(ix) &
               + f16(ix,iz) * work9(ix+1) &
               + f17(ix,iz) * work10(ix,-2) &
               + f18(ix,iz) * work10(ix,-1) &
               + f19(ix,iz) * work10(ix,0) &
               + f20(ix,iz) * work10(ix,1)
       enddo
    enddo
 endif
 !ux(isx,isz) = ux(isx,isz) + fx(isx,isz)
 !uz(isx,isz) = uz(isx,isz) + fz(isx,isz)
 ! swapping u1 & u2
 do iz=2,nz
    do ix=2,nx
       ux2(ix,iz) = ux1(ix,iz)
       ux1(ix,iz) =  ux(ix,iz)
       uz2(ix,iz) = uz1(ix,iz)
       uz1(ix,iz) =  uz(ix,iz)
    enddo
 enddo
 
 return
end subroutine calstep
	




subroutine create_color_image(image_data_2D,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
     NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,field_number)
	
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
  
  integer NX,NY,it,field_number,ISOURCE,JSOURCE,NPOINTS_PML,nrec
  logical USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX
  
  double precision, dimension(NX,NY) :: image_data_2D
  
  integer, dimension(nrec) :: ix_rec,iy_rec
  
  integer :: ix,iy,irec
  
  character(len=100) :: file_name,system_command1,system_command2,system_command3
	
  integer :: R, G, B
  
  double precision :: normalized_value,max_amplitude
  
  !       open image file and create system command to convert image to more convenient format
  !       use the "convert" command from ImageMagick http://www.imagemagick.org
  if(field_number == 1) then
     write(file_name,"('image',i6.6,'_Ux.pnm')") it
     write(system_command1, "('convert image',i6.6,'_Ux.pnm snapshots/image',i6.6,'_Ux.gif')") it,it
     write(system_command2, "('rm image',i6.6,'_Ux.pnm')") it
  else if(field_number == 2) then
     write(file_name,"('image',i6.6,'_Uz.pnm')") it
     write(system_command1,"('convert image',i6.6,'_Uz.pnm snapshots/image',i6.6,'_Uz.gif')") it,it
     write(system_command2,"('rm image',i6.6,'_Uz.pnm')") it
  endif
  
  open(unit=27, file=file_name, status='unknown')
  
  write(27,"('P3')")	! write image in PNM P3 format
  
  write(27,*) NX,NY	! write image size
  write(27,*) '255'	! maximum value of each pixel color
	
  !       compute maximum amplitude
  max_amplitude = maxval(abs(image_data_2D))
  
  !       image starts in upper-left corner in PNM format
  do iy=NY,1,-1
     do ix=1,NX
        
        !       define data as vector component normalized to [-1:1] and rounded to nearest integer
        !       keeping in mind that amplitude can be negative
        normalized_value = image_data_2D(ix,iy) / max_amplitude
        
        !       suppress values that are outside [-1:+1] to avoid small edge effects
        if(normalized_value < -1.d0) normalized_value = -1.d0
        if(normalized_value > 1.d0) normalized_value = 1.d0
        
        !       draw an orange cross to represent the source
        if((ix >= ISOURCE - width_cross .and. ix <= ISOURCE + width_cross .and. &
             iy >= JSOURCE - thickness_cross .and. iy <= JSOURCE + thickness_cross) .or. &
             (ix >= ISOURCE - thickness_cross .and. ix <= ISOURCE + thickness_cross .and. &
             iy >= JSOURCE - width_cross .and. iy <= JSOURCE + width_cross)) then
           R = 255
           G = 157
           B = 0
           
           !       display two-pixel-thick black frame around the image
        else if(ix <= 2 .or. ix >= NX-1 .or. iy <= 2 .or. iy >= NY-1) then
           R = 0
           G = 0
           B = 0
           
           !       display edges of the PML layers
        else if((USE_PML_XMIN .and. ix == NPOINTS_PML) .or. &
             (USE_PML_XMAX .and. ix == NX - NPOINTS_PML) .or. &
             (USE_PML_YMIN .and. iy == NPOINTS_PML) .or. &
             (USE_PML_YMAX .and. iy == NY - NPOINTS_PML)) then
           R = 255
           G = 150
           B = 0
           
           !       suppress all the values that are below the threshold
        else if(abs(image_data_2D(ix,iy)) <= max_amplitude * cutvect) then
           
           !       use a black or white background for points that are below the threshold
           if(WHITE_BACKGROUND) then
              R = 255
              G = 255
              B = 255
           else
              R = 0
              G = 0
              B = 0
           endif
           
           !       represent regular image points using red if value is positive, blue if negative
        else if(normalized_value >= 0.d0) then
           R = nint(255.d0*normalized_value**POWER_DISPLAY)
           G = 0
           B = 0
        else
           R = 0
           G = 0
           B = nint(255.d0*abs(normalized_value)**POWER_DISPLAY)
        endif
        
        !       draw a green square to represent the receivers
        do irec = 1,nrec
           if((ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
                iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square) .or. &
                (ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
                iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square)) then
              !       use dark green color
              R = 30
              G = 180
              B = 60
           endif
        enddo
        
        !       write color pixel
        write(27,"(i3,' ',i3,' ',i3)") R,G,B
        
     enddo
  enddo
  
  !       close file
close(27)

!       call the system to convert image to GIF (can be commented out if "call system" is missing in your compiler)
call system(system_command1)
call system(system_command2)
end subroutine create_color_image

