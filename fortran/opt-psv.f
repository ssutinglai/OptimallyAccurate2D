	program normal5
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computation of the synthetic seismograms in the time domain
c using the normal operators.
c --- heterogeneous medium
c
c						1997.6  N.Takeuchi
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer maxnz
	parameter ( maxnz = 1000 )
c
c parameters for the gridding
	real*8 dt,dx,dz
c parameters for the wavefield
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
c parameter for the structure
	real*8 rrho(6),llam(6),mmu(6)
	real*8 rho(maxnz+1,maxnz+1)
	real*8 lam(maxnz+1,maxnz+1),mu(maxnz+1,maxnz+1)
c parameter for the source
	real*8 tp,ts
c parameter for the receiver
	integer nrx,nrz
c parameter for the waveform
	real*8 t
c
c reading the parameter files
	call pinput( maxnz,nt,nx,nz,dt,dx,dz,rrho,llam,mmu,
     &	             tp,ts,nrx,nrz )
c Initializing the data
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
c computing the intermediate parameters
	call calstruct( maxnz,rrho,dx,dz,nx,nz,rho )
	call calstruct( maxnz,llam,dx,dz,nx,nz,lam )
	call calstruct( maxnz,mmu,dx,dz,nx,nz,mu )
	call cales( maxnz,nx,nz,rho,lam,mu,dt,dx,dz,
     &	             e1, e2, e3, e4, e5, e6, e7, e8,
     &	            e13,e14,e15,e16,e17,e18,e19,e20,
     &	             f1, f2, f3, f4, f5, f6, f7, f8,
     &	            f13,f14,f15,f16,f17,f18,f19,f20 )
c
	call datainit( maxnz,maxnz,lam )
	call datainit( maxnz,maxnz,mu )
	ist = dnint( 2 * tp / dt )
	isx = nx / 2 + 1
	isz = nz / 2 + 1
c
	t = 0.d0
	write(6,*) real(t),real(ux(nrx,nrz)),real(uz(nrx,nrz))
	do 100 it=0,nt
	  call calf( maxnz,it,t,ist,isx,isz,dt,dx,dz,rho(isx,isz),
     &	             tp,ts,lam,mu )
c evaluating the next step
	  call calstep( maxnz,nx,nz,
     &	                 e1, e2, e3, e4, e5, e6, e7, e8,
     &	                e13,e14,e15,e16,e17,e18,e19,e20,
     &	                 f1, f2, f3, f4, f5, f6, f7, f8,
     &	                f13,f14,f15,f16,f17,f18,f19,f20,
     &	                ux,uz,ux1,ux2,uz1,uz2,isx,isz,lam,mu,
     &	                 work(1,1), work(1,5), work(1,9),work(1,13),
     &	                work(1,17),work(1,18),work(1,20),work(1,21),
     &	                work(1,23),work(1,24),work(1,28),work(1,29) )
c increment of t
	  t = t + dt
	  write(6,*) real(t),real(ux(nrx,nrz)),real(uz(nrx,nrz))
  100	continue
c
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine pinput( maxnz,nt,nx,nz,dt,dx,dz,rho,lam,mu,
     &	                   tp,ts,nrx,nrz )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer maxnz,nt,nx,nz,nrx,nrz
	real*8 dt,dx,dz,rho(*),lam(*),mu(*),tp,ts
	character*80 tmpfile,dummy
c
	data tmpfile / '/tmp/work' /
c
c temporary file open
	open( unit=11, file=tmpfile, status='unknown' )
c writing to the temporary file
  100	continue
	  read(5,110) dummy
  110	  format(a80)
	  if ( dummy(1:1).eq.'c' ) goto 100
	  if ( dummy(1:3).eq.'end' ) goto 120
	  write(11,*) dummy
	  goto 100
  120	continue
c temporary file close
	close(11)
c 
c temporary file open
	open( unit=11, file=tmpfile, status='unknown' )
c reading the parameter
	read(11,*) nt,nx,nz
	if ( nx.gt.maxnz ) pause 'nx is too large (pinput).'
	if ( nz.gt.maxnz ) pause 'nz is too large (pinput).'
	read(11,*) dt,dx,dz
	read(11,*) rho(1),rho(2),rho(3),rho(4),rho(5),rho(6)
	read(11,*) lam(1),lam(2),lam(3),lam(4),lam(5),lam(6)
	read(11,*) mu(1),mu(2),mu(3),mu(4),mu(5),mu(6)
	read(11,*) tp,ts
	read(11,*) nrx,nrz
c temporary file close
	close(11)
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine datainit( nx,nz,ux )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nx,nz
	real*8 ux(nx+1,*)
	integer i,j
c
	do 110 j=1,nz+1
	  do 100 i=1,nx+1
	     ux(i,j) = 0.d0
  100	  continue
  110	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine calstruct( maxnz,rrho,dx,dz,nx,nz,rho )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer maxnz,nx,nz
	real*8 rrho(4),dx,dz,rho(maxnz+1,*)
	integer i,j,k,nox(6),noz(6)
	real*8 x,z,xmax,zmax,trho,coef1,coef2
c
	data nox / 0, 1, 0, 2, 1, 0 /
	data noz / 0, 0, 1, 0, 1, 2 /
c
	xmax = dx * nx
	zmax = dz * nz
	do 120 j=1,nz+1
	  z = dble(j-1) * dz
	  coef2 = z / zmax
	  do 110 i=1,nx+1
	    x = dble(i-1) * dx
	    coef1 = x / xmax
            trho = 0.d0
	    do 100 k=1,6
	      trho = trho + rrho(k) * coef1**nox(k)
     &	                            * coef2**noz(k)
  100	    continue
	    rho(i,j) = trho
  110	  continue
  120	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine cales( maxnz,nx,nz,rho,lam,mu,dt,dx,dz,
     &	                   e1, e2, e3, e4, e5, e6, e7, e8,
     &	                  e13,e14,e15,e16,e17,e18,e19,e20,
     &	                   f1, f2, f3, f4, f5, f6, f7, f8,
     &	                  f13,f14,f15,f16,f17,f18,f19,f20 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c
	dt2 = dt * dt
	dx2 = dx * dx
	dz2 = dz * dz
	dxdz = dx * dz
c
	do 120 iz=2,nz
	  do 110 ix=2,nx
	    e1(ix,iz) = dt2 / rho(ix,iz)
     &	                * ( ( lam(ix-1,iz) + lam(ix,iz) )
     &	                    + 2.d0 * ( mu(ix-1,iz) + mu(ix,iz) ) )
     &	                / ( 2.d0 * dx2 )
	    e2(ix,iz) = dt2 / rho(ix,iz)
     &	                * ( ( lam(ix,iz) + lam(ix+1,iz) )
     &	                    + 2.d0 * ( mu(ix,iz) + mu(ix+1,iz) ) )
     &	                / ( 2.d0 * dx2 )
	    e3(ix,iz) = dt2 / rho(ix,iz)
     &	                * ( mu(ix,iz-1) + mu(ix,iz) )
     &	                / ( 2.d0 * dz2 )
	    e4(ix,iz) = dt2 / rho(ix,iz)
     &	                * ( mu(ix,iz) + mu(ix,iz+1) )
     &	                / ( 2.d0 * dz2 )
	    e5(ix,iz) = dt2 / rho(ix,iz) * lam(ix-1,iz)
     &	                / ( 4.d0 * dxdz )
	    e6(ix,iz) = dt2 / rho(ix,iz) * lam(ix+1,iz)
     &	                / ( 4.d0 * dxdz )
	    e7(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz-1)
     &	                / ( 4.d0 * dxdz )
	    e8(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+1)
     &	                / ( 4.d0 * dxdz )
	    e13(ix,iz) = dt2 / rho(ix,iz) * lam(ix-1,iz)
     &	                * ( -5.d0 ) / ( 1728.d0 * dxdz )
	    e14(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz)
     &	                * ( -3.d0 ) / ( 1728.d0 * dxdz )
	    e15(ix,iz) = dt2 / rho(ix,iz) * lam(ix+1,iz)
     &	                * ( +9.d0 ) / ( 1728.d0 * dxdz )
	    if ( ix+2.le.nx+1 ) then
	      e16(ix,iz) = dt2 / rho(ix,iz) * lam(ix+2,iz)
     &	                * ( -1.d0 ) / ( 1728.d0 * dxdz )
	    else
	      e16(ix,iz) = 0.d0
	    endif
	    e17(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz-1)
     &	                * ( -5.d0 ) / ( 1728.d0 * dxdz )
	    e18(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz)
     &	                * ( -3.d0 ) / ( 1728.d0 * dxdz )
	    e19(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+1)
     &	                * ( +9.d0 ) / ( 1728.d0 * dxdz )
	    if ( iz+2.le.nz+1) then
	      e20(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+2)
     &	                * ( -1.d0 ) / ( 1728.d0 * dxdz )
	    else
	      e20(ix,iz) = 0.d0
	    endif
	    f1(ix,iz) = dt2 / rho(ix,iz)
     &	                * ( mu(ix-1,iz) + mu(ix,iz) )
     &	                / ( 2.d0 * dx2 )
	    f2(ix,iz) = dt2 / rho(ix,iz)
     &	                * ( mu(ix,iz) + mu(ix+1,iz) )
     &	                / ( 2.d0 * dx2 )
	    f3(ix,iz) = dt2 / rho(ix,iz)
     &	                * ( ( lam(ix,iz-1) + lam(ix,iz) )
     &	                    + 2.d0 * ( mu(ix,iz-1) + mu(ix,iz) ) )
     &	                / ( 2.d0 * dz2 )
	    f4(ix,iz) = dt2 / rho(ix,iz)
     &	                * ( ( lam(ix,iz) + lam(ix,iz+1) )
     &	                     + 2.d0 * ( mu(ix,iz) + mu(ix,iz+1) ) )
     &	                / ( 2.d0 * dz2 )
	    f5(ix,iz) = dt2 / rho(ix,iz) * mu(ix-1,iz)
     &	                / ( 4.d0 * dxdz )
	    f6(ix,iz) = dt2 / rho(ix,iz) * mu(ix+1,iz)
     &	                / ( 4.d0 * dxdz )
	    f7(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-1)
     &	                / ( 4.d0 * dxdz )
	    f8(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz+1)
     &	                / ( 4.d0 * dxdz )
	    if ( ix-2.ge.1 ) then
	      f13(ix,iz) = dt2 / rho(ix,iz) * mu(ix-2,iz)
     &	                * (  1.d0 ) / ( 1728.d0 * dxdz )
	    else
	      f13(ix,iz) = 0.d0
	    endif
	    f14(ix,iz) = dt2 / rho(ix,iz) * mu(ix-1,iz)
     &	                * ( -9.d0 ) / ( 1728.d0 * dxdz )
	    f15(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz)
     &	                * (  3.d0 ) / ( 1728.d0 * dxdz )
	    f16(ix,iz) = dt2 / rho(ix,iz) * mu(ix+1,iz)
     &	                * (  5.d0 ) / ( 1728.d0 * dxdz )
	    if ( iz-2.ge.1 ) then
	      f17(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-2)
     &	                * (  1.d0 ) / ( 1728.d0 * dxdz )
	    else
	      f17(ix,iz) = 0.d0
	    endif
	    f18(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-1)
     &	                * ( -9.d0 ) / ( 1728.d0 * dxdz )
	    f19(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz)
     &	                * (  3.d0 ) / ( 1728.d0 * dxdz )
	    f20(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz+1)
     &	                * (  5.d0 ) / ( 1728.d0 * dxdz )
  110	  continue
  120	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine calf( maxnz,it,t,ist,isx,isz,dt,dx,dz,rho,
     &	                 tp,ts,fx,fz )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 pi
	parameter ( pi=3.1415926535897932d0 )
	integer maxnz,it,ist,isx,isz
	real*8 t,dt,dx,dz,rho,tp,ts,fx(maxnz+1,*),fz(maxnz+1,*)
	real*8 b
c
	if ( it.le.ist ) then
	  b = pi * ( t - ts ) / tp
	  fx(isx,isz)
     &	      = dsqrt(pi) / 2.d0 * (b*b-0.5d0) * dexp(-b*b)
     &	        / ( dx * dz )
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
c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine calstep( maxnz,nx,nz,
     &	                     e1, e2, e3, e4, e5, e6, e7, e8,
     &	                    e13,e14,e15,e16,e17,e18,e19,e20,
     &	                     f1, f2, f3, f4, f5, f6, f7, f8,
     &	                    f13,f14,f15,f16,f17,f18,f19,f20,
     &	                    ux,uz,ux1,ux2,uz1,uz2,isx,isz,fx,fz,
     &	                    work1,work2,work3,work4,
     &	                    work5,work6,work7,work8,
     &	                    work9,work10,work11,work12 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c
c predicting the wavefield
	do 110 iz=2,nz
	  do 100 ix=2,nx
	    ux(ix,iz) = 2.d0 * ux1(ix,iz) - ux2(ix,iz)
     &	      + e1(ix,iz) * ( ux1(ix-1,iz) - ux1(ix,iz) )
     &	      + e2(ix,iz) * ( ux1(ix+1,iz) - ux1(ix,iz) )
     &	      + e3(ix,iz) * ( ux1(ix,iz-1) - ux1(ix,iz) )
     &	      + e4(ix,iz) * ( ux1(ix,iz+1) - ux1(ix,iz) )
     &	      - e5(ix,iz) * ( uz1(ix-1,iz+1) - uz1(ix-1,iz-1) )
     &	      + e6(ix,iz) * ( uz1(ix+1,iz+1) - uz1(ix+1,iz-1) )
     &	      - e7(ix,iz) * ( uz1(ix+1,iz-1) - uz1(ix-1,iz-1) )
     &	      + e8(ix,iz) * ( uz1(ix+1,iz+1) - uz1(ix-1,iz+1) )
	    uz(ix,iz) = 2.d0 * uz1(ix,iz) - uz2(ix,iz)
     &	      + f1(ix,iz) * ( uz1(ix-1,iz) - uz1(ix,iz) )
     &	      + f2(ix,iz) * ( uz1(ix+1,iz) - uz1(ix,iz) )
     &	      + f3(ix,iz) * ( uz1(ix,iz-1) - uz1(ix,iz) )
     &	      + f4(ix,iz) * ( uz1(ix,iz+1) - uz1(ix,iz) )
     &	      - f5(ix,iz) * ( ux1(ix-1,iz+1) - ux1(ix-1,iz-1) )
     &	      + f6(ix,iz) * ( ux1(ix+1,iz+1) - ux1(ix+1,iz-1) )
     &	      - f7(ix,iz) * ( ux1(ix+1,iz-1) - ux1(ix-1,iz-1) )
     &	      + f8(ix,iz) * ( ux1(ix+1,iz+1) - ux1(ix-1,iz+1) )
  100	  continue
  110	continue
	ux(isx,isz) = ux(isx,isz) + fx(isx,isz)
	uz(isx,isz) = uz(isx,isz) + fz(isx,isz)
c correcting the wavefield
c
	do 120 ix=2,nx
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
  120	continue
	do 130 ix=1,nx+1
	  ix11 = max0( ix-1,1 )
	  ix12 = min0( ix+1,nx+1 )
	  ix21 = max0( ix-2,1 )
	  ix22 = min0( ix+2,nx+1 )
	  work6(ix,0) = 0.d0
	  work6(ix,1) =
     &	    (           ( -work3(ix11,1) )
     &	      + 10.d0 * ( -work3(ix,  1) )
     &	      +         ( -work3(ix12,1) )
     &	    )
	  work8(ix,0) = 0.d0
	  work8(ix,1) =
     &	    (           ( -work4(ix11,1) )
     &	      + 10.d0 * ( -work4(  ix,1) )
     &	      +         ( -work4(ix12,1) )
     &	    )
	  work10(ix,-2) = 0.d0
	  work10(ix,-1) = 0.d0
	  work10(ix,0) = 0.d0
	  work10(ix,1)
     &	    = (          work3(ix21,1) - 9.d0 * work3(ix11,1)
     &	        + 3.d0 * work3(  ix,1) + 5.d0 * work3(ix12,1)
     &	      )
	  work12(ix,-1) = 0.d0
	  work12(ix,0) = 0.d0
	  work12(ix,1)
     &	    = ( - 5.d0 * work4(ix11,1) - 3.d0 * work4(  ix,1)
     &	        + 9.d0 * work4(ix12,1) -        work4(ix22,1)
     &	      )
	  work12(ix,2)
     &	    = ( - 5.d0 * work4(ix11,2) - 3.d0 * work4(  ix,2)
     &	        + 9.d0 * work4(ix12,2) -        work4(ix22,2)
     &	      )
  130	continue
c
	do 180 iz=2,nz
	  iz1 = iz + 1
	  iz2 = min0( iz+2, nz+1 )
	  do 140 ix=2,nx
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
  140	  continue
	  do 150 ix=1,nx+1
	    ix11 = max0( ix-1,1 )
	    ix12 = min0( ix+1,nx+1 )
	    ix21 = max0( ix-2,1 )
	    ix22 = min0( ix+2,nx+1 )
	    work5(ix) =
     &	      (           ( work3(ix11,-1)-work3(ix,-1) )
     &	        + 10.d0 * ( work3(ix11, 0)-work3(ix, 0) )
     &	        +         ( work3(ix11, 1)-work3(ix, 1) )
     &	      )
	    work6(ix,0) = work6(ix,1)
	    work6(ix,1) =
     &	      (           ( work3(ix11,0)-work3(ix11,1) )
     &	        + 10.d0 * ( work3(  ix,0)-work3(ix,  1) )
     &	        +         ( work3(ix12,0)-work3(ix12,1) )
     &	      )
	    work7(ix) =
     &	      (           ( work4(ix11,-1)-work4(ix,-1) )
     &	        + 10.d0 * ( work4(ix11, 0)-work4(ix, 0) )
     &	        +         ( work4(ix11, 1)-work4(ix, 1) )
     &	      )
	    work8(ix,0) = work8(ix,1)
	    work8(ix,1) =
     &	      (           ( work4(ix11,0)-work4(ix11,1) )
     &	        + 10.d0 * ( work4(  ix,0)-work4(  ix,1) )
     &	        +         ( work4(ix12,0)-work4(ix12,1) )
     &	      )
	    work9(ix)
     &	      = (          work3(ix,-2) - 9.d0 * work3(ix,-1)
     &	          + 3.d0 * work3(ix,0)  + 5.d0 * work3(ix,1)
     &	        )
	    work10(ix,-2) = work10(ix,-1)
	    work10(ix,-1) = work10(ix,0)
	    work10(ix,0) = work10(ix,1)
	    work10(ix,1)
     &	      = (          work3(ix21,1) - 9.d0 * work3(ix11,1)
     &	          + 3.d0 * work3(  ix,1) + 5.d0 * work3(ix12,1)
     &	        )
	    work11(ix)
     &	      = ( - 5.d0 * work4(ix,-1)  - 3.d0 * work4(ix,0)
     &	          + 9.d0 * work4(ix, 1)  -        work4(ix,2)
     &	        )
	    work12(ix,-1) = work12(ix,0)
	    work12(ix,0) = work12(ix,1)
	    work12(ix,1) = work12(ix,2)
	    work12(ix,2)
     &	      = ( - 5.d0 * work4(ix11,2) - 3.d0 * work4(  ix,2)
     &	          + 9.d0 * work4(ix12,2) -        work4(ix22,2)
     &	        )
  150	  continue
	  do 170 ix=2,nx
	    ix21 = max0( ix-2,1 )
	    ix22 = min0( ix+2,nx+1 )
	    ux(ix,iz) = ux(ix,iz)
     &	      + (
     &	        - (           (   work1(ix-1,-1) + work1(ix-1,1)
     &	                        + work1(ix+1,-1) + work1(ix+1,1) )
     &	            + 10.d0 * (   work1(ix-1, 0) + work1(  ix,-1)
     &	                        + work1(  ix, 1) + work1(ix+1, 0) )
     &	            + 100.d0 * work1(ix,0) )
     &	        + e1(ix,iz) * work5(ix)
     &	        - e2(ix,iz) * work5(ix+1)
     &	        + e3(ix,iz) * work6(ix,0)
     &	        - e4(ix,iz) * work6(ix,1)
     &	         ) / 144.d0
     &	      + e13(ix,iz) * work11(ix-1)
     &	      + e14(ix,iz) * work11(ix)
     &	      + e15(ix,iz) * work11(ix+1)
     &	      + e16(ix,iz) * work11(ix22)
     &	      + e17(ix,iz) * work12(ix,-1)
     &	      + e18(ix,iz) * work12(ix,0)
     &	      + e19(ix,iz) * work12(ix,1)
     &	      + e20(ix,iz) * work12(ix,2)
	    uz(ix,iz) = uz(ix,iz)
     &	      + (
     &	        - (           (   work2(ix-1,-1) + work2(ix-1,1)
     &	                        + work2(ix+1,-1) + work2(ix+1,1) )
     &	            + 10.d0 * (   work2(ix-1, 0) + work2(  ix,-1)
     &	                        + work2(  ix, 1) + work2(ix+1, 0) )
     &	            + 100.d0 * work2(ix,0) )
     &	        + f1(ix,iz) * work7(ix)
     &	        - f2(ix,iz) * work7(ix+1)
     &	        + f3(ix,iz) * work8(ix,0)
     &	        - f4(ix,iz) * work8(ix,1)
     &	         ) / 144.d0
     &	      + f13(ix,iz) * work9(ix21)
     &	      + f14(ix,iz) * work9(ix-1)
     &	      + f15(ix,iz) * work9(ix)
     &	      + f16(ix,iz) * work9(ix+1)
     &	      + f17(ix,iz) * work10(ix,-2)
     &	      + f18(ix,iz) * work10(ix,-1)
     &	      + f19(ix,iz) * work10(ix,0)
     &	      + f20(ix,iz) * work10(ix,1)
  170	  continue
  180	continue
	ux(isx,isz) = ux(isx,isz) + fx(isx,isz)
	uz(isx,isz) = uz(isx,isz) + fz(isx,isz)
c swapping u1 & u2
	do 230 iz=2,nz
	  do 220 ix=2,nx
	    ux2(ix,iz) = ux1(ix,iz)
	    ux1(ix,iz) =  ux(ix,iz)
	    uz2(ix,iz) = uz1(ix,iz)
	    uz1(ix,iz) =  uz(ix,iz)
  220	  continue
  230	continue
c
	return
	end
