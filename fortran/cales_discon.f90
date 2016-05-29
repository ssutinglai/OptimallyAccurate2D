
subroutine cales_discon( maxnz,nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, & ! hereafter are new variables for cales_discon
     markers,nDiscon,lengthDiscon,dscr)

  ! in the vicinity of boundaries :
  ! we use Operators for discontinuities (Mizutani 2002; Ovcharenko 2015)
  
  integer maxnz,nx,nz
  double precision rho(maxnz+1,*),lam(maxnz+1,*),mu(maxnz+1,*)
  double precision dt,dx,dz
  double precision  e1(maxnz+1,*), e2(maxnz+1,*), e3(maxnz+1,*)
  double precision  e4(maxnz+1,*), e5(maxnz+1,*), e6(maxnz+1,*)
  double precision  e7(maxnz+1,*), e8(maxnz+1,*)
  double precision e13(maxnz+1,*),e14(maxnz+1,*),e15(maxnz+1,*)
  double precision e16(maxnz+1,*),e17(maxnz+1,*),e18(maxnz+1,*)
  double precision e19(maxnz+1,*),e20(maxnz+1,*)
  double precision  f1(maxnz+1,*), f2(maxnz+1,*), f3(maxnz+1,*)
  double precision  f4(maxnz+1,*), f5(maxnz+1,*), f6(maxnz+1,*)
  double precision  f7(maxnz+1,*), f8(maxnz+1,*)
  double precision f13(maxnz+1,*),f14(maxnz+1,*),f15(maxnz+1,*)
  double precision f16(maxnz+1,*),f17(maxnz+1,*),f18(maxnz+1,*)
  double precision f19(maxnz+1,*),f20(maxnz+1,*)
  integer ix,iz
  double precision dt2,dx2,dz2,dxdz
  
  ! interface markers

  integer :: ik,jk,ctr
  integer :: markers(maxnz+1,maxnz+1)
  double precision :: pt0x,pt0z,pt1x,pt1z

  integer :: nDiscon,lengthDiscon ! number of discontinuities
  double precision :: dscr(1:2,1:lengthDiscon,1:nDiscon)

  double precision :: dDiagonal,dDiagonal2 ! sqrt(dx^2+dz^2)
  double precision :: eps ! zero tolerance

  double precision :: xi,zi,distan2 ! intersecion coordinates
  integer :: iLengthDiscon,iDiscon,iInterSection(2),err


  ! Verify that all the coordinates are already with lmargins


  dt2 = dt * dt
  dx2 = dx * dx
  dz2 = dz * dz
  dxdz = dx * dz
  
  dDiagonal2= dx2+dz2
  dDiagonal = sqrt(dDiagonal2)
  
  eps=1.d-8

  ! modified operators for discontinuities
  
  do iz=2,nz
     do ix=2,nx
        if(markers(ix,iz)>0) then

           pt0x=dble(ix-1)*dx
           pt0z=dble(iz-1)*dz


           ! ctr = 1 right-top ix+1,iz+1
           ctr = 1
           distan2 = dDiagonal2
           
           pt1x = pt0x + dx
           pt1z = pt0z + dz

           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,eta,lengthDiscon,nDiscon,iInterSection,err,dscr)

           
           

           ! ctr = 2 right-centre ix+1,iz
           ctr = 2
           distan2 = dx2

           pt1x = pt0x + dx
           pt1z = pt0z

           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,eta,lengthDiscon,nDiscon,iInterSection,err,dscr)

           
           ! ctr = 3 right-bottom ix+1,iz-1
           ctr = 3
           distan2 = dDiagonal2
           
           pt1x = pt0x + dx
           pt1z = pt0z - dz
           
           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,eta,lengthDiscon,nDiscon,iInterSection,err,dscr)
           
           
           


           
           
           ! ctr = 4 centre-top ix,iz+1
           ctr = 4
           distan2 = dz2

           pt1x = pt0x 
           pt1z = pt0z + dz

           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,eta,lengthDiscon,nDiscon,iInterSection,err,dscr)
           
           
           
           ! ctr = 5 centre ix,iz
           ctr = 5
           distan2 = 0.d0

           pt1x = pt0x 
           pt1z = pt0z 
           
           
           ! ctr = 6 centre-bottom ix,iz-1
           ctr = 6
           distan2 = dz2
           
           pt1x = pt0x 
           pt1z = pt0z - dz

           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,eta,lengthDiscon,nDiscon,iInterSection,err,dscr)







           
           
           ! ctr = 7 left-top ix-1,iz+1
           ctr = 7
           distan2 = dDiagonal2
            
           pt1x = pt0x - dx
           pt1z = pt0z + dz

           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,eta,lengthDiscon,nDiscon,iInterSection,err,dscr)

            
           ! ctr = 8 left-centre ix-1,iz
           ctr = 8
           distan2 = dx2
           
           pt1x = pt0x - dx
           pt1z = pt0z 
           
           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,eta,lengthDiscon,nDiscon,iInterSection,err,dscr)

           
            
           ! ctr = 9 left-bottom ix-1,iz-1
           ctr = 9
           distan2 = dDiagonal2

           pt1x = pt0x - dx
           pt1z = pt0z - dz
           
           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,eta,lengthDiscon,nDiscon,iInterSection,err,dscr)




           






                 
                 
                 

                 

























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




        endif
     enddo
  enddo











  return
end subroutine cales_discon






subroutine findNearestPoint(x1,z1,x2,z2,distan2,xi,zi,eta,lengthDiscon,nDiscon,iInterSection,err,dscr)
  implicit none
  double precision :: x1,z1,x2,z2,distan2
  double precision :: xi, zi,eta(0:1,1:2)
  
  integer :: nDiscon,lengthDiscon ! number of discontinuities
  double precision :: dscr(1:2,1:lengthDiscon,1:nDiscon)
  double precision :: distance2(1:lengthDiscon,1:nDiscon)
  integer :: iLengthDiscon,iDiscon,iInterSection(1:2),err
  double precision :: x,z
  
  err=0
  xi=0.d0
  zi=0.d0
  distance2=0.d0

  do iDiscon = 1,nDiscon
     do iLengthDiscon = 1,lengthDiscon
        
        x=dscr(1,iLengthDiscon,iDiscon)
        z=dscr(2,iLengthDiscon,iDiscon)
        
        distance2(iLengthDiscon,iDiscon)=(x-x1)*(x-x1)+(z-z1)*(z-z1)+(x-x2)*(x-x2)+(z-z2)*(z-z2)
        
     enddo
  enddo


  iInterSection(2)=minloc(distance2)
  
  if(distance2(iInterSection(1),iInterSection(2))>distan2) then
     err = 1
  else
     xi = dscr(1,iInterSection(1),iInterSection(2))
     zi = dscr(2,iInterSection(1),iInterSection(2))
     eta(
     
  endif
end subroutine findNearestPoint
     
        
        
  

  
