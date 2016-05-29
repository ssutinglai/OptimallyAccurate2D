
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

  integer :: ik,jk
  integer :: markers(maxnz+1,maxnz+1)
  double precision :: pt0x,pt0z,pt1x,pt1z

  integer :: nDiscon,lengthDiscon ! number of discontinuities
  double precision :: dscr(1:2,1:lengthDiscon,1:nDiscon)

  double precision :: dDiagonal,dDiagonal2 ! sqrt(dx^2+dz^2)
  double precision :: eps ! zero tolerance

  double precision :: xi,zi,distan ! intersecion coordinates
  integer :: iLengthDiscon,iDiscon,nInterSection


  ! Verify that all the coordinates are already with lmargins


  dt2 = dt * dt
  dx2 = dx * dx
  dz2 = dz * dz
  dxdz = dx * dz
  
  dDiagonal2= dx2+dz2
  dDiagonal = sqrt(dDiagonal2)
  
  eps=dDiagonal*1.d-5

  ! modified operators for discontinuities
  
  do iz=2,nz
     do ix=2,nx
        if(markers(ix,iz)>0) then

           pt0x=dble(ix-1)*dx
           pt0z=dble(iz-1)*dz

           ! ctr = 0 right-top ix+1,iz+1

           distan = dDiagonal
           
           ! ctr = 1 right-centre ix+1,iz
           
           distan = dx
           
           ! ctr = 2 right-bottom ix+1,iz-1

            distan = dDiagonal








           ! ctr = 3 centre-top ix,iz+1

            distan = dz

           ! ctr = 5 centre-bottom ix,iz-1


           


           ! ctr = 6 left-top ix-1,iz+1
            distan = dDiagonal
            

            
           ! ctr = 7 left-centre ix-1,iz
           

            
           ! ctr = 8 left-bottom ix-1,iz-1
            distan = dDiagonal

           


           
           do ik=1,-1,-1
              do jk=1,-1,-1
                 
                 ! We are not interested in the same point
                 
                 if((ik.eq.0).and.(jk.eq.0)) cycle
                 
                 pt1x=dble(ix+ik-1)*dx
                 pt1z=dble(iz+jk-1)*dz
                 
                 ! Finding xi,zi (intersection)



                 
                 
                 

                 

























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
