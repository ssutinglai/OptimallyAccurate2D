
subroutine calculate_circle_boundary(nx,nz,centrenx,centrenz,nradius,LBx,RBx,TBz,BBz,minIX,maxIX)
  implicit none
  integer nx,nz
  integer ix,iz
  integer centrenx,centrenz,nradius  !%! Added for the circle
  integer TBz(nx+1),BBz(nx+1),LBx(nz+1),RBx(nz+1) !%! Added for boundary of circle
  !integer TBz_z(nz+1),BBz_z(nz+1),LBx_x(nx+1),RBx_x(nx+1)
  integer identifier4minIX
  integer minIX, maxIX
  integer k

   do iz=1,nz+1
      LBx(iz)=1
      RBx(iz)=-1
      !BBz_z(iz)=1
      !TBz_z(iz)=-1

   enddo

   do ix=1,nx+1
      BBz(ix)=1
      TBz(ix)=-1
      !LBx_x(ix)=1
      !RBx_x(ix)=-1
   enddo

   do ix=1,nx+1
     if(((ABS((centrenx-ix)).le. nradius).or.(ABS(centrenx-ix).eq. nradius)))then

         TBz(ix)=(centrenz+int(sqrt(dble(nradius**2-(centrenx-ix)**2))))+1
         BBz(ix)=(centrenz-int(sqrt(dble(nradius**2-(centrenx-ix)**2))))-1

     endif
   enddo

   do iz=1,nz+1
       if(((ABS((centrenz-iz)).le. nradius).or.(ABS(centrenz-iz).eq. nradius))) then

         LBx(iz)=(centrenx-int(sqrt(dble(nradius**2-(centrenz-iz)**2))))-1
         RBx(iz)=(centrenx+int(sqrt(dble(nradius**2-(centrenz-iz)**2))))+1
       endif
   enddo

    ! Here we want to define the meaningful ix range

   identifier4minIX=1


    do ix=1,nx+1
       if((identifier4minIX*BBz(ix).ne.1).and.(identifier4minIX.eq.1)) then
          minIX=ix
          identifier4minIX=0
       endif
       if((BBz(ix).ne.1).and.(identifier4minIX.eq.0)) then
          maxIX=ix
          identifier4minIX=1
       endif
    enddo

    write(*,*)minIX,maxIX
    !stop




    !do iz=1,nz+1
     !   LBx_x(LBx(iz))=LBx(iz)






  end subroutine calculate_circle_boundary





