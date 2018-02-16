
subroutine calculate_circle_boundary(nx,nz,centrenx,centrenz,nradius,LBx,RBx,TBz,BBz)
  implicit none
  integer nx,nz
  integer ix,iz
  integer centrenx,centrenz,nradius  !%! Added for the circle
  integer TBz(nx+1),BBz(nx+1),LBx(nz+1),RBx(nz+1) !%! Added for boundary of circle



   do iz=1,nz+1
      LBx(iz)=1
      RBx(iz)=-1
   enddo

   do ix=1,nx+1
      BBz(ix)=1
      TBz(ix)=-2
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

  end subroutine calculate_circle_boundary





