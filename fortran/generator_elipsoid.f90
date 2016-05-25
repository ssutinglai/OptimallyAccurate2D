program generator_elipsoid
 implicit none

 integer, parameter                             :: prec = kind (1.0) 
 integer                                        :: i,j

 real, allocatable                              :: fullvp (:,:)
 real, allocatable                              :: fullvs (:,:)
 integer                                        :: NZ, NX
 real                                           :: xc, zc, xr, zr
 real                                           :: dp, ds
 real                                           :: rad
 integer                                        :: NX_TOTAL, NZ_TOTAL
 integer                                        :: recl_size
!************************************************************************	

! PROGRAM GENERATES VELOCITY MODEL (VP & VS)

!************************************************************************
! Model initialization	
 NX_TOTAL = 400
 NZ_TOTAL = 200

!*** Thikness, Vp, Vs *******************************************					
! layer_thikness(1)=20;vp(1)=1900;vs(1)=1100;
! layer_thikness(2)=40;vp(2)=2100;vs(2)=1250;
! layer_thikness(3)=90;vp(3)=2300;vs(3)=1350;
! layer_thikness(4)=140;vp(4)=2500;vs(4)=1470;
! layer_thikness(5)=180;vp(5)=2700;vs(5)=1600;
! layer_thikness(6)=201;vp(6)=3000;vs(6)=1760;

! layer_thikness(1)=20;  vp(1)=4000; vs(1)=2300;
! layer_thikness(2)=40;  vp(2)=4000; vs(2)=2300;
! layer_thikness(3)=90;  vp(3)=4000; vs(3)=2300;
! layer_thikness(4)=140; vp(4)=4000; vs(4)=2300;
! layer_thikness(5)=180; vp(5)=4000; vs(5)=2300;
! layer_thikness(6)=201; vp(6)=4000; vs(6)=2300;

 allocate (fullvp (NX_TOTAL, NZ_TOTAL) )
 allocate (fullvs (NX_TOTAL, NZ_TOTAL) )
 recl_size = prec * NX_TOTAL * NZ_TOTAL

!****************************************************************	    
! Writing velocities to a output-files db_vp.3d and db_vs.3d
 open (1,file='./2d_start.vp',form='unformatted',access='direct',recl=recl_size)
 open (2,file='./2d_start.vs',form='unformatted',access='direct',recl=recl_size)
 open (3,file='./2d_elipsoid.vp',access='direct',form='unformatted',recl=recl_size)
 open (4,file='./2d_elipsoid.vs',access='direct',form='unformatted',recl=recl_size)

!**************************************************************** 
! Perturbation Attributes

! main perturbation
 xc = 200.0E0   ! Position in X direction
 zc = 90.0E0     ! Position in Z direction
 xr = 30.0E0     ! Radius in X direciton
 zr = 8.0E0      ! Radius in Z direction

 dp = 300.0E0    ! Perturbation of Vp
 ds = 150.0E0    ! Perturbation of Vs


!unknown perturbation
! xc = 250.0E0   ! Position in X direction
! zc = 160.0E0     ! Position in Z direction
! xr = 6.0E0     ! Radius in X direciton
! zr = 4.0E0      ! Radius in Z direction

! dp = 30.0E0    ! Perturbation of Vp
! ds = 15.0E0    ! Perturbation of Vs
 
!******************************************************************

 read(1,rec=1) fullvp(:,:)
 read(2,rec=1) fullvs(:,:)


 do j=1, NZ_TOTAL

   do i=1, NX_TOTAL

     rad = (i-xc)*(i-xc)/(xr*xr) + (j-zc)*(j-zc)/(zr*zr)

     if (rad < 1.0 .and. j <= zc) then
       fullvp(i,j) = fullvp(i,j)  - dp
       fullvs(i,j) = fullvs(i,j)  - ds
     elseif (rad <1.0 .and. j > zc) then
       fullvp(i,j) = fullvp(i,j)  + dp
       fullvs(i,j) = fullvs(i,j)  + ds
     endif

   enddo
 enddo
 
 write(3,rec=1) fullvp(:,:)
 write(4,rec=1) fullvs(:,:)

 close (1,status='keep')
 close (2,status='keep')
 close (3,status='keep')
 close (4,status='keep')

 deallocate (fullvp,fullvs)

end program generator_elipsoid
