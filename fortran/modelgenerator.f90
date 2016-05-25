program generator
 implicit none

 integer, parameter     :: prec = kind (1.0) 
 integer                :: i,j
 integer                :: level,tmp


 integer, parameter                             :: nbre_layers = 5
 integer, dimension(nbre_layers)                :: layer_thikness
 real(kind=prec), dimension(nbre_layers)        :: vp
 real(kind=prec), dimension(nbre_layers)        :: vs
 real(kind=prec), allocatable                   :: fullvp (:,:)
 real(kind=prec), allocatable                   :: fullvs (:,:)
 real(kind=prec), allocatable                   :: fullrho(:,:)
 real(kind=prec)                                :: dval,pval,xover1,xover2,grad,const
 real(kind=prec), parameter                     :: vpwater = 1.5e3
 integer                                        :: NZ, NX
 integer                                        :: NX_TOTAL, NZ_TOTAL
 integer                                        :: recl_size
!************************************************************************	

! PROGRAM GENERATES VELOCITY MODEL (VP & VS)

!************************************************************************
! Model initialization	
 NX_TOTAL = 400
 NZ_TOTAL = 200

!*** Thikness, Vp, Vs *******************************************					
! layer_thikness(1)=20;vp(1)=2100;vs(1)=1500;
! layer_thikness(2)=40;vp(2)=2200;vs(2)=1550;
! layer_thikness(3)=90;vp(3)=2300;vs(3)=1600;
! layer_thikness(4)=140;vp(4)=2500;vs(4)=1650;
! layer_thikness(5)=160;vp(5)=2700;vs(5)=1700;
! layer_thikness(6)=180;vp(6)=2800;vs(6)=1750;
! layer_thikness(7)=200;vp(7)=2900;vs(7)=1800;

 layer_thikness(1)=40;vp(1)=2200;vs(1)=1400;
 layer_thikness(2)=90;vp(2)=2300;vs(2)=1450;
 layer_thikness(3)=140;vp(3)=2500;vs(3)=1550;
 layer_thikness(4)=180;vp(4)=2700;vs(4)=1700;
 layer_thikness(5)=200;vp(5)=3000;vs(5)=1900;


! layer_thikness(1)=20;vp(1)=1500;vs(1)=1300;
! layer_thikness(2)=40;vp(2)=2100;vs(2)=1400;
! layer_thikness(3)=90;vp(3)=2300;vs(3)=1450;
! layer_thikness(4)=140;vp(4)=2300;vs(4)=1450;
! layer_thikness(5)=180;vp(5)=2300;vs(5)=1450;
! layer_thikness(6)=201;vp(6)=2300;vs(6)=1450;

! layer_thikness(1)=40;  vp(1)=2700; vs(1)=1700;
! layer_thikness(2)=90;  vp(2)=2700; vs(2)=1700;
! layer_thikness(3)=140;  vp(3)=2700; vs(3)=1700;
! layer_thikness(4)=180; vp(4)=2700; vs(4)=1700;
! layer_thikness(5)=200; vp(5)=2700; vs(5)=1700;

! layer_thikness(1)=20;vp(1)=1500;vs(1)=1300;
! layer_thikness(2)=25;vp(2)=1550;vs(2)=1320;
! layer_thikness(3)=30;vp(3)=1600;vs(3)=1350;
! layer_thikness(4)=35;vp(4)=1700;vs(4)=1370;
! layer_thikness(5)=40;vp(5)=2100;vs(5)=1400;
! layer_thikness(6)=90;vp(6)=2300;vs(6)=1450;
! layer_thikness(7)=140;vp(7)=2500;vs(7)=1550;
! layer_thikness(8)=150;vp(8)=2600;vs(8)=1580;
! layer_thikness(9)=160;vp(9)=2650;vs(9)=1650;
! layer_thikness(10)=170;vp(10)=2670;vs(10)=1670;
! layer_thikness(11)=180;vp(11)=2700;vs(11)=1700;
! layer_thikness(12)=201;vp(12)=3000;vs(12)=1900;



 



 allocate (fullvp (NX_TOTAL, NZ_TOTAL) )
 allocate (fullvs (NX_TOTAL, NZ_TOTAL) )
 allocate (fullrho(NX_TOTAL, NZ_TOTAL) )
 recl_size = prec * NX_TOTAL * NZ_TOTAL

!****************************************************************	   
 
 open (1,file='./2d_start.vp',form='unformatted',access='direct',recl=recl_size)
 open (2,file='./2d_start.vs',form='unformatted',access='direct',recl=recl_size) 
 open (3,file='./2d_start.rho',form='unformatted',access='direct',recl=recl_size)

 tmp=1
 do j=1, NZ_TOTAL
   level=tmp
   if (j .gt. layer_thikness(level)) level=level+1
   tmp=level
   do i=1, NX_TOTAL
     fullvp(i,j) = vp(level)
     fullvs(i,j) = vs(level)
   enddo
 enddo

 write(1,rec=1) fullvp(:,:)
 write(2,rec=1) fullvs(:,:)


 do j=1,NZ_TOTAL
     do i=1,NX_TOTAL
        pval = fullvp(i,j)
        if (pval.le.vpwater) then
           dval = 1000.E0
        elseif (pval > vpwater .and. pval < 2000.E0) then
           dval = 2351.E0-(7497.E0)*(pval/1000.E0)**(-4.656E0)
        elseif (pval >= 2000.E0 .and. pval <= 2150.E0) then
           xover1 = 2351.E0-(7497.E0)*(2000.E0/1000.E0)**(-4.656E0)
           xover2 = 1740.E0*(2150.E0/1000.E0)**(0.25E0)
           grad = 150.E0/(xover2-xover1)
           const = 2000.E0-(xover1*grad)
           dval = (pval-const)/grad
        elseif (pval > 2150) then
           dval = 1740.E0*(pval/1000.E0)**(0.25E0)
        endif
        fullro(i,j) = dval
     enddo
  enddo
  write(3,rec=1) fullro(:,:)

 close (1,status='keep')
 close (2,status='keep')
 close (3,status='keep')
end program generator
