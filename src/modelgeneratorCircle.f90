program generator
 implicit none

 integer, parameter     :: prec = kind (1.0) 
 integer                :: i,j
 integer                :: level,tmp

 integer, parameter :: times = 1

 integer, parameter                             :: nbre_layers = 2
 integer, dimension(nbre_layers)                :: layer_thikness
 real(kind=prec), dimension(nbre_layers)        :: vp
 real(kind=prec), dimension(nbre_layers)        :: vs
 real(kind=prec), allocatable                   :: fullvp (:,:)
 real(kind=prec), allocatable                   :: fullvs (:,:)
 real(kind=prec), allocatable                   :: fullrho(:,:),tmpM(:,:)
 real(kind=prec)                                :: dval,pval,xover1,xover2,grad,const
 real(kind=prec), parameter                     :: vpwater = 1.5e3
 integer                                        :: NZ, NX
 integer ::ix,iz,tmpint
 double precision :: boundary,dx,tmpfloat
 integer                                        :: NX_TOTAL, NZ_TOTAL
 integer                                        :: recl_size
!************************************************************************	

! PROGRAM GENERATES VELOCITY MODEL (VP & VS)

!************************************************************************
! Model initialization	
 NX_TOTAL = (800-1)*times+1
 NZ_TOTAL = (800-1)*times+1

!*** Thikness, Vp, Vs *******************************************		

vp(1)=1500.d0;vs(1)=1000.d0;
vp(2)=340.d0;vs(2)=0.d0;
 !layer_thikness(2)=91;vp(2)=3000;vs(2)=1700;
! layer_thikness(3)=162;vp(3)=2000;vs(3)=1170;
! layer_thikness(4)=180;vp(4)=2700;vs(4)=1700;
! layer_thikness(5)=200;vp(5)=3000;vs(5)=1900;

! layer_thikness(1)=40;vp(1)=2200;vs(1)=1400;
! layer_thikness(2)=90;vp(2)=2200;vs(2)=1400;
! layer_thikness(3)=140;vp(3)=2200;vs(3)=1400;
! layer_thikness(4)=180;vp(4)=2200;vs(4)=1400;
! layer_thikness(5)=200;vp(5)=2200;vs(5)=1400;
 
 
 
 
 



 allocate (fullvp (NX_TOTAL, NZ_TOTAL) )
 allocate (fullvs (NX_TOTAL, NZ_TOTAL) )
 allocate (fullrho(NX_TOTAL, NZ_TOTAL) )





 !allocate(tmpM(NZ_TOTAL,NX_TOTAL))

 recl_size = prec * NX_TOTAL * NZ_TOTAL

!****************************************************************	   
 
 open (1,file='../models/2d_circle.vp',form='unformatted',access='direct',recl=recl_size)
 open (2,file='../models/2d_circle.vs',form='unformatted',access='direct',recl=recl_size) 
 open (3,file='../models/2d_circle.rho',form='unformatted',access='direct',recl=recl_size)


 fullvp(:,:)=vp(2)
 fullvs(:,:)=vs(2)
 fullrho(:,:)=1.2041

 if(0.eq.0) then
    tmp=1
    do j=1, NZ_TOTAL
       do i=1, NX_TOTAL
          !if((i-400)**2+(j-400)**2<=300**2) then
             fullvp(i,j) = vp(1)
             fullvs(i,j) = vs(1)
             fullrho(i,j) = 2.d3
          !endif
       enddo
    enddo
 endif
 write(1,rec=1) fullvp(:,:)
 write(2,rec=1) fullvs(:,:)

 if(0.eq.1) then
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
        fullrho(i,j) = dval
     enddo
  enddo

  endif


  write(3,rec=1) fullrho(:,:)

  

 close (1,status='keep')
 close (2,status='keep')
 close (3,status='keep')
end program generator
