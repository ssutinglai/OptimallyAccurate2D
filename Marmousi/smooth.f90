program smooth

 implicit none 
!local variables:
 integer, parameter             :: prec = kind (1.0) 
 integer                        :: i,j,i_tmp,j_tmp
 integer                        :: topend,botend,lefend,rigend,winlen
 integer                        :: NX_TOTAL, NZ_TOTAL
 real(kind=prec), allocatable   :: fullvp(:,:)
 real(kind=prec), allocatable   :: fullvs(:,:)
 integer                        :: recl_size
 real                           :: pi,sigma
 real                           :: G_sm,G_sm_tmp,fullvp_tmp,fullvs_tmp
!*********************************************************************
 pi = 4.*atan(1.)

! Model initialization	
 NX_TOTAL = 300
 NZ_TOTAL = 100

 recl_size=prec*NX_TOTAL*NZ_TOTAL
 
 topend = 1
 botend = 22
 lefend = 1
 rigend = 300


 winlen = 15 ! half of the window length
 sigma  = 10.E0

! winlen = 2 ! half of the window length
! sigma  = 1.E0

! allocate vp,vs
 allocate (fullvp (NX_TOTAL, NZ_TOTAL) )
 allocate (fullvs (NX_TOTAL, NZ_TOTAL) )

!*******************************************************************
! initiliaze data
 open (1,file='./marmousi_vp',form='unformatted',access='direct',recl=recl_size)
 open (2,file='./marmousi_vs',form='unformatted',access='direct',recl=recl_size)
 open (3,file='./marmousi_vp_sm',access='direct',form='unformatted',recl=recl_size)
 open (4,file='./marmousi_vs_sm',access='direct',form='unformatted',recl=recl_size)


 read(1,rec=1) fullvp(:,:)
 read(2,rec=1) fullvs(:,:)


!*******************************************************************
! executable code:
! Gaussian Smoothing

! if (topend < winlen+1)           topend = winlen + 1
! if (botend > NZ_TOTAL-winlen-1)  botend = NZ_TOTAL - winlen - 1
! if (lefend < winlen+1)           lefend = winlen + 1
! if (rigend > NX_TOTAL-winlen-1)  rigend = NX_TOTAL - winlen - 1

 write (*,*) 'topend=', topend
 write (*,*) 'botend=', botend
 write (*,*) 'lefend=', lefend
 write (*,*) 'rigend=', rigend
 write (*,*) 'winlen=', winlen


!*** ***********************
 
 do j = topend, botend
   do i = lefend,rigend
 
     G_sm =0.0E0 
     fullvp_tmp = 0.0E0
     fullvs_tmp = 0.0E0
     do j_tmp =  -winlen,winlen
       do i_tmp =  -winlen,winlen 
         if (i+i_tmp <1 .or. i+i_tmp > NX_TOTAL .or. &
                 j+j_tmp <1 .or. j+j_tmp > NZ_TOTAL) then
           G_sm_tmp = 0.0E0
         else  
           G_sm_tmp = 0.5E0 * pi / (sigma**2)*exp(-(i_tmp**2+j_tmp**2)/2/sigma**2)
           fullvp_tmp = fullvp_tmp + fullvp(i+i_tmp,j+j_tmp) * G_sm_tmp
           fullvs_tmp = fullvs_tmp + fullvs(i+i_tmp,j+j_tmp) * G_sm_tmp
         endif
         G_sm = G_sm + G_sm_tmp
       enddo
     enddo
     fullvp(i,j) = fullvp_tmp / G_sm
     fullvs(i,j) = fullvs_tmp / G_sm          
   enddo
 enddo


!**********************************************************************
! write the data

 write(3,rec=1) fullvp(:,:)
 write(4,rec=1) fullvs(:,:)

 close (1,status='keep')
 close (2,status='keep')
 close (3,status='keep')
 close (4,status='keep')


!**********************************************************************
 deallocate (fullvp,fullvs)

end program smooth
