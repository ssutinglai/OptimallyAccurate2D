program multipleSourcesOPT2D

  ! Computation of the synthetic seismograms in the time domain
  ! using the optimally accurate operators.
  ! 2D PSV heterogeneous medium
  ! CPML or Cerjan boundary conditions
  !
  !					originally from 1997.6  N. Takeuchi
  !                                                     2016.5. N. Fuji
  !                         discon operators (matlab) : 2016.5. O. Ovcharenko
  !                                         colorbars : 2016.5. K. Okubo
  !  
  !                                          cleaning : 2016.6. N. Fuji   

  use parameters
  implicit none

  ! Cerjan boundary
  lmargin(1)=NPOINTS_PML
  rmargin(1)=NPOINTS_PML
  lmargin(2)=NPOINTS_PML
  rmargin(2)=NPOINTS_PML

  call paramMultiReader
  
  call vectorAllocate
  
  call disconConfig ! discontinuity configuration
  
  call ReceiverSourcePositions 

  !!!! for each source we calculate synthetics
  
  ! reading intermediate parameters (vp,vs,rho)
  
  call calstruct( maxnx,maxnz,rhofile,nx,nz,rho )
  call calstruct( maxnx,maxnz,vpfile, nx,nz,vp )
  call calstruct( maxnx,maxnz,vsfile, nx,nz,vs )
    
  call freeConfig

  ! calculate lamda and mu
  call calstruct2(maxnx,maxnz,nx,nz,rho,vp,vs,lam,mu,liquidmarkers)
  
  !write(12,*) rho
  !write(13,*) vp
  !write(14,*) vs


  call calstructBC(maxnx, maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)


  ! Smoothed version of CONV/OPT operators

  call cales( nx,nz,rho,lam,mu,dt,dx,dz, &
       e1, e2, e3, e4, e5, e6, e7, e8, &
       e13,e14,e15,e16,e17,e18,e19,e20, &
       f1, f2, f3, f4, f5, f6, f7, f8, &
       f13,f14,f15,f16,f17,f18,f19,f20 )

  ! discontinuities
  
  ee12 = 0.d0
  ee34 = 0.d0
  ee56 = 0.d0
  ee65 = 0.d0
  ee78 = 0.d0
  ee87 = 0.d0

  ff12 = 0.d0
  ff34 = 0.d0
  ff56 = 0.d0
  ff65 = 0.d0
  ff78 = 0.d0
  ff87 = 0.d0
  
  if(nDiscon.ne.0) then
 
     ! changing dscr by putting lmargin(1) and (2)
     tmpvaluex=dble(lmargin(1))*dx
     tmpvaluez=dble(lmargin(2))*dz
     do ix=1,nDiscon
        do iz=1,lengthDiscon
           dscr(1,iz,ix)=dscr(1,iz,ix)+tmpvaluex
           dscr(2,iz,ix)=dscr(2,iz,ix)+tmpvaluez
        enddo
     enddo
          
     call cales_discon( nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, & 
     ! hereafter are new variables for cales_discon
     ee12,ee34,ee56,ee65,ee78,ee87, &
     ff12,ff34,ff56,ff65,ff78,ff87, &
     markers,nDiscon,lengthDiscon,dscr)

  endif
  
  if(lengthFreeSurface.ne.0) then
          
     call cales_free( maxnx,nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, & 
     ! hereafter are new variables for cales_discon
     ee12,ee34,ee56,ee65,ee78,ee87, &
     ff12,ff34,ff56,ff65,ff78,ff87, &
     zerodisplacement,lengthFreeSurface,free)

     
  endif



  if(lengthFreeSurface.ne.0) then
     ! changing free by putting lmargin(1) and (2)
     tmpvaluex=dble(lmargin(1))*dx
     tmpvaluez=dble(lmargin(2))*dz
     
     do ix=1,lengthFreeSurface
           free(1,ix)=free(1,ix)+tmpvaluex
           free(2,ix)=free(2,ix)+tmpvaluez           
     enddo
  endif

  ! for Cerjan absorbing boundary


  weightBC=1.d0
     
  call compNRBCpre(weightBC(1:nx+1,1:nz+1),CerjanRate,lmargin,rmargin,nx+1,nz+1)
  

  do ir= 1, nReceiver
     nrx(ir)=nrx(ir)+lmargin(1)
     nrz(ir)=nrz(ir)+lmargin(2)
  enddo
  

  do iSource = 1, nSource
    
     isx=iisx(iSource)+lmargin(1)
     isz=iisz(iSource)+lmargin(2)
     
     ist=nt/4
    

     ! for video (without boundary)
     recl_size=(nx+1)*(nz+1)*kind(0e0)
    
     
     ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0) ! from Festa and Vilotte
     
     ! Initializing the data
     call datainit( maxnz,maxnz,ux )
     call datainit( maxnz,maxnz,uz )
     call datainit( maxnz,maxnz,ux1 )
     call datainit( maxnz,maxnz,uz1 )
     call datainit( maxnz,maxnz,ux2 )
     call datainit( maxnz,maxnz,uz2 )
     

     call datainit( maxnz,31,work )
     
     
 
     ! R. Courant et K. O. Friedrichs et H. Lewy (1928)
     cp=maxval(vp)
     Courant_number = cp * dt * sqrt(1.d0/dx**2 + 1.d0/dz**2)
     print *, 'Courant number is', Courant_number
     
     


     
     call datainit( maxnz,maxnz,fx)
     call datainit( maxnz,maxnz,fz)
     

     ! ist = dnint( 2 * tp / dt )
     ! isx = nx / 2 + 1
     ! isz = nz / 2 + 1
     
   
 
     
     t=0.d0
     time(0)=t
     do it=0,nt
        print *, isx,isz
        call calf2( nx,nz,it,t,ist,isx,isz,dt,dx,dz,rho(isx,isz),f0,t0,fx,fz )
        t=t+dt
        !write(13,*) t, fx(isx,isz),fz(isx,isz)
        
     enddo
     !print *, maxnz,it,t,ist,isx,isz,dt,dx,dz,rho(isx,isz),f0,t0
     !stop

     
     t = 0.d0
     !write(14,*) real(t),real(ux(nrx,nrz)),real(uz(nrx,nrz))
     do ir = 1,nReceiver
        synx(0,ir)=ux(nrx(ir),nrz(ir))
        synz(0,ir)=uz(nrx(ir),nrz(ir))
     enddo

     

     do it=0,nt
        call calf2( maxnz,it,t,ist,isx,isz,dt,dx,dz,rho(isx,isz),f0,t0,fx,fz )
        ! evaluating the next step
        
        !if(nDiscon.eq.0) then
        !   call calstep( maxnz,nx,nz, &
        !        e1, e2, e3, e4, e5, e6, e7, e8, &
        !        e13,e14,e15,e16,e17,e18,e19,e20, &
        !        f1, f2, f3, f4, f5, f6, f7, f8, &
        !        f13,f14,f15,f16,f17,f18,f19,f20, &
        !        ux,uz,ux1,ux2,uz1,uz2,isx,isz,fx,fz, &
        !        work(1,1), work(1,5), work(1,9),work(1,13), &
        !        work(1,17),work(1,18),work(1,20),work(1,21), &
        !        work(1,23),work(1,24),work(1,28),work(1,29), optimise)
           
        !else
           call calstep_discon( nx,nz, &
                e1, e2, e3, e4, e5, e6, e7, e8, &
                e13,e14,e15,e16,e17,e18,e19,e20, &
                f1, f2, f3, f4, f5, f6, f7, f8, &
                f13,f14,f15,f16,f17,f18,f19,f20, &
                ux,uz,ux1,ux2,uz1,uz2,isx,isz,fx,fz, &
                work(1,1), work(1,5), work(1,9),work(1,13), &
                work(1,17),work(1,18),work(1,20),work(1,21), &
                work(1,23),work(1,24),work(1,28),work(1,29), optimise, & 
                ! Hereafter are new variables for cales_discon
                ee12,ee34,ee56,ee65,ee78,ee87, &
                ff12,ff34,ff56,ff65,ff78,ff87)
           

        !endif


           if(lengthFreeSurface.ne.0) then
              do iz=1,nz+1
                 do ix=1,nx+1
                    if(zerodisplacement(ix,iz).eq.1) then
                       ux(ix,iz) = 0.d0
                       uz(ix,iz) = 0.d0
                    endif
                 enddo
              enddo
           endif


           ! increment of t
        t = t + dt
        time(it)=t
        !write(14,*) real(t),real(ux(nrx,nrz)),real(uz(nrx,nrz))
        do ir = 1,nReceiver
           synx(it,ir)=ux(nrx(ir),nrz(ir))
           synz(it,ir)=uz(nrx(ir),nrz(ir))
        enddo
        
     
        
        ! applying Cerjan boundary
        
        do iz=1,nz+1
           do ix=1,nx+1
              uz(ix,iz)=uz(ix,iz)*weightBC(ix,iz)
              ux1(ix,iz)=ux1(ix,iz)*weightBC(ix,iz)
              uz1(ix,iz)=uz1(ix,iz)*weightBC(ix,iz)
              ux2(ix,iz)=ux2(ix,iz)*weightBC(ix,iz)
              uz2(ix,iz)=uz2(ix,iz)*weightBC(ix,iz)
           enddo
        enddo
        
        
        
        ! calculating strains
        
        
        if(writingStrain.and.(mod(it,IT_DISPLAY).eq.0)) then
           singleStrainDiagonal=0.e0
           tmpsingleStrain=0.e0
           call calStrainDiagonal(nx,nz,ux,uz,lmargin,rmargin,singleStrainDiagonal)
          


           if(optimise) then
              write(outfile,'("strain",I5,".",I5,".",I5,".OPT_dat") ') it,isx-lmargin(1),isz-lmargin(2)
           else
              write(outfile,'("strain",I5,".",I5,".",I5,".CON_dat") ') it,isx-lmargin(1),isz-lmargin(2)
           endif
           do j=1,24
              if(outfile(j:j).eq.' ') outfile(j:j)='0'
           enddo
           
           outfile = './strains/'//trim(modelname)//'/'//outfile
           open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)

           tmpsingleStrain(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2)) = &
                singleStrainDiagonal(lmargin(1)+1:nx+1-rmargin(1),lmargin(2)+1:nz+1-rmargin(2))
           write(1,rec=1)  tmpsingleStrain
           close(1,status='keep')
        endif


        
        !write(*,*) it, ' of ', nt
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
           
           
           
           !call create_color_image(ux(1:nx+1,1:nz+1),nx+1,nz+1,it,isx,isz,ix_rec,iz_rec,1,0, &
           !     dummylog,dummylog,dummylog,dummylog,1)
           !call create_color_image(ux(1:nx+1,1:nz+1),nx+1,nz+1,it,isx,isz,ix_rec,iz_rec,1,&
           !    NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,1)
           if(videoornot) then
              call create_color_image(uz(1:nx+1,1:nz+1),nx+1,nz+1,it,isx,isz, &
                   nrx(1:nReceiver),nrz(1:nReceiver),nReceiver, &
                   NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,2)
              
              
           endif
           
           !if(optimise) then
           !   write(outfile,'("video",I5,".",I5,".",I5,".OPT_UX") ') it,isx-lmargin(1),isz-lmargin(2)
           !else
           !   write(outfile,'("video",I5,".",I5,".",I5,".CON_UX") ') it,isx-lmargin(1),isz-lmargin(2)
           !endif
           !do j=1,24
           !   if(outfile(j:j).eq.' ') outfile(j:j)='0'
           !enddo
           
           !outfile = './synthetics/'//trim(modelname)//'/'//outfile
           !open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
           !video(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2))= &
           !     ux(lmargin(1)+1:nx+1-rmargin(1),lmargin(2)+1:nz+1-rmargin(2))
           !write(1,rec=1)  video(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2))
           !close(1,status='keep')
           
           
           
           !if(optimise) then
           !   write(outfile,'("video",I5,".",I5,".",I5,".OPT_UX") ') it,isx-lmargin(1),isz-lmargin(2)
           !else
           !   write(outfile,'("video",I5,".",I5,".",I5,".CON_UX") ') it,isx-lmargin(1),isz-lmargin(2)
           !endif
           !do j=1,24
           !   if(outfile(j:j).eq.' ') outfile(j:j)='0'
           !enddo
           
           !outfile = './synthetics/'//trim(modelname)//'/'//outfile
           !open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
           !video(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2))= &
           !     ux(lmargin(1)+1:nx+1-rmargin(1),lmargin(2)+1:nz+1-rmargin(2))
           !write(1,rec=1)  video(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2))
           !close(1,status='keep')
           
        endif

        
        
        !call compNRBC2(ux(1:nx+1,1:nz+1),ux1(1:nx+1,1:nz+1),ux2(1:nx+1,1:nz+1), &
        !     uz(1:nx+1,1:nz+1),uz1(1:nx+1,1:nz+1),uz2(1:nx+1,1:nz+1), CerjanRate, lmargin, rmargin,nx+1,nz+1)
        
     enddo

     !write(18,*) singleStrainDiagonal(:,:)


     if(videoornot) then
        
        if(optimise) then
           write(outfile,'("video",".",I5,".",I5,".OPT.mp4") ') isx-lmargin(1),isz-lmargin(2)
        else
           write(outfile,'("video",".",I5,".",I5,".CON.mp4") ') isx-lmargin(1),isz-lmargin(2)
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './videos/'//trim(modelname)//'/'//outfile
        
        
        commandline="ffmpeg -framerate 5 -pattern_type glob -i 'snapshots/*.png' -c:v libx264 -pix_fmt yuv420p "//outfile
     
     endif
     
     
     call system(commandline)
  
  
     do ir = 1,nReceiver
        if(optimise) then
           write(outfile,'(I5,".",I5,".",I5,".",I5,".OPT_UX") ') nrx(ir)-lmargin(1),nrz(ir)-lmargin(2), &
                isx-lmargin(1),isz-lmargin(2)
        else
           write(outfile,'(I5,".",I5,".",I5,".",I5,".CON_UX") ') nrx(ir)-lmargin(1),nrz(ir)-lmargin(2), &
                isx-lmargin(1),isz-lmargin(2)
        endif
        
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './synthetics/'//trim(modelname)//'/'//outfile
        open(1, file=outfile,status='unknown',form='formatted')
        do it=0,nt
           write (1,*) time(it),synx(it,ir)
        enddo
        close(1)
        !open(1,file=outfile,form='unformatted',access='direct',recl=kind(0e0)*(nt+1))
        !write(1,rec=1) synx(0:nt,ir)
        !close(1,status='keep')

        
        if(optimise) then
           write(outfile,'(I5,".",I5,".",I5,".",I5,".OPT_UZ") ') nrx(ir)-lmargin(1),nrz(ir)-lmargin(2), &
                isx-lmargin(1),isz-lmargin(2)
        else
           write(outfile,'(I5,".",I5,".",I5,".",I5,".CON_UZ") ') nrx(ir)-lmargin(1),nrz(ir)-lmargin(2), &
                isx-lmargin(1),isz-lmargin(2)
        endif
        
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './synthetics/'//trim(modelname)//'/'//outfile
        open(1, file=outfile,status='unknown',form='formatted')
        do it=0,nt
           write (1,*) time(it), synz(it,ir)
        enddo
        close(1)

        !open(1,file=outfile,form='unformatted',access='direct',recl=kind(0e0)*(nt+1))
        !write(1,rec=1) synz(0:nt,ir)
        !close(1,status='keep')


        
     enddo
     
  enddo




end program multipleSourcesOPT2D




subroutine datainit( nx,nz,ux )
  
  integer nx,nz
  double precision ux(nx+1,*)
  integer i,j
  
  do j=1,nz+1
     do i=1,nx+1
        ux(i,j) = 0.d0
     enddo
  enddo
  
  return
end subroutine datainit



 




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
  double precision ux(maxnz+1,*),ux1(maxnz+1,*),ux2(maxnz+1,*)
  double precision uz(maxnz+1,*),uz1(maxnz+1,*),uz2(maxnz+1,*)
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
  double precision fx(maxnz+1,*),fz(maxnz+1,*)
  double precision work1(maxnz+1,-2:1),work2(maxnz+1,-1:2)
  double precision work3(maxnz+1,-2:1),work4(maxnz+1,-1:2)
  double precision work5(*),work6(maxnz+1,0:1)
  double precision work7(*),work8(maxnz+1,0:1)
  double precision work9(*),work10(maxnz+1,-2:1)
  double precision work11(*),work12(maxnz+1,-1:2)
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
     ux(isx,isz) = ux(isx,isz) + fx(isx,isz)
     uz(isx,isz) = uz(isx,isz) + fz(isx,isz)
 endif

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
  double precision, parameter :: POWER_DISPLAY = 1.d0
  
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
  logical, parameter :: kurama =.true.
  integer :: ix,iy,irec
  
  character(len=100) :: file_name,system_command1,system_command2,system_command3
	
  integer :: R, G, B
  double precision :: R1, G1, B1
  
  double precision :: normalized_value,max_amplitude
  
  !       open image file and create system command to convert image to more convenient format
  !       use the "convert" command from ImageMagick http://www.imagemagick.org
  if(field_number == 1) then
     write(file_name,"('image',i6.6,'_Ux.pnm')") it
     write(system_command1, "('convert image',i6.6,'_Ux.pnm snapshots/imageUx',i6.6,'.png')") it,it
     write(system_command2, "('rm image',i6.6,'_Ux.pnm')") it
  else if(field_number == 2) then
     write(file_name,"('image',i6.6,'_Uz.pnm')") it
     write(system_command1,"('convert image',i6.6,'_Uz.pnm snapshots/imageUz',i6.6,'.png')") it,it
     write(system_command2,"('rm image',i6.6,'_Uz.pnm')") it
  endif
  
  open(unit=27, file=file_name, status='unknown')
  
  write(27,"('P3')")	! write image in PNM P3 format
  
  write(27,*) NX,NY	! write image size
  write(27,*) '255'	! maximum value of each pixel color
	
  !       compute maximum amplitude
  max_amplitude = maxval(abs(image_data_2D))
  
  !       image starts in upper-left corner in PNM format
  do iy=1,NY,1
     do ix=1,NX
        
        !       define data as vector component normalized to [-1:1] and rounded to nearest integer
        !       keeping in mind that amplitude can be negative
        normalized_value = image_data_2D(ix,iy) / max_amplitude *1.5d0
        
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
        elseif(kurama) then
           normalized_value=normalized_value**POWER_DISPLAY
           call plotcolor(normalized_value,R1,G1,B1)
           R = nint(R1)
           G = nint(G1)
           B = nint(B1)
        else
           if(normalized_value >= 0.d0) then
              R = 255
              G = nint(255.d0-255.d0*normalized_value**POWER_DISPLAY)
              B = nint(255.d0-255.d0*normalized_value**POWER_DISPLAY)
           else
              R = nint(255.d0-255.d0*abs(normalized_value)**POWER_DISPLAY)
              G = nint(255.d0-255.d0*abs(normalized_value)**POWER_DISPLAY)
              B = 255
           endif
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



subroutine  compNRBC2(ux,ux1,ux2,uz,uz1,uz2, rrate, lmargin, rmargin,nnx,nnz)

  ! Cerjan boundary conditions (2D)
  implicit none
  integer :: nnx, nnz

  real*8, intent(inout) :: ux2(nnx,nnz),uz2(nnx,nnz)
  real*8, intent(inout) :: ux1(nnx,nnz),uz1(nnx,nnz)
  real*8, intent(inout) :: ux(nnx,nnz),uz(nnx,nnz)
  real*8, intent(in) :: rrate
  integer, dimension(3), intent(in) :: lmargin, rmargin
  integer ix, iy, iz
  integer i, j, k
  real*8 r
  
  do iz = 1, nnz
     !do iy = 1, nny
     do ix = 1, nnx
        
        i = 0
        j = 0
        k = 0
           
        if (ix < lmargin(1) + 1) i = lmargin(1) + 1 - ix
        !   if (iy < lmargin(2) + 1) j = lmargin(2) + 1 - iy
        if (iz < lmargin(2) + 1) k = lmargin(2) + 1 - iz
        if (nnx - rmargin(1) < ix) i = ix - nnx + rmargin(1)
        !if (nny - rmargin(2) < iy) j = iy - nny + rmargin(2)
        if (nnz - rmargin(2) < iz) k = iz - nnz + rmargin(2)
           
        if (i == 0 .and. j == 0 .and. k == 0) cycle
        
        r = rrate * rrate * dble( i * i + j * j + k * k )
        r = exp( - r )
        

        ux2(:,:) = ux2(:,:) * r
        ux1(:,:) = ux1(:,:) * r
        ux(:,:) = ux(:,:) * r

        uz2(:,:) = uz2(:,:) * r
        uz1(:,:) = uz1(:,:) * r
        uz(:,:) = uz(:,:) * r
        
     enddo
     !enddo
  enddo
  
end subroutine compNRBC2



subroutine  compNRBCpre(r,rrate, lmargin, rmargin,nnx,nnz)

  implicit none
  integer :: nnx, nnz
  ! Cerjan boundary conditions (2D)
  double precision :: r(nnx,nnz)
  real*8, intent(in) :: rrate
  integer, dimension(3), intent(in) :: lmargin, rmargin
  integer ix, iy, iz
  integer i, j, k
  double precision :: rr
  
  
  do iz = 1, nnz
     !do iy = 1, nny
     do ix = 1, nnx
        
        i = 0
        j = 0
        k = 0
           
        if (ix < lmargin(1) + 1) i = lmargin(1) + 1 - ix
        !   if (iy < lmargin(2) + 1) j = lmargin(2) + 1 - iy
        if (iz < lmargin(2) + 1) k = lmargin(2) + 1 - iz
   
        if (nnx - rmargin(1) < ix) i = ix - nnx + rmargin(1)
        !if (nny - rmargin(2) < iy) j = iy - nny + rmargin(2)
        if (nnz - rmargin(2) < iz) k = iz - nnz + rmargin(2)
           
        if (i == 0 .and. j == 0 .and. k == 0) cycle
        
        rr = rrate * rrate * dble( i * i + j * j + k * k )
        r(ix,iz) = exp( - rr )
        
        !if(r(ix,iz).ne.1.d0) then
        !   print *, ix,iz,r(ix,iz)
        !endif
        
        !print *, ix,iy,r
        !ux2(:,:) = ux2(:,:) * r
        !ux1(:,:) = ux1(:,:) * r
        !ux(:,:) = ux(:,:) * r

        !uz2(:,:) = uz2(:,:) * r
        !uz1(:,:) = uz1(:,:) * r
        !uz(:,:) = uz(:,:) * r
        
     enddo

     
     
     !enddo
  enddo
  
end subroutine compNRBCpre






subroutine plotcolor(v,r1,g1,b1)
  !======================================================================
  ! Interpolate r1 g1 b1
  !               2016. 5. Kurama OKUBO IPGP
  !======================================================================
  implicit none
  real(8) v,v1,nl
  double precision r1, g1, b1
  real(8),allocatable,dimension(:) :: x,r,g,b
  
  integer n,i
  
  !read colormap data
  open (17, file='../colormap/colormap.dat', status='old')
  read (17, *) nl
  
  n = int(nl)
  allocate( x(n) )
  allocate( r(n) )
  allocate( g(n) )
  allocate( b(n) )
  
  do i = 1, int(nl)
     read (17, *) x(i), r(i), g(i), b(i)
  end do
  close(17)
  !rescale normalized value from [-1 1] to [0 1]
  !Depending on colormap.dat
  
  !when colormap domain is [0 1]
  !v1 = (1.0d0+v)/2.0d0
  
  !when colormap domain is [-1 1]
  v1 = v
  
  call colormap(v1,n,x,r,g,b,r1,g1,b1)
  
end subroutine plotcolor

subroutine colormap(v,n,x,r,g,b,r1,g1,b1)
  !======================================================================
  !v: normalized value [-1 1]
  !n: number of the data in colormap
  !x: ampritude in dataset
  !r: dataset of R
  !g: dataset of G
  !b: dataset of B
  !r1: interpolated R
  !g1 interpolated G
  !b1: interpolated B
  !
  !                  2016.5. Kurama OKUBO IPGP
  !----------------------------------------------------------------------
  
  implicit none
  integer i,n
  real(8) v,x(n),r(n),g(n),b(n)
  double precision ispline, r1, g1, b1
  real(8) b_spline(n),c_spline(n),d_spline(n)

  !for R
  call spline (x, r, b_spline, c_spline, d_spline,n)
  r1 = ispline(v, x, r, b_spline, c_spline, d_spline, n)

  !for G
  call spline (x, g, b_spline, c_spline, d_spline,n)
  g1 = ispline(v, x, g, b_spline, c_spline, d_spline, n)

  !for B
  call spline (x, b, b_spline, c_spline, d_spline,n)
  b1 = ispline(v, x, b, b_spline, c_spline, d_spline, n)


end subroutine colormap


subroutine spline (x, y, b, c, d, n)
  !======================================================================
  !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
  !  for cubic spline interpolation
  !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
  !  for  x(i) <= x <= x(i+1)
  !  Alex G: January 2010
  !----------------------------------------------------------------------
  !  input..
  !  x = the arrays of data abscissas (in strictly increasing order)
  !  y = the arrays of data ordinates
  !  n = size of the arrays xi() and yi() (n>=2)
  !  output..
  !  b, c, d  = arrays of spline coefficients
  !  comments ...
  !  spline.f90 program is based on fortran version of program spline.f
  !  the accompanying function fspline can be used for interpolation
  !======================================================================
  implicit none
  integer n
  double precision x(n), y(n), b(n), c(n), d(n)
  integer i, j, gap
  double precision h
  
  gap = n-1
  ! check input
  if ( n < 2 ) return
  if ( n < 3 ) then
     b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
     c(1) = 0.
     d(1) = 0.
     b(2) = b(1)
     c(2) = 0.
     d(2) = 0.
     return
  end if
  !
  ! step 1: preparation
  !
  d(1) = x(2) - x(1)
  c(2) = (y(2) - y(1))/d(1)
  do i = 2, gap
     d(i) = x(i+1) - x(i)
     b(i) = 2.0*(d(i-1) + d(i))
     c(i+1) = (y(i+1) - y(i))/d(i)
     c(i) = c(i+1) - c(i)
  end do
  !
  ! step 2: end conditions
  !
  b(1) = -d(1)
  b(n) = -d(n-1)
  c(1) = 0.0
  c(n) = 0.0
  if(n /= 3) then
     c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
     c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
     c(1) = c(1)*d(1)**2/(x(4)-x(1))
     c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
  end if
  !
  ! step 3: forward elimination
  !
  do i = 2, n
     h = d(i-1)/b(i-1)
     b(i) = b(i) - h*d(i-1)
     c(i) = c(i) - h*c(i-1)
  end do
  !
  ! step 4: back substitution
  !
  c(n) = c(n)/b(n)
  do j = 1, gap
     i = n-j
     c(i) = (c(i) - d(i)*c(i+1))/b(i)
  end do
  !
  ! step 5: compute spline coefficients
  !
  b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
  do i = 1, gap
     b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
     d(i) = (c(i+1) - c(i))/d(i)
     c(i) = 3.*c(i)
  end do
  c(n) = 3.0*c(n)
  d(n) = d(n-1)
end subroutine spline

function ispline(u, x, y, b, c, d, n)
  !======================================================================
  ! function ispline evaluates the cubic spline interpolation at point z
  ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
  ! where  x(i) <= u <= x(i+1)
  !----------------------------------------------------------------------
  ! input..
  ! u       = the abscissa at which the spline is to be evaluated
  ! x, y    = the arrays of given data points
  ! b, c, d = arrays of spline coefficients computed by spline
  ! n       = the number of data points
  ! output:
  ! ispline = interpolated value at point u
  !=======================================================================
  implicit none
  double precision ispline
  integer n
  double precision  u, x(n), y(n), b(n), c(n), d(n)
  integer i, j, k
  double precision dx
  
  ! if u is ouside the x() interval take a boundary value (left or right)
  if(u <= x(1)) then
     ispline = y(1)
     return
  end if
  if(u >= x(n)) then
     ispline = y(n)
     return
  end if
  
  !*
  !  binary search for for i, such that x(i) <= u <= x(i+1)
  !*
  i = 1
  j = n+1
  do while (j > i+1)
     k = (i+j)/2
     if(u < x(k)) then
        j=k
     else
        i=k
     end if
  end do
  !*
  !  evaluate spline interpolation
  !*
  dx = u - x(i)
  ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline
