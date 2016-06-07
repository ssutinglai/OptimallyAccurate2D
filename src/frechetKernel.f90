program frechetKernel

  

  ! Computation of Frechet derivatives for Vp and Vs
  ! Strain wavefields should be calculated beforehand with OPT2D.f90
  !

  !
  !
  !                                            2016.6. N. Fuji (IPGP)
  !
  !


  use parameters
  use paramFrechet
  implicit none

  call paramFrechetReader
  call vectorAllocateFrechet
  call ReceiverSourcePositions

  isx1 = iisx(i1Source)
  isz1 = iisz(i1Source)
  isx2 = iisx(i2Source)
  isz2 = iisz(i2Source)

  do it = 0, nt
     time(it)=dt*dble(it)
  enddo

  recl_size=kind(1.0)*(nx+1)*(nz+1)
  
  do it = 0,nt,IT_DISPLAY

     kernelP = 0.d0
     
     do it1 = 0,nt,IT_DISPLAY
     
        it2 = it1+it
        
        if(it2.gt.nt) cycle

        if(optimise) then
           write(outfile,'("strain",I5,".",I5,".",I5,".OPT_dat") ') it1,isx1,isz1
        else
           write(outfile,'("strain",I5,".",I5,".",I5,".CON_dat") ') it1,isx1,isz1
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './strains/'//trim(modelname)//'/'//outfile
        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
        read(1,rec=1) singleStrainForward
        
        StrainForward(:,:) = singleStrainForward(:,:)
        
  
        
        if(optimise) then
           write(outfile,'("strain",I5,".",I5,".",I5,".OPT_dat") ') it2,isx2,isz2
        else
           write(outfile,'("strain",I5,".",I5,".",I5,".CON_dat") ') it2,isx2,isz2
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './strains/'//trim(modelname)//'/'//outfile
        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
        read(1,rec=1) singleStrainBack
        
        StrainBack(:,:) = singleStrainBack(:,:)

        kernelP= kernelP+IT_DISPLAY*dble(dt)*(StrainForward*StrainBack)

     enddo   
     singleKernelP = kernelP

     
     
  enddo


end program frechetKernel
  
