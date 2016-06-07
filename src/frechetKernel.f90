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

  call paramMultiReader

  ! nx and nz here are not considering boundary conditions
  allocate(singleStrainForward(1:nx+1,1:nz+1))
  allocate(singleStrainBack(1:nx+1,1:nz+1))
  allocate(strainForward(1:nx+1,1:nz+1))
  allocate(strainBack(1:nx+1,1:nz+1))

  
  
  t=0.d0
  time(0)=t
  
  do it = 0, nt
     if(writingStrain.and.(mod(it,IT_DISPLAY).eq.0)) then
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
        read(1,rec=1) tmpsingleStrain
        
  enddo


end program frechetKernel
  
