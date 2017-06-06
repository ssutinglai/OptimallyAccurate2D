subroutine approximatedHessian

  

  ! Computation of Frechet derivatives for approximated Hessian for Vp and Vs
  !                                    and gradient direction for Vp and Vs
  ! Strain wavefields should be calculated beforehand
  !

  !
  !
  !                                            2017.6. N. Fuji (IPGP)
  !
  !


  use parameters
  use paramFWI
  !use paramFrechet
  implicit none

  !double complex :: tmpfrechet(0:nFreq-1,1:2,1:nReceiver,1:nSource,&
  !     1:nx+1-rmargin(1)-lmargin(1),1:nz+1-rmargin(2)-lmargin(2))
  !double complex :: deltad(0:nFreq-1,1:nReceiver,1:nSource) ! for vertical components for the moment
  integer :: iTypeParam,jTypeParam ! 1 for Vp and 2 for Vs
  integer :: iFreq
  integer :: ixz
  integer :: jxz,jx,jz
  integer :: jxzlocal
  double complex :: tmpfrechet1,tmpfrechet2

  tmpfrechet1=cmplx(0.d0)
  tmpfrechet2=cmplx(0.d0)
  !deltad=cmplx(0.d0)
  
  ata=cmplx(0.d0)
  atd=cmplx(0.d0)



  print *, "start approximated Hessian"

  synFieldZ(:,:,:)=obsFieldZ(:,:,:)-synFieldZ(:,:,:)



     
   do iSource=1,nSource
      isx = iisx(iSource)
      isz = iisz(iSource)
      if(optimise) then
         write(outfile,'("FourierStrainD",".",I5,".",I5,".OPT_dat") ') isx,isz
      else
         write(outfile,'("FourierStrainD",".",I5,".",I5,".CON_dat") ') isx,isz
      endif
      do j=1,26
         if(outfile(j:j).eq.' ') outfile(j:j)='0'
      enddo
      
      outfile = './strains/'//trim(modelname)//'/'//outfile
      open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_fft)
      read(1,rec=1)  singleStrainFieldD(0:nnFreq-1,1:boxnx+1,1:boxnz+1)
      close(1,status='keep')
      
      
      if(optimise) then
         write(outfile,'("FourierStrainS",".",I5,".",I5,".OPT_dat") ') isx,isz
      else
         write(outfile,'("FourierStrainS",".",I5,".",I5,".CON_dat") ') isx,isz
      endif
      do j=1,26
         if(outfile(j:j).eq.' ') outfile(j:j)='0'
      enddo
      
      outfile = './strains/'//trim(modelname)//'/'//outfile
      open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_fft)
      read(1,rec=1)  singleStrainFieldS(0:nnFreq-1,1:boxnx+1,1:boxnz+1)
      close(1,status='keep')
      
      strainFieldD=cmplx(0.d0)
      strainFieldS=cmplx(0.d0)
        
      strainFieldD(0:nnFreq-1,1:boxnx,1:boxnz)=singleStrainFieldD(0:nnFreq-1,1:boxnx,1:boxnz)
      strainFieldS(0:nnFreq-1,1:boxnx,1:boxnz)=singleStrainFieldS(0:nnFreq-1,1:boxnx,1:boxnz)
      
      do iReceiver=1,nReceiver
         
         if(optimise) then
            write(outfile,'("FourierStrainD",".",I5,".",I5,".OPT_dat") ') & 
                 nrx(iReceiver)-lmargin(1),nrz(iReceiver)-lmargin(2)
         else
            write(outfile,'("FourierStrainD",".",I5,".",I5,".CON_dat") ') &
                 nrx(iReceiver)-lmargin(1),nrz(iReceiver)-lmargin(2)
         endif
         do j=1,26
            if(outfile(j:j).eq.' ') outfile(j:j)='0'
         enddo
         
         outfile = './strains/'//trim(modelname)//'/'//outfile
         open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_fft)
         read(1,rec=1)  singleStrainFieldD(0:nnFreq-1,1:boxnx+1,1:boxnz+1)
         close(1,status='keep')
         
         print *, "iSource, iReceiver=", iSource, iReceiver

         if(optimise) then
            write(outfile,'("FourierStrainS",".",I5,".",I5,".OPT_dat") ') & 
                 nrx(iReceiver)-lmargin(1),nrz(iReceiver)-lmargin(2) 
         else
            write(outfile,'("FourierStrainS",".",I5,".",I5,".CON_dat") ')  &
                 nrx(iReceiver)-lmargin(1),nrz(iReceiver)-lmargin(2)
         endif
         do j=1,26
            if(outfile(j:j).eq.' ') outfile(j:j)='0'
         enddo
         
         outfile = './strains/'//trim(modelname)//'/'//outfile
         open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_fft)
         read(1,rec=1)  singleStrainFieldS(0:nnFreq-1,1:boxnx+1,1:boxnz+1)
         close(1,status='keep')
         
         backStrainFieldD=cmplx(0.d0)
         backStrainFieldS=cmplx(0.d0)
         
         backStrainFieldD(0:nnFreq-1,1:boxnx,1:boxnz)=singleStrainFieldD(0:nnFreq-1,1:boxnx,1:boxnz)
         backStrainFieldS(0:nnFreq-1,1:boxnx,1:boxnz)=singleStrainFieldS(0:nnFreq-1,1:boxnx,1:boxnz)
           

         !write(14,*) singleStrainFieldS(:,boxnx/2,boxnz/2)
        
         
         do ixz=1,(boxnx+1)*(boxnz+1)
            iz=(ixz-1)/(boxnx+1)+1
            ix=mod(ixz-1,boxnx+1)+1
            do iTypeParam=1,2
               do iFreq=0,nnFreq-1 
                  call frechet1point(iFreq,iTypeParam,ix,iz,tmpfrechet1)
                  atd(2*(ixz-1)+iTypeParam)= &
                       atd(2*(ixz-1)+iTypeParam)+ &
                       conjg(tmpfrechet1)* &   
                       synFieldZ(iFreq,iReceiver,iSource)
                  
                  do jz=max(iz-nNeighbours/2,1),min(iz+nNeighbours/2,boxnz+1)
                     do jx=max(ix-nNeighbours/2,1),min(ix+nNeighbours/2,boxnx+1)
                        jxz=(jz-1)*(boxnx+1)+jx
                        jxzlocal=(jz-iz+nNeighbours/2)*nNeighbours+(jx-ix+nNeighbours/2+1)
                        !print *, ix,iz,jx,jz,jxzlocal
                        
                        
                        do jTypeParam=1,2
                           
                           call frechet1point(iFreq,jTypeParam,jx,jz,tmpfrechet2)
                           
                          
                           ata(2*(jxzlocal-1)+jTypeParam,2*(ixz-1)+iTypeParam)= &
                                ata(2*(jxzlocal-1)+jTypeParam,2*(ixz-1)+iTypeParam)+ &
                                tmpfrechet1*conjg(tmpfrechet2)
                           !! AtA here is the conjugate transpose !!

                        enddo
                        
                     enddo
                     
                  enddo
                 
              enddo
           enddo
        enddo
     enddo
  enddo
  
  print *, "end Hessian calculation"
  
  open(unit=1,form="unformatted",file='ataatd')
  write(1) ata,atd
  close(1)
  return

end subroutine approximatedHessian

subroutine frechet1point(iFreq,iTypeParam,indexx,indexz,tmpfrechet)
  use parameters
  use paramFWI
  implicit none
  integer :: iFreq, iTypeParam,indexx,indexz
  double complex :: tmpfrechet


  if(iTypeParam.eq.1) then

     tmpfrechet= &
          2.d0*rho(indexx+lmargin(1),indexz+lmargin(2))*vp(indexx+lmargin(1),indexz+lmargin(2))* &
          strainFieldD(iFreq,indexx,indexz)* &
          conjg(backStrainFieldD(iFreq,indexx,indexz))
     
  elseif(iTypeParam.eq.2) then
  
     tmpfrechet= &
          2.d0*rho(indexx+lmargin(1),indexz+lmargin(2))*vs(indexx+lmargin(1),indexz+lmargin(2))* &
          strainFieldS(iFreq,indexx,indexz)* &
          conjg(backStrainFieldS(iFreq,indexx,indexz))
  endif
  
end subroutine frechet1point
