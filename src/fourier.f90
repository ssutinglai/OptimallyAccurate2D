subroutine FourierDeallocate
  use paramFWI

  deallocate(strainFieldD)
  deallocate(strainFieldS)
  deallocate(synFieldX)
  deallocate(synFieldZ)
  return
end subroutine FourierDeallocate

subroutine FourierAllocate
  use paramFWI
  use parameters
  implicit none
  integer :: iFreq
    
  recl_size=kind(1.e0)*(boxnx+1)*(boxnz+1)
  recl_size_syn=(maxnt+1)*(nReceiver+1)*kind(0e0)


  nFreq = 1
 
  do while (nFreq < (maxnt+1))
     nFreq = nFreq*2
  enddo
  
  tlen = dt*dble(nFreq)

  ! to stablise the deconvolution

  nFreqStart=int(8.d-1*f0*tlen)
  nFreqStop=int(3.2d0*f0*tlen)

  ! faut enlever !!
  !nnFreq=nFreq

  nFreqStep=(nFreqStop-nFreqStart)/nnFreq
  
  allocate(nFreqSample(0:nnFreq-1))

  do iFreq=0,nnFreq-1
     nFreqSample(iFreq)=nFreqStart+iFreq*nFreqStep
  enddo

  recl_size_fft=2*kind(1.e0)*nnFreq*(boxnx+1)*(boxnz+1)

  allocate(singleStrainFieldD(0:nnFreq-1,1:boxnx+1,1:boxnz+1))
  allocate(singleStrainFieldS(0:nnFreq-1,1:boxnx+1,1:boxnz+1))

  allocate(strainFieldD(0:nnFreq-1,1:boxnx+1,1:boxnz+1))
  allocate(strainFieldS(0:nnFreq-1,1:boxnx+1,1:boxnz+1))
  allocate(backStrainFieldD(0:nnFreq-1,1:boxnx+1,1:boxnz+1))
  allocate(backStrainFieldS(0:nnFreq-1,1:boxnx+1,1:boxnz+1))

  allocate(synFieldX(0:nnFreq-1,1:nReceiver,1:nSource))
  allocate(synFieldZ(0:nnFreq-1,1:nReceiver,1:nSource))
  allocate(obsFieldX(0:nnFreq-1,1:nReceiver,1:nSource))
  allocate(obsFieldZ(0:nnFreq-1,1:nReceiver,1:nSource))

  return
end subroutine FourierAllocate

subroutine FourierAll
  use parameters
  use paramFWI
  implicit none
  integer :: iFreq
  double precision :: angfreq0
  double precision, allocatable :: angfreq(:)
  double complex, allocatable :: sourceFreq(:)
  double complex, allocatable :: tmptensorD(:,:,:),tmptensorS(:,:,:)
  double complex, allocatable :: tmpSynX(:,:),tmpSynZ(:,:)
  ! determination of frequency numbers

  allocate(sourceFreq(0:nnFreq-1))
  allocate(angfreq(0:nnFreq-1))

  angfreq0  = 2.d0*pi*f0

  ! Ricker wavelet in frequency domain

  do iFreq=0,nnFreq-1
      
     angfreq(iFreq)=pi/tlen*dble(nFreqSample(iFreq))
    
     sourceFreq(iFreq)=2.d0*angfreq(iFreq)**2/sqrt(pi)/angfreq0**3* &
          exp(-(angfreq(iFreq)/angfreq0)**2)*exp(cmplx(0.d0,-t0*angfreq(iFreq))) 
     !write(13,*) real(sourceFreq(iFreq)),imag(sourceFreq(iFreq))
  enddo
  

  

  synFieldX=cmplx(0.d0)
  synFieldZ=cmplx(0.d0)


  do iSource = 1,nSource
     isx = iisx(iSource)
     isz = iisz(iSource)

     allocate(tmptensorD(0:2*nFreq-1,1:boxnx+1,1:boxnz+1))
     allocate(tmptensorS(0:2*nFreq-1,1:boxnx+1,1:boxnz+1))
     
     tmptensorD=cmplx(0.d0)
     tmptensorS=cmplx(0.d0)
     
     strainFieldD=cmplx(0.d0)
     strainFieldS=cmplx(0.d0)
     singleStrainFieldD=cmplx(0.e0)
     singleStrainFieldS=cmplx(0.e0)
     

     do it=0,nt,IT_DISPLAY
        tmpsingleStrain=0.e0     
        if(optimise) then
           write(outfile,'("strainD",I5,".",I5,".",I5,".OPT_dat") ') it,isx,isz
        else
           write(outfile,'("strainD",I5,".",I5,".",I5,".CON_dat") ') it,isx,isz
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './strains/'//trim(modelname)//'/'//outfile
        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
        read(1,rec=1)  tmpsingleStrain(1:boxnx+1,1:boxnz+1)
        close(1,status='keep')
        
        tmptensorD(it,1:boxnx+1,1:boxnz+1) &
             = tmpsingleStrain(1:boxnx+1,1:boxnz+1)
        
        
        tmpsingleStrain=0.e0     
        if(optimise) then
           write(outfile,'("strainS",I5,".",I5,".",I5,".OPT_dat") ') it,isx,isz
        else
           write(outfile,'("strainS",I5,".",I5,".",I5,".CON_dat") ') it,isx,isz
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './strains/'//trim(modelname)//'/'//outfile
        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
        read(1,rec=1)  tmpsingleStrain(1:boxnx+1,1:boxnz+1)
        close(1,status='keep')
        
        tmptensorS(it,1:boxnx+1,1:boxnz+1) &
             = tmpsingleStrain(1:boxnx+1,1:boxnz+1)
        
     enddo
     
     

     ! Deconvolution of Ricker wavelet in frequency


     !do iz=0,maxnt
     !write(17,*) real(tmptensorD(iz,1,1))
     !enddo


     do iz=1,nz+1-rmargin(2)-lmargin(2)
        do ix=1,nx+1-rmargin(1)-lmargin(1)

          
           call FFT_double(nFreq,tmptensorD(0:2*nFreq-1,ix,iz),tlen)
           call FFT_double(nFreq,tmptensorS(0:2*nFreq-1,ix,iz),tlen)
           
           do iFreq = 0,nnFreq-1
              strainFieldD(iFreq,ix,iz)=tmptensorD(nFreqSample(iFreq),ix,iz)/sourceFreq(iFreq)
              strainFieldS(iFreq,ix,iz)=tmptensorS(nFreqSample(iFreq),ix,iz)/sourceFreq(iFreq)
           enddo
           
        enddo
     enddo


        
     deallocate(tmptensorD)
     deallocate(tmptensorS)

     
     !! Write strains

     singleStrainFieldD(0:nnFreq-1,1:boxnx,1:boxnz)=strainFieldD(0:nnFreq-1,1:boxnx,1:boxnz)
     singleStrainFieldS(0:nnFreq-1,1:boxnx,1:boxnz)=strainFieldS(0:nnFreq-1,1:boxnx,1:boxnz)
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
     write(1,rec=1)  singleStrainFieldD(0:nnFreq-1,1:boxnx+1,1:boxnz+1)
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
     write(1,rec=1)  singleStrainFieldS(0:nnFreq-1,1:boxnx+1,1:boxnz+1)
     close(1,status='keep')

     
     allocate(tmpSynX(0:2*nFreq-1,1:nReceiver))
     allocate(tmpSynZ(0:2*nFreq-1,1:nReceiver))
     
     if(iterationIndex.eq.0) then
        
        ! Reading OBS data
        
        
        if(trim(extentionOBSx).ne."9999") then 
           write(outfile,'(I5,".",I5,".") ')  &
                isx,isz
           
           do j=1,12
              if(outfile(j:j).eq.' ') outfile(j:j)='0'
           enddo
           
           outfile = trim(outfile)//trim(extentionOBSx)       
           outfile = trim(obsdir)//'/'//trim(outfile)
           
           
           open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
           read(1,rec=1) obsx(0:maxnt,1:nReceiver)
           close(1)
           
        endif

        if(trim(extentionOBSz).ne."9999") then
           
           
           write(outfile,'(I5,".",I5,".") ') &
                isx,isz
           
           
           do j=1,12
              if(outfile(j:j).eq.' ') outfile(j:j)='0'
           enddo
           
           outfile = trim(outfile)//trim(extentionOBSz)
           outfile = trim(obsdir)//'/'//trim(outfile)
           
           open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
           read(1,rec=1) obsz(0:maxnt,1:nReceiver)
           close(1)
           
        endif
        
        ! Reading OBS data done
        
         
        tmpSynX(0:maxnt,1:nReceiver)=obsx(0:maxnt,1:nReceiver)
      
        tmpSynZ(0:maxnt,1:nReceiver)=obsz(0:maxnt,1:nReceiver)
     
        do iReceiver=1,nReceiver
           call FFT_double(nFreq,tmpSynX(0:2*nFreq-1,iReceiver),tlen)
           call FFT_double(nFreq,tmpSynZ(0:2*nFreq-1,iReceiver),tlen)
           do iFreq = 0,nnFreq-1
              obsFieldX(iFreq,iReceiver,iSource)= &
                   tmpSynX(nFreqSample(iFreq),iReceiver)/sourceFreq(iFreq)
              obsFieldZ(iFreq,iReceiver,iSource)= &
                   tmpSynZ(nFreqSample(iFreq),iReceiver)/sourceFreq(iFreq)
           enddo
        enddo        
     endif


     if(iterationIndex.eq.0) then
        if(optimise) then
           write(outfile,'(I5,".",I5,".OPT_UX") ') isx,isz
        else
           write(outfile,'(I5,".",I5,".CON_UX") ') isx,isz
        endif
        
        do j=1,12
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './synthetics/'//trim(modelname)//'/'//outfile
     else
        if(optimise) then
           write(outfile,'(I5,".",I5,".OPT_UX.it",I3.3) ') isx,isz,iterationIndex
        else
           write(outfile,'(I5,".",I5,".CON_UX,it",I3.3) ') isx,isz,iterationIndex
        endif
        
        do j=1,12
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './synthetics/'//trim(modelname)//'/'//outfile
        
     endif
     
     synx=0.e0
     
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
     read(1,rec=1) synx(0:maxnt,1:nReceiver)
     close(1)
     
     
     tmpSynX(0:maxnt,1:nReceiver)=synx(0:maxnt,1:nReceiver)
     
     
     
     if(iterationIndex.eq.0) then
        if(optimise) then
           write(outfile,'(I5,".",I5,".OPT_UZ") ') isx,isz
        else
           write(outfile,'(I5,".",I5,".CON_UZ") ') isx,isz
        endif
        
        do j=1,12
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './synthetics/'//trim(modelname)//'/'//outfile
     else
        if(optimise) then
           write(outfile,'(I5,".",I5,".OPT_UZ.it",I3.3) ') isx,isz,iterationIndex
        else
           write(outfile,'(I5,".",I5,".CON_UZ,it",I3.3) ') isx,isz,iterationIndex
        endif
        
        do j=1,12
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './synthetics/'//trim(modelname)//'/'//outfile
        
     endif
     
     synz=0.e0
     
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
     read(1,rec=1) synz(0:maxnt,1:nReceiver)
     close(1)
     
     
     tmpSynZ(0:maxnt,1:nReceiver)=synz(0:maxnt,1:nReceiver)
     
 
     
     ! Deconvolution of Ricker wavelet in frequency
     
     
     !do iz=0,nFreq-1
     !write(16,*) real(tmpSynZ(iz,1))
     !enddo


     do iReceiver=1,nReceiver
        call FFT_double(nFreq,tmpSynX(0:2*nFreq-1,iReceiver),tlen)
        call FFT_double(nFreq,tmpSynZ(0:2*nFreq-1,iReceiver),tlen)
        do iFreq = 0,nnFreq-1
           synFieldX(iFreq,iReceiver,iSource)= &
                tmpSynX(nFreqSample(iFreq),iReceiver)/sourceFreq(iFreq)
           synFieldZ(iFreq,iReceiver,iSource)= &
                tmpSynZ(nFreqSample(iFreq),iReceiver)/sourceFreq(iFreq)
        enddo
     enddo
   

     !do iz=0,nnFreq-1
     !write(17,*) synFieldZ(:,1,1)
     !enddo
     


     !call IFFT_double(nFreq,tmpSynZ(0:2*nFreq-1,1),tlen)
     !do iz=0,nFreq-1
     !write(18,*) real(tmpSynZ(iz,1))
     !enddo

     !stop


     

     deallocate(tmpSynX)
     deallocate(tmpSynZ)

  enddo

  print *, "FFT done"
  deallocate(sourceFreq)
  deallocate(angfreq)
  !deallocate(strainFieldD)
  !deallocate(strainFieldS)
  !deallocate(synFieldX)
  !deallocate(synFieldZ)

end subroutine FourierAll


subroutine FFT_double(nFreq,cvec,tlen)
  ! this subroutine particularly performs IFFT of the given tensor and make a double tensor
  
  implicit none

  integer :: i,j,nFreq,n1,m1
  
  complex(kind(0d0)) :: cvec(0:2*nFreq-1)

  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)), parameter :: sqrtpi2piinverse=8.979356106d-2
  real(kind(0d0)) :: tlen,samplingHz
  
  samplingHz = dble(2*nFreq)/tlen
 
  call cdft(4*nFreq,cos(pi/(2*nFreq)),-sin(pi/(2*nFreq)), cvec(0:2*nFreq-1))

  cvec(:)=cvec(:)/dble(2*nFreq)
  return
end subroutine FFT_double


subroutine IFFT_double(nFreq,cvec,tlen)
  ! this subroutine particularly performs IFFT of the given tensor and make a double tensor
  
  implicit none

  integer :: i,j,nFreq,n1,m1
  
  complex(kind(0d0)) :: cvec(0:2*nFreq-1)

  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)), parameter :: sqrtpi2piinverse=8.979356106d-2
  real(kind(0d0)) :: tlen,samplingHz
  
 

  samplingHz = dble(2*nFreq)/tlen
  
  do i = 0, nFreq-1
     n1 = nFreq +i
     m1 = nFreq -i
     cvec(n1) = conjg(cvec(m1))
  enddo
  
  call cdft(4*nFreq,cos(pi/(2*nFreq)),sin(pi/(2*nFreq)), cvec(0:2*nFreq-1))
     
  !cvec(:)=cvec(:)*sqrtpi2piinverse*samplingHz

  return
end subroutine IFFT_double


subroutine cdft(n, wr, wi, c)
      
  integer :: n, i, j, k, l, m
  real(kind(0d0)) :: wr, wi, a(0 : n - 1), wmr, wmi, wkr, wki 
  real(kind(0d0)) ::wdr, wdi, ss, xr, xi
  complex(kind(0d0)) :: c(0:n/2-1)
  
  do i = 0, n/2-1
     a(2*i) = dble(c(i))
     a(2*i+1) = imag(c(i))
  enddo
  

  wmr = wr
  wmi = wi
  m = n
  do while (m .gt. 4)
     l = m / 2
     wkr = 1
     wki = 0
     wdr = 1 - 2 * wmi * wmi
     wdi = 2 * wmi * wmr
     ss = 2 * wdi
     wmr = wdr
     wmi = wdi
     do j = 0, n - m, m
        i = j + l
        xr = a(j) - a(i)
        xi = a(j + 1) - a(i + 1)
        a(j) = a(j) + a(i)
        a(j + 1) = a(j + 1) + a(i + 1)
        a(i) = xr
        a(i + 1) = xi
        xr = a(j + 2) - a(i + 2)
        xi = a(j + 3) - a(i + 3)
        a(j + 2) = a(j + 2) + a(i + 2)
        a(j + 3) = a(j + 3) + a(i + 3)
        a(i + 2) = wdr * xr - wdi * xi
        a(i + 3) = wdr * xi + wdi * xr
     enddo
     do k = 4, l - 4, 4
        wkr = wkr - ss * wdi
        wki = wki + ss * wdr
        wdr = wdr - ss * wki
        wdi = wdi + ss * wkr
        do j = k, n - m + k, m
           i = j + l
           xr = a(j) - a(i)
           xi = a(j + 1) - a(i + 1)
           a(j) = a(j) + a(i)
           a(j + 1) = a(j + 1) + a(i + 1)
           a(i) = wkr * xr - wki * xi
           a(i + 1) = wkr * xi + wki * xr
           xr = a(j + 2) - a(i + 2)
           xi = a(j + 3) - a(i + 3)
           a(j + 2) = a(j + 2) + a(i + 2)
           a(j + 3) = a(j + 3) + a(i + 3)
           a(i + 2) = wdr * xr - wdi * xi
           a(i + 3) = wdr * xi + wdi * xr
        enddo
     enddo
     m = l
  enddo
  if (m .gt. 2) then
     do j = 0, n - 4, 4
        xr = a(j) - a(j + 2)
        xi = a(j + 1) - a(j + 3)
        a(j) = a(j) + a(j + 2)
        a(j + 1) = a(j + 1) + a(j + 3)
        a(j + 2) = xr
        a(j + 3) = xi
     enddo
  endif
  if (n .gt. 4) call bitrv2(n, a)

  
  do i = 0, n/2-1
     c(i) = dcmplx(a(2*i), a(2*i+1))
  enddo
  
  
end subroutine cdft



subroutine bitrv2(n, a)
  integer :: n, j, j1, k, k1, l, m, m2, n2
  real(kind(0d0)) :: a(0 : n - 1), xr, xi
  
  m = n / 4
  m2 = 2 * m
  n2 = n - 2
  k = 0
  do j = 0, m2 - 4, 4
     if (j .lt. k) then
        xr = a(j)
        xi = a(j + 1)
        a(j) = a(k)
        a(j + 1) = a(k + 1)
        a(k) = xr
        a(k + 1) = xi
     else if (j .gt. k) then
        j1 = n2 - j
        k1 = n2 - k
        xr = a(j1)
        xi = a(j1 + 1)
        a(j1) = a(k1)
        a(j1 + 1) = a(k1 + 1)
        a(k1) = xr
        a(k1 + 1) = xi
     endif
     k1 = m2 + k
     xr = a(j + 2)
     xi = a(j + 3)
     a(j + 2) = a(k1)
     a(j + 3) = a(k1 + 1)
     a(k1) = xr
     a(k1 + 1) = xi
     l = m
     do while (k .ge. l)
        k = k - l
        l = l / 2
     enddo
     k = k + l
  enddo
end subroutine bitrv2
