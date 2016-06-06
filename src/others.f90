subroutine paramMultiReader
  use parameters
  implicit none

  
  ! Reading Inf File
110 format(a80)
  read(5,110) modelname
  read(5,110) vpfile
  read(5,110) vsfile
  read(5,110) rhofile
  read(5,'(L1)') optimise
  read(5,'(L1)') videoornot
  read(5,'(L1)') writingStrain
  read(5,*) IT_DISPLAY
  read(5,*) iSourceStart,iSourceInterval,nSource
  read(5,*) iReceiverStart,iReceiverInterval,nReceiver
  read(5,*) nt,nx,nz
  read(5,*) dt,dx,dz
  read(5,*) f0, t0


  nt=nt*times
  nx=nx*times
  nz=nz*times
  dt=dt/dble(times)
  dx=dx/dble(times)
  dz=dz/dble(times)
  
  maxnt = nt
  maxnx = nx+(lmargin(1)+rmargin(1))*times
  maxnz = nz+(lmargin(2)+rmargin(2))*times

  
 
  call system('mkdir ./inffile')
   
  commandline="mkdir synthetics"
  call system(commandline)
  commandline="mkdir snapshots"
  call system(commandline)
  commandline="mkdir videos"
  call system(commandline)
  commandline="mkdir synthetics/"//trim(modelname)
  call system(commandline)
  commandline="mkdir videos/"//trim(modelname)
  call system(commandline)
  commandline="mkdir strains"
  call system(commandline)
  commandline="mkdir strains/"//trim(modelname)
  call system(commandline)
end subroutine paramMultiReader



subroutine ReceiverSourcePositions

  use parameters
  implicit none
  

  
  do iReceiver = 1, nReceiver
     nrx(iReceiver)=(iReceiverStart-1)*times+1+iReceiverInterval*times*(iReceiver-1)
     nrz(iReceiver)=(130-1)*times+1
  enddo
  
  do iSource = 1, nSource
     iisx(iSource)=(iSourceStart-1)*times+1+iSourceInterval*times*(iSource-1)
     iisz(iSource)=(100-1)*times+1
     write(filename, '(I5,".",I5,".inf")') iisx(iSource),iisz(iSource)
     do j=1, 12
        if(filename(j:j).eq.' ') filename(j:j)='0'
     enddo
    filename= './inffile/'//trim(modelname)//'.'//filename
    
       

     open(1, file=filename, form='formatted')
     write(1,'(a)') 'c grids'
     write(1,*) nt, nx, nz
     write(1,*) dt, dx, dz
     write(1,'(a)') 'c modelname'
     write(1,'(a)') trim(modelname)
     write(1,'(a)') trim(vpfile)
     write(1,'(a)') trim(vsfile)
     write(1,'(a)') trim(rhofile)

     write(1,'(a)') 'c source position (in grids) '
     write(1,*) iisx(iSource),iisz(iSource)
     write(1,'(a)') 'c source time function (Ricker wavelet)'
     write(1,*) f0,t0     
     write(1,'(a)') 'c receivers information'
     write(1,*) nReceiver
     do iReceiver = 1, nReceiver
        write(1,*) nrx(iReceiver), nrz(iReceiver)
     enddo
     

     !write(1,*) 'c'
     write(1,'(a)') 'end'     
  enddo

end subroutine ReceiverSourcePositions
