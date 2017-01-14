program multipleSourcesFWI2D

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
  use paramFWI
  implicit none


  ! Cerjan boundary
  lmargin(1)=NPOINTS_PML
  rmargin(1)=NPOINTS_PML
  lmargin(2)=NPOINTS_PML
  rmargin(2)=NPOINTS_PML
  
  call paramFWIReader

  call vectorAllocate
  
  call vectorAllocateFWI

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

  ! structuring absorbing boundary

  call calstructBC(maxnx, maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)


  ! first forward modelling

  call forwardmodelling
     
  iterationIndex=1
  
  if(iterationIndex<numberIteration) then
     iterationIndex=iterationIndex+1
     call backpropagation
     
     call gradientCalculation
     
     write(outfile,'("./iteratedModels/",I5,".vpgrad")'),iterationIndex
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
     write(1,rec=1) kernelP
     close(1)

 
     write(outfile,'("./iteratedModels/",I5,".vsgrad")'),iterationIndex
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
     write(1,rec=1) kernelS
     close(1)

     vp(:,:) = vp(:,:) + steplengthVp * kernelP(:,:)
     vs(:,:) = vs(:,:) + steplengthVs * kernelS(:,:)
     call calstruct2(maxnx,maxnz,nx,nz,rho,vp,vs,lam,mu,liquidmarkers)
     call calstructBC(maxnx, maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)
     call forwardmodelling

     synx(:,:) = obsx(:,:)-synx(:,:)
     synz(:,:) = obsz(:,:)-synz(:,:)

     ! here, syn is no more syn !!!

     numeratorG = sum(synx(:,:)*delx(:,:))+sum(synz(:,:)*delx(:,:))
     denominatorG = sum(synx(:,:)*synx(:,:))+sum(synz(:,:)*synz(:,:))

     alphaVp = numeratorG/denominatorG*steplengthVp
     alphaVs = numeratorG/denominatorG*steplengthVs


     ! new model construction
     
     vp(:,:) = vp(:,:) + (alphaVp-steplengthVp) * kernelP(:,:)
     vs(:,:) = vs(:,:) + (alphaVs-steplengthVs) * kernelS(:,:)


     write(outfile,'("./iteratedModels/",I5,".vpmodel")'),iterationIndex
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
     write(1,rec=1) vp
     close(1)

 
     write(outfile,'("./iteratedModels/",I5,".vsmodel")'),iterationIndex
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
     write(1,rec=1) vs
     close(1)


     call calstruct2(maxnx,maxnz,nx,nz,rho,vp,vs,lam,mu,liquidmarkers)
     call calstructBC(maxnx, maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)
     call forwardmodelling
     
     

  endif
     

     
     
  

end program multipleSourcesFWI2D




