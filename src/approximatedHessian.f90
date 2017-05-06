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
  double complex :: tmpfrechet(0:nFreq-1,1:nReceiver,1:nSource,&
       1:nx-rmargin(1)-lmargin(1),1:nz+1-rmargin(2)-lmargin(2),1:2)
  integer :: iTypeObservable ! 1 for Vp and 2 for Vs
  integer :: iFreq

  
  do iTypeObservable=1,2
     do iz=1,nz+1-rmargin(2)-lmargin(2)
        do ix=1,nx+1-rmargin(1)-lmargin(1)
           do iSource=1,nSource
              do iRecever=1,nReceiver
                 do iFreq= 0,nFreq-1 
                    
  
  
  

  

  return

end subroutine approximatedHessian

