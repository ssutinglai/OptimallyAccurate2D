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
  double complex :: tmpfrechet1,tmpfrechet2

  tmpfrechet1=cmplx(0.d0)
  tmpfrechet2=cmplx(0.d0)
  !deltad=cmplx(0.d0)
  
  ata=0.d0
  atd=0.d0


  open(unit=1,form="unformatted",file="./tmpbinary")
  read(1) strainFieldS, strainFieldD, obsFieldZ,synFieldZ
  close(1)


  print *, "start approximated Hessian"


  do ixz=1,(boxnx+1)*(boxnz+1)
     iz=(ixz-1)/(boxnx+1)+1
     ix=mod(ixz-1,boxnx+1)+1
     
     
     do iSource=1,nSource
        do iReceiver=1,nReceiver
           do iFreq=0,nFreq-1
              do iTypeParam=1,2
                 call frechet1point(iFreq,iTypeParam,ix,iz,tmpfrechet1)
                 atd(2*(ixz-1)+iTypeParam)= &
                      atd(2*(ixz-1)+iTypeParam)+ &
                      conjg(tmpfrechet1)* &   
                      (obsFieldZ(iFreq,iReceiver,iSource)-synFieldZ(iFreq,iReceiver,iSource))
                      !deltad(iFreq,iReceiver,iSource)
                 do jxz=1,(boxnx+1)*(boxnz+1)
                    jz=(jxz-1)/(boxnx+1)+1
                    jx=mod(jxz-1,boxnx+1)+1
                    
                    do jTypeParam=1,2
                       
                       call frechet1point(iFreq,jTypeParam,jx,jz,tmpfrechet2)
                       
                       
                       ata(2*(ixz-1)+iTypeParam,2*(jxz-1)+jTypeParam)= &
                            ata(2*(ixz-1)+iTypeParam,2*(jxz-1)+jTypeParam)+ &
                            conjg(tmpfrechet1)*tmpfrechet2
              
                    enddo
                    
                 enddo

              enddo
           enddo
        enddo
     enddo
  enddo
  
  print *, "end Hessian calculation"
  

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
          strainFieldD(iFreq,indexx,indexz,iSource)* &
          conjg(strainFieldD(iFreq,indexx,indexz,iReceiver))
     
  elseif(iTypeParam.eq.2) then
  
     tmpfrechet= &
          2.d0*rho(indexx+lmargin(1),indexz+lmargin(2))*vs(indexx+lmargin(1),indexz+lmargin(2))* &
          strainFieldS(iFreq,indexx,indexz,iSource)* &
          conjg(strainFieldS(iFreq,indexx,indexz,iReceiver))
  endif
  
end subroutine frechet1point
