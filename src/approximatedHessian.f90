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
  double complex :: tmpfrechet(0:nFreq-1,1:2,1:nReceiver,1:nSource,&
       1:nx-rmargin(1)-lmargin(1),1:nz+1-rmargin(2)-lmargin(2))
  double complex :: deltad(0:nFreq-1,1:nReceiver,1:nSource) ! for vertical components for the moment
  integer :: iTypeParam,jTypeParam ! 1 for Vp and 2 for Vs
  integer :: iFreq
  integer :: ixz
  integer :: jxz,jx,jz

  tmpfrechet=cmplx(0.d0)
  deltad=cmplx(0.d0)

  do iSource=1,nSource
     do iReceiver=1,nReceiver
        do iFreq=0,nFreq-1
           deltad(iFreq,iRecever,iSource)= &               
                obsFieldZ(iFreq,iReceiver,iSource)-synFieldZ(iFreq,iReceiver,iSource)
        enddo
     enddo
  enddo


  do iz=lmargin(2)+1,nz+1-rmargin(2)
     do ix=lmargin(1)+1,nx+1-rmargin(1)
        do iSource=1,nSource
           do iRecever=1,nReceiver
              do iFreq=0,nFreq-1 

                 tmpfrechet(iFreq,1,iReceiver,iSource,ix-lmargin(1),iz-lmargin(2))= &
                      !tmpfrechet(iFreq,1,iReceiver,iSource,ix-lmargin(1),iz-lmargin(2))+ &
                      2.d0*rho(ix,iz)*vp(ix,iz)* &
                      strainFieldD(iFreq,ix-lmargin(1),iz-lmargin(2),iSource)* &
                      conjg(strainFieldD(iFreq,ix-lmargin(1),iz-lmargin(2),iReceiver)
                 
                 
                 tmpfrechet(iFreq,2,iReceiver,iSource,ix-lmargin(1),iz-lmargin(2))= &
                      !tmpfrechet(iFreq,1,iReceiver,iSource,ix-lmargin(1),iz-lmargin(2))+ &
                      2.d0*rho(ix,iz)*vs(ix,iz)* &
                      strainFieldS(iFreq,ix-lmargin(1),iz-lmargin(2),iSource)* &
                      conjg(strainFieldS(iFreq,ix-lmargin(1),iz-lmargin(2),iReceiver)
              enddo
           enddo
        enddo
     enddo
  enddo

  ata=0.d0
  atd=0.d0


  do ixz=1,(nx+1)*(nz+1)
     iz=(ixz-1)/(nx+1)+1
     ix=mod(ixz-1,nx+1)+1
     
     do iTypeParam=1,2
     


        do jxz=1,(nx+1)*(nz+1)
           jz=(jxz-1)/(nx+1)+1
           jx=mod(jxz-1,nx+1)+1

           do jTypeParam=1,2
             



              

           enddo

        enddo


     enddo

  enddo
  

  

  return

end subroutine approximatedHessian

