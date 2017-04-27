subroutine FourierAll
  use parameters
  use paramFWI
  implicit none
  
  ! determination of frequency numbers
  
  nFreq = 1
  
  do while (nFreq < (maxnt+1))
     nFreq = nFreq*2
  enddo
     

  ! record size
  recl_size=(nx+1-rmargin(1)-lmargin(1))*(nz+1-rmargin(2)-lmargin(2))*kind(0e0)
  recl_size_syn=(maxnt+1)*(nReceiver+1)*kind(0e0)
  recl_size_strain=recl_size*nFreq*kind(cmplx(0d0))


  allocate(strainFieldP(0:nFreq,1:nx+1-rmargin(1)-lmargin(1), &
       1:nz+1-rmargin(2)-lmargin(2)))
