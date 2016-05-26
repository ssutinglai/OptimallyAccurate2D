

	program analitic
! c
! c	program to calculate an analitical solution
! c	for accoustic homogenious medium with a free surface
! c	in 2D
! c
        implicit none
	character*1 :: ans
	character*60 :: gfile,wfile,rfile
        real, dimension(10000) :: g,wvlt
	real :: vel,sdepth,rdepth,srdist,si,recl
	real :: r,rc,rc2,t,tmrc
	integer :: i,nsamp,ns2,ns2b,lwf,np2

! c
! cccc	get information
! c
! c
! cccc	velocity
! c
	write(*,'(''Enter velocity : ''$)')
	read(*,*) vel
!c
!cccc	geometry of source and receiver
!c
	write(*,'(''Enter depth of source (meters) : ''$)')
	read(*,*) sdepth
	write(*,'(''Enter depth of receiver (meters): ''$)')
	read(*,*) rdepth
	write(*,'(''Enter horizontal distance src-rec (meters): ''$)')
	read(*,*)srdist
	write(*,'(''Enter sampling rate in sec. : ''$)')
	read(*,*)si
	write(*,'(''Enter total recording lenght in sec. : ''$)')
	read(*,*)recl
!c
!cccc	output file
!c
	write(*,'(''Enter name of output green function : ''$)')
	read(*,'(a60)')gfile

!c
!cccc	calculate distance source-receiver for direct arrival
!c
	r = sqrt(srdist**2. + (rdepth - sdepth)**2.)
	write(*,'(''Direct arrival src-rec : '',f6.1,'' meters'')')r
	do 10 i = 1,10000
	g(i) = 0.
 10	continue
	nsamp = int(recl / si )+1
	rc = r / vel
	rc2 = rc**2.
!c
	do 20 i = 1,nsamp
	t = si * float(i-1)
	tmrc = t - rc
	if (tmrc .ge. si) then
	g(i) =  1. / sqrt( t**2 - rc2)
	end if
 20	continue
!c
!cccc	free surface
!c
	write(*,'(''Do you want a free surface (y/n) : ''$)')
	read(*,'(a1)')ans
	if (ans .eq. 'y') then
!c
!cccc	calculate distance virt. src. - rec.
!c
	r = sqrt(srdist**2 + (rdepth + sdepth)**2)
	write(*,'(''Refl. arrival src-rec : '',f6.1,'' meters'')')r
	rc = r / vel
	rc2 = rc**2
!c
	do 30 i = 1,nsamp
	t = si * float(i-1)
	tmrc = t - rc
	if (tmrc .gt. 0.) then
	g(i) = g(i) - (1. / sqrt( t**2 - rc2))
	end if
 30	continue
	end if
!c
!c	write green's function to file
!c
	open(unit=10,file=gfile,status='unknown', &
       form='unformatted',access='direct',recl=nsamp)
	write(10,rec=1)(g(i),i=1,nsamp)
	write(11,*)(g(i),i=1,nsamp)
	close(unit=10)
!c
	write(*,'(''Green function written to file '',a60)')gfile
	write(*,'(''Number of samples : '',i4)')nsamp
!c
!cccc	convlove with source  ???
!c
	write(*,'(''Convolve with a source ? ''$)')
	read(*,'(a1)')ans
	if (ans .eq. 'y') then
	write(*,'(''Enter name of src. funct. file : ''$)')
	read(*,'(a60)')wfile
	write(*,*) 'Enter number of samples'
        read(*,*)lwf
	open(unit=10,file=wfile,status='unknown', &
      form='unformatted',access='direct',recl=lwf)
	read(10,rec=1)(wvlt(i),i=1,lwf)
	close(unit=10)
!c
!cccc	find power of two
!c
	do 50 i=2,13
	np2=i
	ns2b = 2 ** np2
	if (ns2b .gt. nsamp+nsamp) go to 60
 50	continue
 60	continue
	ns2 = ns2b /2
!c------ fft of green
	call fft(g,np2,1,1)
	call fft(wvlt,np2,1,1)
	call cmul(g,wvlt,g,ns2)
	call fft(g,np2,-1,1)
	write(*,'(''Enter name of record file : ''$)')
	read(*,'(a60)')rfile
	open(unit=10,file=rfile,status='unknown', &
       form='unformatted',access='direct',recl=nsamp)
	write(10,rec=1)(g(i),i=1,nsamp)
	close(unit=10)
	end if

!c
!c	that's all
!c
	stop
	end  			
