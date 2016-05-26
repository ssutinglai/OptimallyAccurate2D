	program gauss
	character*30 file1,file2,file3
	dimension wvlt(10000)
c
	write(*,'(''Enter name of wvlt file : ''$)')
	read(*,'(a30)')file1
	write(*,'(''Enter name of deri file : ''$)')
	read(*,'(a30)')file2
	write(*,'(''Enter name of 2 deri file : ''$)')
	read(*,'(a30)')file3

	write(*,'(''Enter central frequency : ''$)')
	read(*,*)f1
	write(*,'(''Enter amplitude : ''$)')
	read(*,*)A
	write(*,'(''Enter number of samples out : ''$)')
	read(*,*)nsamp
	write(*,'(''Enter sampling interval in sec : ''$)')
	read(*,*)si
c
	pi=3.14159
	pi2=pi*pi
	alpha = pi2 * f1 * f1
	t0 = 8. / sqrt(2. * alpha)
c
	do 10 i=1,nsamp
	t = si * (float(i-1))
	wvlt(i) = A * exp(-alpha*(t0-t)**2)
 10	continue
c
	write(*,'(''First sample is : '',f10.6)')wvlt(1)
	write(*,'(''Last sample is : '',f10.6)')wvlt(nsamp)
	wvlt(1)=0.
	wvlt(nsamp)=0.
	open(unit=10,file=file1,form='unformatted',
     *  access='direct',status='unknown',recl=4*nsamp)
	write(10,rec=1)(wvlt(i),i=1,nsamp)
	close(unit=10)
c
	do 30 i=1,nsamp
	t = si *(float(i-1))
	wvlt(i)=A*(-2.)*alpha*(t-t0)*exp(-alpha*(t-t0)**2)
 30	continue
	write(*,'(''First sample is : '',f10.6)')wvlt(1)
	write(*,'(''Last sample is : '',f10.6)')wvlt(nsamp)
	wvlt(1)=0.
        wvlt(nsamp)=0.
        open(unit=10,file=file2,form='unformatted',
     *  access='direct',status='unknown',recl=4*nsamp)
        write(10,rec=1)(wvlt(i),i=1,nsamp)
        close(unit=10)
c
	do 50 i=1,nsamp
	t = si *(float(i-1))
	wvlt(i)=A*(1.-2.*alpha*(t-t0)*(t-t0))*exp(-alpha*
     *  (t0-t)*(t0-t))
 50	continue
	write(*,'(''First sample is : '',f15.6)')wvlt(1)
	write(*,'(''Last sample is : '',f15.6)')wvlt(nsamp)
	wvlt(1)=0.
        wvlt(nsamp)=0.
        open(unit=10,file=file3,form='unformatted',
     *  access='direct',status='unknown',recl=4*nsamp)
        write(10,rec=1)(wvlt(i),i=1,nsamp)
        close(unit=10)

	stop
	end
