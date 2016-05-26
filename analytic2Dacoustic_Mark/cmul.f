* Copyright (c) Colorado School of Mines, 1994.
* All rights reserved.

	subroutine cmul(a,b,c,n)
c
	complex		a(n),b(n),c(n)
c
		do i = 1,n
		c(i) = a(i) * b(i)
		end do
c
	return
	end
c
c
