subroutine snapshot2d (flag, n, sn)
        implicit none
! arguments
	integer				:: flag, n, sn
! local variables
	character			:: nnum*5, name*256, snnum*3
! MPI variables
        integer                         :: recl_size
!**********************************************************************************************************
! executable code:

        write (nnum,'(i5.5)') n
        write (snnum,'(i2.2)') sn

!**********************************************************************************************************
! Forward modeling output
     do c=1,snanum
       snap = snasta + ((c-1)*snainc)
       if (snap == n) then


    recl_size = prec * NX * NZ

    name = './v_2D'//'.snap'//nnum//'.shot'//snnum

    open(unit=50,file=name,form='unformatted',access='direct',recl=recl_size)
    write(50,rec=1) v(:,:)
    close(50,status='keep')

  endif

enddo


end subroutine snapshot2d
