program test
integer :: i
do i=1, 100
    open(11, file='Output'//trim(str(i))//'.txt')
    write (11, *) i
    close (11)
end do

contains 
    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
            write (str, *) k
            str = adjustl(str)
    end function str

end program test

