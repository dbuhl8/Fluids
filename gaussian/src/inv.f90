

!This code is to check the validity of the psuedo inv code

program inv

    use gaussian_mod
    implicit none

    real(8) :: x(10), mat(10, 10), inv_mat(10, 10), s, kscale
    
    call random_number(x)
    kscale = sqrt(.1)
    mat = esqk(x, x, 10, 10, kscale)
    inv_mat = isqed(mat, 10)
    s = sum(matmul(inv_mat, mat))/10

    print *, "If this number is 1 then the psuedo inverse is a good approximation :", s

end program inv
