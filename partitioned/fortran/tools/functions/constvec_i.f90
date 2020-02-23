function constvec_i(n, val) result(v)
! Returns constant vector of val
! === Parameters ===
! n - vector dimenson
! val - (COMPLEX) 2x1 array containing inteval
! === Output ===
! c - (array) nx1 array with all entries set to val
integer, intent(in) :: n
complex(dp), intent(in) :: val
complex(dp) :: v(n)
integer :: j
do j=1,n
    v(j) = val
enddo
end function constvec_i