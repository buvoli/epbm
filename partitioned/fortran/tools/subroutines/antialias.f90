subroutine antialias(yh)
! performs 1d antialias 2/3 filter on vector yh
! zeros out modes (Nx/2 - b) ... (Nx/2) and (-Nx/2 + 1) ... (-Nx / 2 + b) where b = Nx/6
! === Parameters ===
! yh - 1d solution in fourier space. Expected wave number ordering is [0:params.Nx/2, -params.Nx/2+1:-1]
! === Output ===
! yh - anti aliased 1d solution in fourier space
    complex(dp), intent(inout)  :: yh(:)
    integer :: a, b, Np
    
    Np = size(yh)
    a  = Np / 2
    b  = Np / 6 ! one-half filter length
    yh( a + 1 - b : a + 1 + b) = 0.0_dp;
end subroutine antialias