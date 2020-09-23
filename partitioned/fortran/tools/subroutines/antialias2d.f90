subroutine antialias2d(yh)
! performs 2d antialias 2/3 filter on vector yh
! zeros out modes (Nx/2 - b) ... (Nx/2) and (-Nx/2 + 1) ... (-Nx / 2 + b)
! === Parameters ===
! yh - 1d solution in fourier space. Wave
! === Output ===
! yh - anti aliased 1d solution in fourier space
    complex(dp), intent(inout)  :: yh(:,:)
    integer :: NpX, aX, bX, NpY, aY, bY
    
    NpX = size(yh, 1)
    NpY = size(yh, 2)
    
    aX  = NpX / 2
    bX  = NpX / 6

    aY  = NpY / 2
    bY  = NpY / 6
      
    yh( aX + 1 - bX : aX + 1 + bX, : ) = 0.0_dp;
    yh( : , aY + 1 - bY : aY + 1 + bY ) = 0.0_dp;

end subroutine antialias2d