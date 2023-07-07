! ============================ Module Description ===========================
! PHI_MOD: Initilizes Phi functions for scaler arguments. Contains:
!   1. phi - initilizes phi functions using contour integral & recursion
!            relation.
! ===========================================================================

module phi_mod

    ! Module Parameters
    use tools_mod, only: dp, norm, PI, II
    implicit none

    ! Cauchy Integral Settings
    real(dp),   parameter :: c_tol = 1.0_dp         ! Smallest Lambda for contour integral
    real(dp),   parameter :: c_R   = 2.0_dp * c_tol ! Contour Radius
    integer,    parameter :: c_M   = 32             ! Number of Points

    ! Scaling and Squaring Settings
    integer,    parameter :: t_M   = 20 ! Number of Taylor Series Terms to sum

    contains

    ! =======================================================================
    ! PHI   Evaluates \phi_i(L) for i=0,...,n and scaler/vector L using
    !       Recursion relation and Cauchy Integral Formula.
    !
    ! Arguments
    !
    !   L   (input) COMPLEX*16 array, dimensions(n)
    !       array of Lambda values cooresponding to PDE linear component
    !
    !   n   (input) INTEGER
    !       highest phi function to initialize.
    !
    !   P   (output) COMPLEX*16, dimensions(n,size(L))
    !       array where on exit P(i,j) = \phi_{i-1}(L(j))
    !
    !   algo_in (input) character, length=1
    !       Optional argument to specify algorithm
    !       "s" - scaling and squaring (default)
    !       "c" - contour integral
    ! =======================================================================

    subroutine phi(L, n, P, algo_in)
        
        ! Arguments
        complex(dp), intent(in)  :: L(:)
        integer,     intent(in)  :: n
        complex(dp), intent(out) :: P(n+1,size(L))
        character(len=1), optional, intent(in) :: algo_in
        
        ! Local Variables
        character(len=1) :: algo

        if(present(algo_in)) then
            algo = algo_in
        else
            algo = "s"
        endif

        if (algo == "c") then
            call phi_contour(L, N, P)
        else if(algo == "s") then
            call phi_scale_square(L, N, P)
        else 
            write (*,*) "phi_mod : invalid algorithm"
            call exit(1)
        endif

    end subroutine phi
    
    subroutine phi_scale_square(L, n, P)
        
        ! Arguments
        complex(dp), intent(in)  :: L(:)
        integer,     intent(in)  :: n
        complex(dp), intent(out) :: P(n+1,size(L))

        ! Local variables
        integer :: i

        do i = 1, size(L)
            call phi_scale_square_z(L(i), n, P(:, i))
        end do

    end subroutine phi_scale_square
    
    subroutine phi_scale_square_z(L, n, P)
        
        ! Arguments
        complex(dp), intent(in)  :: L
        integer,     intent(in)  :: n
        complex(dp), intent(out) :: P(n+1)
        
        ! Local variables
        integer :: s, k, i, j
        complex(dp), allocatable :: LS, PS(:)
        real(dp), allocatable :: factorial(:)
        
        ! Compute factorial
        allocate(factorial(t_M + n + 1))
        factorial(1) = 1.0_dp ! 0!
        do i = 1, t_M + n
            factorial(i + 1) = factorial(i) * real(i, dp)
        end do

        ! compute scaling factor
        s = max(0, ceiling(log(abs(L)) / log(2.0_dp)));
        LS = L / (2.0_dp ** real(s, dp));

        ! Initial Taylor Series Defenition
        do i = 0, n  
            ! Horner's Scheme
            P(i + 1) = LS / factorial(t_M + i + 1)
            do k = (t_M - 1), 1, -1
                P(i + 1) = LS / factorial(k + i + 1) + LS * P(i + 1);
            end do
            P(i + 1) = 1 / factorial(i + 1) + P(i + 1);
        end do

        ! Scale up by Powers of Two
        allocate(PS(n+1))
        do i = 1, s
            PS = P
            do k = 0, n
                if(k == 0) then
                    PS(1) = PS(1) ** 2.0_dp
                else
                    PS(k + 1) = ( 1.0_dp / 2.0_dp ** k) * P(1) * P(k + 1)
                    do j = 1, k
                         PS(k+1) = PS(k+1) + (1.0_dp / 2.0_dp ** k) * P(j+1) / factorial(k - j + 1)
                    end do
                end if
            end do
            P = PS
        end do
 
    end subroutine phi_scale_square_z

    subroutine phi_contour(L,n,P)
        ! Arguments
        complex(dp), intent(in)  :: L(:)
        integer,     intent(in)  :: n
        complex(dp), intent(out) :: P(n+1,size(L))
        ! Local Variables
        integer     :: i,j,k,nL
        real(dp)    :: f
        complex(dp) :: z(c_M)
        complex(dp) :: Li,Lzi,pp

        nL = size(L)
        P = 0;
        ! Set contour points
        do i=1,c_M
            z(i) = c_R * exp(2.0_dp*PI*II*(i-1.0_dp)/real(c_M,dp))
        enddo
        ! Compute Phi
        do i=1,nL
            Li = L(i)
            if(abs(Li) >= c_tol) then
                ! Direct Formula
                P(1,i) = exp(Li)
                f = 1.0_dp;
                do j=2,n+1
                    P(j,i) = (P(j-1,i) - 1.0_dp/f)/Li
                    f = f*(j-1)
                enddo
            else
                ! Cauchy Integral Formula
                do k=1,c_M
                    Lzi = Li + z(k)
                    pp = exp(Lzi)
                    P(1,i) = P(1,i) + pp/c_M
                    f = 1.0_dp;
                    do j=2,n+1
                        pp = (pp - 1.0_dp/f)/Lzi
                        P(j,i) = P(j,i) + pp/c_M;
                        f = f*(j-1)
                    enddo
                enddo
                ! remove imaginary roundoff if L(i) is real
                if(aimag(Li) == 0.0_dp) then
                  P(:,i) = REALPART(P(:,i))
                endif
            end if
        end do
    end subroutine phi_contour

    ! =========================================================================
    ! REMARK: the following function treats all elements in L with the same
    ! scaling and squaring parameters and can be used directly in place of 
    ! phi_scale_square. However, if L contains entries of largly varying 
    ! magnitude, then the phi functions for small L will be inaccurate
    ! since the entire vector will be scaled and squared as many times as 
    ! required for the largest entries of L.
    ! =========================================================================

    subroutine phi_scale_square_zvec(L, n, P)
        
        ! Arguments
        complex(dp), intent(in)  :: L(:)
        integer,     intent(in)  :: n
        complex(dp), intent(out) :: P(n+1,size(L))
        
        ! Local variables
        integer :: s, k, i, j
        complex(dp), allocatable :: LS(:), PS(:,:)
        real(dp), allocatable :: factorial(:)
        
        ! Compute factorial
        allocate(factorial(t_M + n + 1))
        factorial(1) = 1.0_dp ! 0!
        do i = 1, t_M + n
            factorial(i + 1) = factorial(i) * real(i, dp)
        end do

        ! compute scaling factor
        s = max(0, ceiling(log(maxval(abs(L))) / log(2.0_dp)));
        LS = L / (2.0_dp ** real(s, dp));

        ! Initial Taylor Series Defenition
        do i = 0, n  
            ! Horner's Scheme
            P(i + 1, :) = LS / factorial(t_M + i + 1)
            do k = (t_M - 1), 1, -1
                P(i + 1, :) = LS / factorial(k + i + 1) + LS * P(i + 1, :);
            end do
            P(i + 1, :) = 1 / factorial(i + 1) + P(i + 1, :);
        end do

        ! Scale up by Powers of Two
        allocate(PS(n+1, size(L)))
        do i = 1, s
            PS = P
            do k = 0, n
                if(k == 0) then
                    PS(1, :) = PS(1, :) ** 2.0_dp
                else
                    PS(k + 1, :) = ( 1.0_dp / 2.0_dp ** k) * P(1, :) * P(k + 1, :)
                    do j = 1, k
                         PS(k+1, :) = PS(k+1, :) + (1.0_dp / 2.0_dp ** k) * P(j+1, :) / factorial(k - j + 1)
                    end do
                end if
            end do
            P = PS
        end do
 
    end subroutine phi_scale_square_zvec

end module phi_mod