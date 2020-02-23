! ============================ Module Description ===========================
! etdmp_mod Implementation of ETDMP PESE method. Contains a 2 subroutines:
!   1. etdmp    - ETDMP timestepping code
!   2. initW    - initializes ETD coefficients
! ===========================================================================

module etdab_mod

    use tools_mod,  only: dp, tic, toc, weights, linspace, constvec_r
    use phi_mod,    only: phi

    implicit none
    ! Custom Settings Data Type
    type etdab_settings
        integer               :: n
    end type etdab_settings

    contains

    ! =======================================================================
    ! SOLVE   Implements ETDPBM Numerical Time Integrator with PMFO_mS AII
    !
    ! Arguments
    !
    !   L       (input) COMPLEX*16 array (square)
    !           vector cooresponding to PDE linear operator
    !
    !   N       (input) SUBROUTINE:
    !           cooresponds to PDE nonlinear operator. Must be of the form:
    !               subroutine N(t,y_in,N_out)
    !                   real(kind=8),    intent(in)  :: t
    !                   complex(kind=8), intent(in)  :: y_in(:)
    !                   complex(kind=8), intent(out) :: N_out(:)
    !               end subroutine N
    !
    !   y_in    (input) COMPLEX*16 array, dimension (N)
    !           Initial condition.
    !
    !   tspan   (input) DOUBLE array, dimension(2)
    !           contains left and right integration bounds
    !
    !   options (input) ETDMP_SETTINGS
    !           derived data type that contains one field:
    !               options%z : DOUBLE(:) (nodes)
    !               options%b : DOUBLE(:) (expansion points)
    !               options%alpha : DOUBLE (extrapolation parameter)
    !               options%mS : DOUBLE(:) setminus
    !               options%m : DOUBLE number of correction iterations
    !
    !   y_out   (output) COMPLEX*16 array, dimension (N)
    !           Initial condition.
    !
    !   times   (output) DOUBLE array, dimension(4)
    !           times(1:2) = [cputime,clocktime] for timestepping loop
    !           times(3:4) = [cputime,clocktime] for initializing coefficients,
    ! =======================================================================

    subroutine etdab(L, N, tspan, y_in, Nt, settings, y_out, times)

        ! Arguments
        complex(dp),           intent(in)  :: L(:)
        real(dp),              intent(in)  :: tspan(2)
        complex(dp),           intent(in)  :: y_in(:)
        integer,               intent(in)  :: Nt
        type(etdab_settings), intent(in)  :: settings
        complex(dp),           intent(out) :: y_out(:)
        real(dp),              intent(out) :: times(4)
        ! Define interface subroutine function N
        interface
            subroutine N(t,yh_in, N_out, thread_id)
            import :: dp
            real(dp),    intent(in)  :: t
            complex(dp), intent(in)  :: yh_in(:)
            complex(dp), intent(out) :: N_out(:)
            integer, intent(in), optional :: thread_id
            end subroutine N
        end interface

        ! Local Variables
        integer  :: i, j, k, m, q, nL, LN_AII_len, lAII
        real(dp) :: r, h, t
        complex(dp), allocatable :: W_m(:,:), WA_m(:,:,:), W_i(:,:,:)
        complex(dp), allocatable :: E_m(:), EA_m(:,:), E_i(:,:)
        complex(dp), allocatable :: yA(:,:), y_n(:), N_n(:,:)
        real(dp),    allocatable :: z(:)

        ! Initialize Method Parameters
        q  = settings%n
        z  = linspace(0.0_dp, real(q-1,dp), q)
        h  = (tspan(2) - tspan(1)) / (real(Nt + q - 1, dp))
        nL = size(L)

        allocate(E_m(nL), W_m(nl, q), EA_m(nL, 1), WA_m(nL, q, 1))
        allocate(E_i(nL, q), W_i(nL, q, q))
        allocate(y_n(nL), yA(nL, q), N_n(nL, q))

        ! === Initial ETD Coefficients =======================================================================================================================================
        call tic()
        call initWMP(L, h, z, constvec_r(1, real(q - 1, dp)), constvec_r(1, real(q, dp)), WA_m, EA_m)
        call initWMP(L, h, z, constvec_r(q, 0.0_dp), z,  W_i, E_i)
        E_m(:) = EA_m(:, 1)  ! eliminate extra dimension
        W_m(:,:) = WA_m(:,:,1) ! eliminate extra dimension
        call toc(times(3:4))
        deallocate(EA_m)
        deallocate(WA_m)

        ! === Initial Conditions : Iterator ==================================================================================================================================
        t = tspan(1)
        do i = 1, q
            yA(:, i) = y_in;
        enddo
        do i = 1, q
            do j = 1, q
                call N(t + h * z(j), yA(:, j), N_n(:, j))
            enddo
            do j = 2, q
                yA(:,j) = E_i(:, j) * yA(:, 1) + sum(W_i(:, :, j) * N_n, 2)
            enddo
        enddo

        ! === Time-Stepping Loop =============================================================================================================================================
        call tic()
        do j = 1, q
            call N(t + h * z(j), yA(:, j), N_n(:, j))
        enddo
        y_n = yA(:, q)
        t   = h * q
        do i = 1, Nt
            y_n = E_m * y_n + sum(W_m * N_n, 2)
            N_n(:, 1 : q - 1) = N_n(:, 2 : q)
            call N(t + h, y_n, N_n(:, q))
            t = t + h
        enddo
        call toc(times(1:2))
        y_out = y_n;

    end subroutine etdab

    ! =======================================================================
    ! INITW   Initializes ETDMP W functions for matrix L using weights
    !         function by Fornberg.
    !
    ! Arguments
    !
    !   L   (input) COMPLEX*16 array, dimension(M)
    !       vector cooresponding to PDE linear operator
    !
    !   h   (input) DOUBLE
    !       timestep
    !
    !   tau (input) DOUBLE array, dimensions(n,1)
    !       normalized quadrature points
    !
    !   a   (input) DOUBLE array, dimension(n,1)
    !       left integration bounds
    !
    !   b   (input) DOUBLE array, dimension(n,1)
    !       right integration bounds
    !
    !   W   (output) COMPLEX*16 array, dimension(size(tau), size(tau), size(L))
    !       contains the ETD Coefficients so that W(:,:,j) is the ETD integration
    !       matrix cooresponding to lambda=L(j)
    !
    !   E   (output) COMPLEX*16 array, dimension(size(tau),size(L))
    !       contains Matrix Exponentials where E(i,:) = exp(h_i L)
    ! =======================================================================

    subroutine initWMP(L, h, z, a, b, W, E)

        ! Argument
        complex(dp), intent(in)  :: L(:)
        real(dp),    intent(in)  :: h
        real(dp),    intent(in)  :: z(:), a(:), b(:)
        complex(dp), intent(out) :: W(size(L), size(z), size(a))
        complex(dp), intent(out) :: E(size(L), size(a))

        ! Declare variables
        real(dp) :: eta(size(a))
        real(dp) :: sz(size(z))
        real(dp) :: FD(size(z), size(z))
        complex(dp), allocatable :: p(:, :)
        integer :: i, n, q

        q     = size(z)
        n     = size(a)
        eta   = b - a
        allocate(p(q + 1, size(L)))

        ! Form W Functions
        do i = 1, n
            if(eta(i) .ne. 0.0_dp) then
                sz = (z - a(i)) / eta(i)                                        ! scaled nodes
                call weights(0.0_dp, sz, q - 1, FD)                             ! finite difference matrix
                call phi(eta(i) * h * L, q, p)                                  ! phi functions 0 to n
                W(:, :, i) = transpose(h * eta(i) * matmul(FD, p(2:q + 1,:)))   ! store ith row of W matrix
                E(:, i)    = p(1, :)                                            ! store exp(h_i L)
            else
                W(:, :, i) = 0.0_dp
                E(:, i)    = 1.0_dp
            endif
        enddo
        deallocate(p)

    end subroutine initWMP

end module etdab_mod
