! ============================ Module Description ===========================
! quasigeostrophic_mod stores experiment parameters and subroutines for
! evaluating Linear and Nonlinear Operators for the quasigeostrophic equation
! when solving in Fourier space (Diagonal Lambda).
!
! Contains 5 subroutines:
!   1. init         - initializes spatial grid, fft plans, experiment 
!                     parameters.
!   2. ic           - sets initial condition in physical space
!   3. L            - returns diagonal entries of linear operator Lambda
!   4. N            - subroutine for evaluating nonlinear operator
!   5. buildDiffMat - build relevant differentiation matrices
!   6. error_filter - defines error function used for numerical experiments
! ===========================================================================

module quasigeostrophic_mod

    ! Module Parameters
    use tools_mod, only: dp, PI, II, logspace, plinspace, fourierwavenum, meshgrid, relerror_c, antialias2d
    implicit none
    ! Numerical Parameters
    integer,    parameter :: Nx                   = 2**8                          ! Number of x spatial points (Must be Even)
    integer,    parameter :: Ny                   = 2**8                          ! Number of y spatial points (Must be Even)
    integer,    parameter :: Np                   = Nx * Ny                       ! Total Number of spatial points
    real(dp),   parameter :: tspan(2)             = [ 0.0_dp, 5.0_dp ]            ! Time integration window
    logical,    parameter :: reference_methods(3) = [ .true., .true., .false. ]   ! Methods for Reference Solution (ETDSDC,IMEXSDC,ETDRK)
    integer,    parameter :: num_tests            = 16                            ! Number of Numerical Tests
    real,       parameter :: smallest_F           = 1.5e3_dp                      ! Smallest Number of Function Evaluations
    real,       parameter :: largest_F            = 1.0e5_dp                      ! Maximum Number of Function Evaluations
    logical,    parameter :: antialiasing_enabled = .true.                        ! determines if 2/3 antialiasing rule should be applied.
    ! Storage Arrays
    real(dp), parameter :: Lx       = 2.0_dp * PI
    real(dp), parameter :: Ly       = 2.0_dp * PI
    real(dp), parameter :: epsilon  = 1.0e-2_dp
    real(dp), parameter :: v        = 1.0e-14_dp
    real(dp), parameter :: beta     = 10.0_dp
    character(len=*), parameter :: eqn_name = "quasigeostrophic"                   ! Used to save results to correct folder
    ! Storage Arrays
    complex(dp), dimension(:), allocatable   :: y0
    real(dp),    dimension(:,:), allocatable :: xs,ys
    complex(dp), dimension(:,:,:), allocatable :: P_H, P_LAP, N_temp
    complex(dp), dimension(:, :, :), allocatable :: y, yh
    ! Diff Matrices
    complex(dp), dimension(:,:), allocatable :: DX, DY, LAP, ILAP
    ! FFT Parameters
    integer(kind=8), allocatable :: plan_backward(:), plan_forward(:)
    include "fftw3.f90"
    ! Additional Settings
    real(dp), allocatable :: Fs(:) ! Function Counts to test (set in init function)

    contains

    subroutine init(max_num_threads)
        ! arguments
        integer, intent(in), optional :: max_num_threads
        ! variables
        integer     :: i, nt
        real(dp)    :: k,s
        real(dp)    :: h = Lx/Nx
        real(dp), allocatable :: xsr(:), ysr(:)
        if(present(max_num_threads)) then
            nt = max_num_threads
        else
            nt = 1
        endif
        ! Init FFT plans
        allocate(y(Nx, Ny, nt), yh(Nx, Ny, nt), plan_backward(nt), plan_forward(nt))
        do i = 1, nt
            call dfftw_plan_dft_2d_ (plan_forward(i), Nx, Ny, y(:, :, i), yh(:, :, i), FFTW_FORWARD, FFTW_MEASURE)
            call dfftw_plan_dft_2d_ (plan_backward(i), Nx, Ny, yh(:, :, i), y(:, :, i), FFTW_BACKWARD, FFTW_MEASURE)
        end do
        ! allocate arrays
        allocate(xs(Nx,Ny), ys(Nx,Ny), P_H(Nx,Ny,nt), P_LAP(Nx,Ny,nt), N_temp(Nx,Ny,nt))
        allocate(xsr(Nx),ysr(Ny),y0(Np))
        ! setup spatial domain
        xsr = plinspace(-Lx/2.0_dp,Lx/2.0_dp,Nx)
        ysr = plinspace(-Ly/2.0_dp,Ly/2.0_dp,Ny)
        call meshgrid(xsr,ysr,xs,ys)
        ! Set Differentiation Matrices
        call buildDiffMat()
        ! Set Initial condition
        call ic()
        ! Set Function Evalutations for tests
        Fs = logspace(log(smallest_f)/log(10.0_dp),log(largest_F)/log(10.0_dp),num_tests)
    end subroutine

    ! Form Initial Condition
    subroutine ic()
        y(:,:,1) = (1.0_dp/8.0_dp) * exp(-8.0_dp*(2.0_dp*ys**2.0_dp + 0.5_dp*xs**2.0_dp - PI/4.0_dp)**2.0_dp)
        call dfftw_execute(plan_forward(1))
        if(antialiasing_enabled) then
            call antialias2d(yh(:,:,1))
        endif
        y0 = reshape(LAP * yh(:,:,1), [ Nx * Ny ] )
    end subroutine ic

    subroutine L(lambda)
        complex(dp), dimension(Np), intent(out) :: lambda
        complex(dp), allocatable :: l2d(:,:)

        allocate(l2d(Nx, Ny))
        l2d = -1.0_dp*beta*DX*ILAP - epsilon - v*(DX**8 + DY**8)
        if(antialiasing_enabled) then
            call antialias2d(l2d)
        endif

        lambda = reshape(l2d, [ Np  ])
    end subroutine L

    ! Nonlinear Operator
    subroutine N(t,LAP_P_H_in, N_out, thread_id)
        complex(dp),    intent(in) :: t
        complex(dp), dimension(:), intent(in)  :: LAP_P_H_in !Fourier transform of Laplacian if Phi
        complex(dp), dimension(:), intent(out) :: N_out
        integer, intent(in), optional :: thread_id
        integer :: i, j, tid
        if(present(thread_id)) then
            tid = thread_id
        else
            tid = 1
        endif
        yh(:,:,tid) = reshape(LAP_P_H_in, [Nx, Ny])
        ! Fourier Transform of Phi
        P_H(:,:,tid) = ILAP * yh(:,:,tid)
        ! Laplacian of Phi
        call dfftw_execute(plan_backward(tid))
        P_LAP(:,:,tid) = (1.0_dp/Np) * y(:,:,tid)
        ! Compute N1
        yh(:,:,tid) = -Dy * P_H(:,:,tid)
        call dfftw_execute(plan_backward(tid))
        y(:,:,tid)  = (1.0_dp/Np) * y(:,:,tid) * P_LAP(:,:,tid)
        call dfftw_execute(plan_forward(tid))
        N_temp(:,:,tid)  = -DX * yh(:,:,tid)
        ! Compute N2
        yh(:,:,tid) = Dx * P_H(:,:,tid)
        call dfftw_execute(plan_backward(tid))
        y(:,:,tid)  = (1.0_dp/Np) * y(:,:,tid) * P_LAP(:,:,tid)
        call dfftw_execute(plan_forward(tid))
        N_temp(:,:,tid) = N_temp(:,:,tid) - (DY * yh(:,:,tid))
        if(antialiasing_enabled) then
            call antialias2d(N_temp(:,:,tid))
        endif
        N_out  = reshape(N_temp(:,:,tid),[Np])
    end subroutine N

    ! Differention Matrices
    subroutine buildDiffMat()
        real(dp), allocatable :: k_x(:), k_y(:)
        integer :: i
        allocate(DX(NX,NY), DY(NX,NY), LAP(NX,NY), ILAP(NX,NY))
        ! Get Fourier wave numbers
        k_x = fourierWaveNum(Nx,Lx)
        k_y = fourierWaveNum(Ny,Ly)
        ! Form Discrete D/DX Diff matrix
        do i=1,nx
            DX(:,i) = II * k_x(i)
        enddo
        ! Form Discrete D/DY Diff matrix
        do i=1,nx
            DY(i,:) = II * k_y(i)
        enddo
        LAP  = DX**2 + DY**2        ! Laplacian operator
        ILAP = 1.0_dp/LAP           ! Inverse Laplacian
        ILAP(1,1) = (0.0_dp,0.0_dp) ! Set (0,0) mode to zero
    end subroutine buildDiffMat


    ! Error Filter (Inf norm in physical space)
    function error_filter(estimate,exact) result(error)

        complex(dp), intent(in)  :: estimate(Np)
        complex(dp), intent(in)  :: exact(Np)

        real(dp)  :: error
        complex(dp), allocatable :: a(:),b(:)

        allocate(a(Np),b(Np))

        ! IFFT of estimate
        yh(:,:,1) = (1.0_dp / Np) * reshape(estimate, [Nx, Ny])
        call dfftw_execute(plan_backward(1))
        a(:) = reshape(y(:,:,1),[Np])

        ! IFFT of exact
        yh(:,:,1) = (1.0_dp / Np) * reshape(exact, [Nx, Ny])
        call dfftw_execute(plan_backward(1))
        b(:) = reshape(y(:,:,1),[Np])

        error = relerror_c(a,b)
    end function error_filter

end module quasigeostrophic_mod
