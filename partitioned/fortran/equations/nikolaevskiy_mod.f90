! ============================ Module Description ===========================
! nikolaevskiy_mod stores experiment parameters and subroutines for evaluating 
! Linear and Nonlinear Operators for the Nikolaevskiy equation when solving
! in Fourier space (Diagonal Lambda).
!
! Contains 5 subroutines:
!   1. init         - initializes spatial grid, fft plans, experiment 
!                     parameters.
!   2. ic           - sets initial condition in physical space
!   3. L            - returns diagonal entries of linear operator Lambda
!   4. N            - subroutine for evaluating nonlinear operator
!   5. error_filter - defines error function used for numerical experiments
! ===========================================================================

module nikolaevskiy_mod

    ! Module Parameters
    use tools_mod, only: dp, PI, II, logspace, plinspace, fourierWaveNum, relerror_c, antialias
    implicit none
    ! Numerical Parameters
    integer,    parameter :: Np                   = 2**12                         ! Number of spatial points (Must be Even)
    real(dp),   parameter :: tspan(2)             = [0.0_dp, 50.0_dp]             ! Time integration window
    logical,    parameter :: reference_methods(3) = [.true., .true., .false.]     ! Methods for Reference Solution (ETDSDC,IMEXSDC,ETDRK)
    integer,    parameter :: num_tests            = 16                            ! Number of Numerical Tests
    real,       parameter :: smallest_F           = 1e3_dp                        ! Smallest Number of Function Evaluations
    real,       parameter :: largest_F            = 1e5_dp                        ! Maximum Number of Function Evaluations
    logical,    parameter :: antialiasing_enabled = .true.                        ! determines if 2/3 antialiasing rule should be applied.
    ! Equation parameters
    real(dp), parameter     :: Lx    = 150.0_dp * PI          ! Spatial Domain Size
    real(dp), parameter     :: r     = 0.25_dp
    real(dp), parameter     :: alpha = 2.1_dp
    real(dp), parameter     :: beta  = 0.77_dp
    character(len=*), parameter :: eqn_name = "nikolaevskiy"    ! Used to save results to correct folder
    ! Storage Arrays
    complex(dp), dimension(:), allocatable :: y0,ks
    complex(dp), dimension(:,:), allocatable :: y,yh
    real(dp),    dimension(:), allocatable :: xs
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
        integer :: i, nt
        if(present(max_num_threads)) then
            nt = max_num_threads
        else
            nt = 1
        endif
        ! Init FFT plans
        allocate(y(Np, nt), yh(Np, nt), plan_backward(nt), plan_forward(nt))
        do i = 1, nt
            call dfftw_plan_dft_1d_ (plan_forward(i), Np, y(:, i), yh(:, i), FFTW_FORWARD, FFTW_MEASURE)
            call dfftw_plan_dft_1d_ (plan_backward(i), Np, yh(:, i), y(:, i), FFTW_BACKWARD, FFTW_MEASURE)
        end do
        ! allocate arrays
        allocate(ks(Np), xs(Np))
        ! set x domain
        xs = plinspace(-Lx/2.d0,Lx/2.d0,Np)
        ! set initial condition
        call ic()
        ! Set Fourier wave numbers
        ks = fourierWaveNum(Np,Lx)
        ! Set Function Evalutations for tests
        Fs = logspace(log(smallest_f)/log(10.0_dp),log(largest_F)/log(10.0_dp),num_tests)
    end subroutine

    ! Initial Condition
    subroutine ic()
        y(:,1) = sin(xs) + 0.1_dp * sin(xs/25.0_dp) ! NOTE: xs global var set in init
        call dfftw_execute(plan_forward)
        y0 = yh(:,1)
        if(antialiasing_enabled) then
            call antialias(y0)
        endif
    end subroutine ic

    ! Linear Operator
    subroutine L(lambda)
        complex(dp), dimension(Np), intent(out) :: lambda
        lambda = (r - 1.0_dp)*(ks**2) - (alpha*II*ks**3) + (2*ks**4) + beta*(II*ks**5) - (ks**6)
        if(antialiasing_enabled) then
            call antialias(lambda)
        endif
    end subroutine L

    ! Nonlinear Operator
    subroutine N(t,yh_in, N_out, thread_id)
        complex(dp),    intent(in) :: t
        complex(dp), dimension(:), intent(in)  :: yh_in
        complex(dp), dimension(:), intent(out) :: N_out
        integer, intent(in), optional :: thread_id
        integer :: tid
        if(present(thread_id)) then
            tid = thread_id
        else
            tid = 1
        endif

        yh(:,tid) = yh_in
        call dfftw_execute(plan_backward(tid))
        y(:,tid) = (y(:,tid)/Np)**2;
        call dfftw_execute(plan_forward(tid))
        N_out = -0.5_dp * II * ks * yh(:,tid);
        if(antialiasing_enabled) then
            call antialias(N_out)
        endif
    end subroutine N

    ! Error Filter (Inf norm in physical space)
    function error_filter(estimate,exact) result(error)

        complex(dp), intent(in)  :: estimate(Np)
        complex(dp), intent(in)  :: exact(Np)

        real(dp)  :: error
        complex(dp), allocatable :: a(:),b(:)

        allocate(a(Np), b(Np))

        yh(:, 1) = estimate / Np
        call dfftw_execute(plan_backward(1))
        a = y(:, 1)

        yh(:, 1) = exact / Np
        call dfftw_execute(plan_backward(1))
        b = y(:, 1)

        error = relerror_c(a,b)
    end function error_filter

end module nikolaevskiy_mod