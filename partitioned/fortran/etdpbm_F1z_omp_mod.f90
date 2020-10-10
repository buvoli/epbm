! ============================ Module Description ==============================
! etdpbm_omp - OpenMP Implementation of ETD PBM method with:
!               Nodes: Left Sweeping real
!               Active Index Set (AIS): PMFO
!               L_n AIS: PMFO_mS
!               L_y AIS: PMFO
!               Endpoints: Fixed, n=1
!               Parallel Iterator (PMFO)
!
!   Note: this code will not produce correct results if it is compiled without
!         OPENMP
!
! Contains subroutines:
!   1. etdpbm_omp - ETDPBM timestepping code
!   2. initW      - initializes requesite ETD coefficients
!   3. getLN_AII  - forms L_n AII set
! ==============================================================================

module etdpbm_F1z_omp_mod

    use omp_lib
    use tools_mod,  only: dp, tic, toc, weights_z, constvec_z
    use phi_mod,    only: phi

    implicit none
    ! Custom Settings Data Type
    type etdpbm_F1z_omp_settings
        complex(dp), allocatable :: z(:)
        complex(dp), allocatable :: b(:)
        real(dp)              :: alpha
        integer               :: m
        integer,  allocatable :: mS(:)
    end type etdpbm_F1z_omp_settings

    contains

    ! =======================================================================
    ! ETDPBM_OMP   Implements EPBM Numerical Time Integrator
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

    subroutine etdpbm_F1z_omp(L, N, tspan, y_in, Nt, settings, y_out, times)

        ! Arguments
        complex(dp),           intent(in)  :: L(:)
        real(dp),              intent(in)  :: tspan(2)
        complex(dp),           intent(in)  :: y_in(:)
        integer,               intent(in)  :: Nt
        type(etdpbm_F1z_omp_settings), intent(in)  :: settings
        complex(dp),           intent(out) :: y_out(:)
        real(dp),              intent(out) :: times(4)
        ! Define interface subroutine function N
        interface
            subroutine N(t,yh_in, N_out, thread_id)
            import :: dp
            complex(dp),    intent(in)  :: t
            complex(dp), intent(in)  :: yh_in(:)
            complex(dp), intent(out) :: N_out(:)
            integer, intent(in), optional :: thread_id
            end subroutine N
        end interface

        ! Local Variables
        integer  :: i, j, k, m, q, nL, LN_AII_len, thread_id, ind
        real(dp) :: r, h
        complex(dp) :: t
        complex(dp), allocatable :: W_m(:,:,:), W_i(:,:,:)
        complex(dp), allocatable :: E_m(:,:),   E_i(:,:)
        complex(dp), allocatable :: y_n(:),     y_np1(:, :), N_n(:,:)
        complex(dp), allocatable :: zh(:), bh(:), zh_LN(:)
        integer,     allocatable :: LN_AII(:)

        ! verify valid nodes
        if(settings%z(1) .ne. -1.0_dp) then
            write (*,*) "Invalid Nodes provided to ETDPBM_OMP; z(1) must be equal to -1."
            call exit(1)
        endif

        ! verify valid enpoints (This implementation only supports Fixed (\ell = 1)
        if(size(settings%b .eq. -1.0_dp) .ne. size(settings%z)) then
            write (*,*) "Invalid Expansion provided to ETDPBM_OMP; b(j) must be equal to -1."
            call exit(1)
        endif

        ! verify valid sM (This implementation only supports Fixed (sM = 1)
        if((size(settings%mS) .ne. 1) .or. (settings%mS(1) .ne. 1)) then
            write (*,*) "Invalid Set SM provided to ETDPBM_OMP; Sm must be equal to [1]."
            call exit(1)
        endif

        ! Initialize Method Parameters
        q = size(settings%z)
        m = settings%m
        ! scale node vector z and endpoint vector b relative to the stepsize h
        zh = settings%z / settings%alpha
        bh = settings%b / settings%alpha

        call getLN_AII(q, settings%mS, LN_AII)
        LN_AII_len = size(LN_AII)
        zh_LN = (settings%z(LN_AII) + 1.0_dp) / settings%alpha ! scaled z used for nonlinear function evaluations. shift by 1 so that we do not need to integrate backwards in time

        nL = size(L)
        h  = (tspan(2) - tspan(1)) / (real(Nt, dp))

        allocate(E_m(nL, q), W_m(nL, LN_AII_len, q))
        allocate(E_i(nL, q), W_i(nL, LN_AII_len, q))
        allocate(y_n(nL),    y_np1(nL, q), N_n(nL, LN_AII_len))

        ! === Initial ETD Coefficients =======================================================================================================================================
        call tic()
        call initWMP(L, h, zh(LN_AII), bh, zh + 1.0_dp, W_m, E_m)
        call initWMP(L, h, zh(LN_AII), constvec_z(q, zh(1)), zh,  W_i, E_i)
        call toc(times(3:4))

        ! === Initial Conditions : Iterator ==================================================================================================================================
        t   = tspan(1)
        y_n = y_in;
        call N(t + h * zh_LN(1), y_n, N_n(:, 1))
        do i = 2, LN_AII_len
            N_n(:, i) = N_n(:,1);
        enddo

        do i = 1, q
            do j = 2, q
                y_np1(:,j) = E_i(:, j) * y_n + sum(W_i(:,:,j) * N_n, 2)
            enddo
            do j = 1, LN_AII_len
                call N(t + h * zh_LN(j), y_np1(:, LN_AII(j)), N_n(:, j))
            enddo
        enddo

        ! === Time-Stepping Loop =============================================================================================================================================
        call tic()
        thread_id = 1

        !$ call omp_set_num_threads(q)
        !$omp parallel private(thread_id, t, ind)
        !$ thread_id = omp_get_thread_num() + 1
        do i = 1, Nt

            ! --- propagator -------------------------------------------------------------------------------------------------------------------------------------------------
            !$omp barrier
            y_np1(:,thread_id) = E_m(:, thread_id) * y_n + sum(W_m(:,:,thread_id) * N_n, 2)
            !$omp barrier
            if(thread_id .le. LN_AII_len) then
                call N(t + h * zh_LN(thread_id), y_np1(:,LN_AII(thread_id)), N_n(:, thread_id), thread_id)
                t   = t + h
            else
                y_n = y_np1(:,1)
                t   = t + h
            end if

             ! --- iterator ---------------------------------------------------------------------------------------------------------------------------------------------------
             do k = 1, m
                !$omp barrier
                if(thread_id .le. LN_AII_len) then
                    ind = LN_AII(thread_id);
                    y_np1(:, ind) = E_i(:, ind) * y_n + sum(W_i(:,:, ind) * N_n, 2)
                endif
                !$omp barrier
                if(thread_id .le. LN_AII_len) then
                    call N(t + h * zh_LN(thread_id), y_np1(:, ind), N_n(:, thread_id), thread_id)
                endif
             enddo

        enddo
        !$omp end parallel

        call toc(times(1:2))
        y_out = y_n;

    end subroutine etdpbm_F1z_omp

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
        complex(dp), intent(in)  :: z(:), a(:), b(:)
        complex(dp), intent(out) :: W(size(L), size(z), size(a))
        complex(dp), intent(out) :: E(size(L), size(a))

        ! Declare variables
        complex(dp) :: eta(size(a))
        complex(dp) :: sz(size(z))
        complex(dp) :: FD(size(z), size(z))
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
                call weights_z((0.0_dp, 0.0_dp), sz, q - 1, FD)                           ! finite difference matrix
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


    ! =======================================================================
    ! GETLN_AII   Creates array with active input indices for polynomial L_N.
    !             This function simply computes setdiff(1 : q, sM)
    !
    ! Arguments
    !
    !   q   (input) INTEGER
    !       number of nodes
    !
    !   sM  (input) INTEGER array
    !       array indices to remove
    !
    !   AII (input) INTEGER array
    !       active input indices for L_n
    ! =======================================================================

    subroutine getLN_AII(q, sM, AII)
        ! Arguments
        integer, intent(in)  :: q
        integer, intent(in)  :: sM(:)
        integer, allocatable, intent(inout) :: AII(:)

        ! Variables
        integer :: i
        integer :: rcount, count
        logical, allocatable :: AII_flag(:)

        ! Search valid Indices (AII_flag(i) = true \implies i is valid)
        allocate(AII_flag(q))
        AII_flag = .true.
        rcount   = q
        do i = 1, size(sM)
            if((sM(i) .le. q) .and. (SM(i) .ge. 1)) then
                AII_flag(sM(i)) = .false.
                rcount = rcount - 1
            end if
        end do

        ! populate array with valid indices
        allocate(AII(rcount))
        count = 1
        do i = 1, q
            if(AII_flag(i) .eqv. .true.) then
                AII(count) = i
                count = count + 1
            end if
        end do

        deallocate(AII_flag)

    end subroutine

end module etdpbm_F1z_omp_mod
