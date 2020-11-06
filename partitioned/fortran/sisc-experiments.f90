! ============================ Program Description ==========================
! EXPERIMENTRUN_MP  Generates results for:
! 	Buvoli T, "Exponential Polynomial Time Integrators"
! ===========================================================================

program SiscExperiment

    ! ===================== Specify Equation Module Here ====================
    use kuramoto_mod, only: L,N,init,y0,Np,tspan,Fs,reference_methods,eqn_name,error_filter ! quasigeostrophic, nikolaevskiy, kuramoto, kdv, nls
    ! =======================================================================

    use tools_mod,   only: dp, linspace, chebpts, legpts, print_cvector, isfinite, relerror_c, &
        write_rmatrix, write_rvector, write_cvector, constvec_r
    use etdrk4_mod,  only: etdrk4
    use etdsdc_mod,  only: etdsdc,  etdsdc_settings
    use imexsdc_mod, only: imexsdc, imexsdc_settings
    use etdab_mod,   only: etdab,   etdab_settings
    use etdpbm_mod,  only: etdpbm,  etdpbm_settings
    !use etdpbm_omp_mod,   only: etdpbm_omp,   etdpbm_omp_settings
    use etdpbm_F1_omp_mod,  only: etdpbm_F1_omp, etdpbm_F1_omp_settings
    use etdpbm_F1z_omp_mod,  only: etdpbm_F1z_omp, etdpbm_F1z_omp_settings
    use etdpbm_F1_mod,  only: etdpbm_F1, etdpbm_F1_settings

    implicit none

    ! Local Variables
    complex(dp),allocatable     :: y_reference(:), y_out(:),Lambda(:)
    integer                     :: F,relcost,order,Nt,solution_count,i,j,n_imex, n_etd, n_etd_mp, n_Fs, sdc_reference_n, num_methods
    integer                     :: nn, q, shift
    TYPE(etdsdc_settings)       :: etdsdc_s
    TYPE(imexsdc_settings)      :: imexsdc_s
    TYPE(etdab_settings)        :: etdab_s
    TYPE(etdpbm_settings)       :: etdpbm_s
    !TYPE(etdpbm_omp_settings)   :: etdpbm_omp_s
    TYPE(etdpbm_F1_omp_settings):: etdpbm_F1_omp_s
    TYPE(etdpbm_F1z_omp_settings) :: etdpbm_F1z_omp_s
    TYPE(etdpbm_F1_settings)    :: etdpbm_F1_s
    integer, allocatable        :: ETD_SDC_orders(:), ETD_AB_orders(:), ETD_PBM_qs(:)
    real(dp)                    :: time(4)
    real(dp),allocatable        :: times(:,:), errors(:,:), Nts(:,:), hs(:,:)
    character(len=*), parameter :: results_dir = "../results/"
    LOGICAL                     :: omp_active = .FALSE.

    ! check for omp
    !$ omp_active = .TRUE.

    ! ETDSDC Methods to test
    ETD_SDC_orders  = [2, 4, 6, 8]
    ETD_AB_orders   = [2, 4, 6, 8]
    ETD_PBM_qs      = [2, 4, 6, 8]

    ! Init Equation
    allocate(Lambda(Np))
    call init(8)
    call L(Lambda)

    write (*,"(a,a,a)") "Running experiment for ", eqn_name, " equation"
    write (*,"(a,$)") "Computing reference solution... "

    ! === Compute Reference Solution ===
    allocate(y_out(Np),y_reference(Np))
    F = 4*Fs(size(Fs))
    y_reference = (0.d0,0.d0)
    solution_count = 0
    sdc_reference_n = 32

    if(reference_methods(1)) then ! Use ETDSDC_N^N-1
        etdsdc_s%tau = chebpts(sdc_reference_n,(/0.d0, 1.d0 /))
        etdsdc_s%m   = sdc_reference_n-1
        relcost      = size(etdsdc_s%tau) * (etdsdc_s%m + 1)
        Nt = F/relcost + 1
        call etdsdc(Lambda,N,tspan,y0,Nt,etdsdc_s,y_out,time)
        if(isfinite(y_out)) then
            y_reference = y_reference + y_out
            solution_count = solution_count + 1
        endif
    endif

    if(reference_methods(2)) then ! Use IMEXSDC_N^N-1
        imexsdc_s%tau = chebpts(sdc_reference_n,(/0.d0, 1.d0 /))
        imexsdc_s%m   = sdc_reference_n-1
        relcost       = size(imexsdc_s%tau) * (imexsdc_s%m + 1)
        Nt = F/relcost + 1
        call imexsdc(lambda,N,tspan,y0,Nt,imexsdc_s,y_out,time)
        if(isfinite(y_out)) then
            y_reference = y_reference + y_out
            solution_count = solution_count + 1
        endif
    endif

    if(reference_methods(3)) then ! Use ETDRK4
        relcost = 4
        Nt = F/relcost + 1
        call etdrk4(Lambda,N,tspan,y0,Nt,y_out,time)
        if(isfinite(y_out)) then
            y_reference = y_reference + y_out
            solution_count = solution_count + 1
        endif
    endif

    y_reference = y_reference/real(solution_count,8)
    write (*,"(a)") "done."

    call write_cvector(y_reference, results_dir//eqn_name//"/reference.txt")

    ! === Run Numerical tests ==================================================
    write (*,"(a)") "Running Numerical Tests... "
    num_methods = size(ETD_SDC_orders) + size(ETD_AB_orders) + 1 + 5 * size(ETD_PBM_qs)
    n_Fs = size(Fs)
    allocate(times(n_Fs,num_methods), errors(n_Fs,num_methods), Nts(n_Fs,num_methods), hs(n_Fs,num_methods))
    do i=1,n_Fs
        F = Fs(i)
        shift = 0
        write (*,"(a,I4,a,I4,a,I10)") "     F (",i,"/",n_Fs,") = ",F
        ! Run ESDC Methods
        do j=1,size(ETD_SDC_orders)
            nn = ETD_SDC_orders(j)
            etdsdc_s%tau  = chebpts(nn, [0.d0, 1.d0])
            etdsdc_s%m    = nn - 1
            relcost       = nn + (nn - 1) ** 2
            Nt            = F / (relcost) + 1
            call etdsdc(Lambda, N, tspan, y0, Nt, etdsdc_s, y_out, time)
            Nts(i,    j + shift) = Nt
            errors(i, j + shift) = error_filter(y_out, y_reference)
            times(i,  j + shift) = time(2)
        enddo
        shift = shift + size(ETD_SDC_orders)
        ! Run EAB Methods
        do j=1,size(ETD_AB_orders)
            nn = ETD_AB_orders(j)
            etdab_s%n = nn
            relcost   = 1
            Nt        = F / (relcost) + 1
            call etdab(lambda, N, tspan, y0, Nt, etdab_s, y_out, time)
            Nts(i,    j + shift) = Nt
            errors(i, j + shift) = error_filter(y_out, y_reference)
            times(i,  j + shift) = time(2)
        enddo
        shift = shift + size(ETD_AB_orders)
        ! Run ERK4
        relcost = 4
        Nt = F/(relcost) + 1
        call etdrk4(Lambda,N,tspan,y0,Nt,y_out,time)
        Nts(i,    1 + shift) = Nt
        errors(i, 1 + shift) = error_filter(y_out, y_reference)
        times(i,  1 + shift) = time(2)
        shift = shift + 1
        if(omp_active .eqv. .TRUE.) then
            ! Run OMP EPBM Methods (Parallel - alpha = 2, m = 0)
            do j=1,size(ETD_PBM_qs)
                q = ETD_PBM_qs(j)
                etdpbm_F1_omp_s%z     = [-1.0_dp, legpts(q - 1, [-1.0_dp, 1.0_dp])]
                etdpbm_F1_omp_s%b     = constvec_r(q, -1.0_dp)
                etdpbm_F1_omp_s%alpha = 2.0_dp
                etdpbm_F1_omp_s%m     = 0
                etdpbm_F1_omp_s%mS    = [1]
                relcost            = q
                Nt = F/(relcost) + 1
                call etdpbm_F1_omp(lambda, N, tspan, y0, Nt, etdpbm_F1_omp_s, y_out, time)
                Nts(i,    j + shift) = Nt
                errors(i, j + shift) = error_filter(y_out, y_reference)
                times(i,  j + shift) = time(2)
            enddo
            shift = shift + size(ETD_PBM_qs)
            ! Run OMP EPBM Methods (Parallel - alpha = 2, m = 1)
            do j=1,size(ETD_PBM_qs)
                q = ETD_PBM_qs(j)
                etdpbm_F1_omp_s%z     = [-1.0_dp, legpts(q - 1, [-1.0_dp, 1.0_dp])]
                etdpbm_F1_omp_s%b     = constvec_r(q, -1.0_dp)
                etdpbm_F1_omp_s%alpha = 2.0_dp
                etdpbm_F1_omp_s%m     = 1
                etdpbm_F1_omp_s%mS    = [1]
                relcost            = q
                Nt = F/(relcost) + 1
                call etdpbm_F1_omp(lambda, N, tspan, y0, Nt, etdpbm_F1_omp_s, y_out, time)
                Nts(i,    j + shift) = Nt
                errors(i, j + shift) = error_filter(y_out, y_reference)
                times(i,  j + shift) = time(2)
            enddo
            shift = shift + size(ETD_PBM_qs)
            ! Run OMP EPBM Methods (Parallel - alpha = 1, m = 0)
            do j=1,size(ETD_PBM_qs)
                q = ETD_PBM_qs(j)
                etdpbm_F1_omp_s%z     = [-1.0_dp, legpts(q - 1, [-1.0_dp, 1.0_dp])]
                etdpbm_F1_omp_s%b     = constvec_r(q, -1.0_dp)
                etdpbm_F1_omp_s%alpha = 1.0_dp
                etdpbm_F1_omp_s%m     = 0
                etdpbm_F1_omp_s%mS    = [1]
                relcost            = q
                Nt = F/(relcost) + 1
                call etdpbm_F1_omp(lambda, N, tspan, y0, Nt, etdpbm_F1_omp_s, y_out, time)
                Nts(i,    j + shift) = Nt
                errors(i, j + shift) = error_filter(y_out, y_reference)
                times(i,  j + shift) = time(2)
            enddo
            shift = shift + size(ETD_PBM_qs)
            ! Run EPBM Methods (Serial - alpha = 2, m = 0)
            do j=1,size(ETD_PBM_qs)
                q = ETD_PBM_qs(j)
                etdpbm_F1_s%z     = [-1.0_dp, legpts(q - 1, [-1.0_dp, 1.0_dp])]
                etdpbm_F1_s%b     = constvec_r(q, -1.0_dp)
                etdpbm_F1_s%alpha = 2.0_dp
                etdpbm_F1_s%m     = 0
                etdpbm_F1_s%mS    = [1]
                relcost = q
                Nt = F/(relcost) + 1
                call etdpbm_F1(lambda, N, tspan, y0, Nt, etdpbm_F1_s, y_out, time)
                Nts(i,    j + shift) = Nt
                errors(i, j + shift) = error_filter(y_out, y_reference)
                times(i,  j + shift) = time(2)
            enddo
            shift = shift + size(ETD_PBM_qs)
            ! Run EPBM Methods - Imaginary Nodes (Serial - alpha = 1, m = 0)
            do j=1,size(ETD_PBM_qs)
                q = ETD_PBM_qs(j)
                etdpbm_F1z_omp_s%z     = [ (-1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp) * linspace(-1.0_dp, 1.0_dp, q - 1)]
                etdpbm_F1z_omp_s%b     = constvec_r(q, -1.0_dp)
                etdpbm_F1z_omp_s%alpha = 2.0_dp
                etdpbm_F1z_omp_s%m     = 0
                etdpbm_F1z_omp_s%mS    = [1]
                relcost = q
                Nt = F/(relcost) + 1
                call etdpbm_F1z_omp(lambda, N, tspan, y0, Nt, etdpbm_F1z_omp_s, y_out, time)
                Nts(i,    j + shift) = Nt
                errors(i, j + shift) = error_filter(y_out, y_reference)
                times(i,  j + shift) = time(2)
            enddo
        else
            ! Run EPBM Methods (Serial)
            do j=1,size(ETD_PBM_qs)
                q = ETD_PBM_qs(j)
                etdpbm_F1_s%z     = [-1.0_dp, legpts(q - 1, [-1.0_dp, 1.0_dp])]
                etdpbm_F1_s%b     = constvec_r(q, -1.0_dp)
                etdpbm_F1_s%alpha = 2.0_dp
                etdpbm_F1_s%m     = 0
                etdpbm_F1_s%mS    = [1]
                relcost = q
                Nt = F/(relcost) + 1
                call etdpbm_F1(lambda, N, tspan, y0, Nt, etdpbm_F1_s, y_out, time)
                Nts(i,    j + shift) = Nt
                errors(i, j + shift) = error_filter(y_out, y_reference)
                times(i,  j + shift) = time(2)
            enddo
        endif

        ! Save Partial Results
        call write_rmatrix(Nts,         results_dir//eqn_name//"/Nts.txt")
        call write_rmatrix(hs,          results_dir//eqn_name//"/hs.txt")
        call write_rmatrix(errors,      results_dir//eqn_name//"/errors.txt")
        call write_rmatrix(times,       results_dir//eqn_name//"/times.txt")
        call write_rvector(Fs,          results_dir//eqn_name//"/Fs.txt")

    enddo
    write (*,"(a)") "done."
    ! Store stepsizes
    hs = (tspan(2) - tspan(1))/Nts

    ! Save Full Results
    call write_rmatrix(Nts,         results_dir//eqn_name//"/Nts.txt")
    call write_rmatrix(hs,          results_dir//eqn_name//"/hs.txt")
    call write_rmatrix(errors,      results_dir//eqn_name//"/errors.txt")
    call write_rmatrix(times,       results_dir//eqn_name//"/times.txt")
    call write_rvector(Fs,          results_dir//eqn_name//"/Fs.txt")

end program SiscExperiment
