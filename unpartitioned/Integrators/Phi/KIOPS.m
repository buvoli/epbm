classdef KIOPS < PhiEvaluator
    %KIOPS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        stats = PhiStats() % IMPLEMENT STATS
        tolerance;
        krylov_max_dimension;
        krylov_min_dimension;
        krylov_start_dimension;
        min_tau;
    end
    
    methods
        
        function this = KIOPS(options)
            %KIOPS Construct an instance of this class
            %   Detailed explanation goes here
            if(nargin == 0)
            	options = struct();
            end
            options = this.DefaultOptions(options);
            % set options
            this.tolerance              = options.tolerance;
            this.krylov_max_dimension   = options.krylov_max_dimension;
            this.krylov_min_dimension   = options.krylov_min_dimension;
            this.krylov_start_dimension = options.krylov_start_dimension;
            this.min_tau                = options.min_tau;
        end
        
        function [y, exit_flag] = phi(this, n, t, A) 
            error('KIOPS does not support phi');
        end
        
        function [y, exit_flag] = phis(this, n, t, A) 
            error('KIOPS does not support phis');
        end
        
        function [y, exit_flag] = phib(this, n, t, A, b)
            % Initialize phi_n(t*A) b
            [y, ~, ~, exit_flag] = this.kiops_raw(transpose(t(:)), A, [zeros(size(b,1), n) b], this.tolerance, this.krylov_start_dimension, this.krylov_min_dimension, this.krylov_max_dimension, true); 
        end
        
        function [y, exit_flag] = phibs(this, t, A, bs)
            % Initialize phi_0(t*A) b(:, 1) + phi_1(t * A) b(:, 2) + ... + phi_{length(b)-1}(t * A) b(:,end)
            [y, ~, ~, exit_flag] = this.kiops_raw(transpose(t(:)), A, bs, this.tolerance, this.krylov_start_dimension, this.krylov_min_dimension, this.krylov_max_dimension, true); 
        end
        
        function [y, exit_flag] = solve(this, t, A, bs)
            % compute solution to the differential equation
            %       y' = A y + bs(:, 1) + t * bs(:, 2) / 1! + t^2 * bs(:, 3) / 2! + ...
            % by initialize phi_0(t*A) b(:, 1) + t phi_1(t * A) b(:, 2) + ... + tphi_{length(b)-1}(t * A) b(:,end)
            [y, ~, ~, exit_flag] = this.kiops_raw(transpose(t(:)), A, bs, this.tolerance, this.krylov_start_dimension, this.krylov_min_dimension, this.krylov_max_dimension, false);
        end
        
    end
    
    methods(Access = private)
        
        function options = DefaultOptions(this, options)
            if(~isfield(options, 'tolerance'))
                options.tolerance = 1e-10;
            end
            if(~isfield(options, 'krylov_max_dimension'))
                options.krylov_max_dimension = 128;
            end
            if(~isfield(options, 'krylov_min_dimension'))
                options.krylov_min_dimension = 8;
            end
            if(~isfield(options, 'krylov_start_dimension'))
                options.krylov_start_dimension = 8;
            end
            if(~isfield(options, 'min_tau'))
                options.min_tau = eps;
            end
        end 
        
        function [w, m_ret, stats, clean_exit] = kiops_raw(this, t, A, u, tol, m_init, mmin, mmax, trimcoef)
            % Evaluates a linear combinaton of the phi functions
            % evaluated at tA acting on vectors from u, that is
            %
            % w = phi_0(tA) u(:, 1) + phi_1(tA) u(:, 2) + phi_2(tA) u(:, 3) + ...
            %
            % The size of the Krylov subspace is changed dynamically
            % during the integration. The Krylov subspace is computed
            % using the incomplete orthogonalization method.
            %
            % PARAMETERS:
            %   t          - constant value represent time.
            %   A          - the matrix argument of the phi functions.
            %   u          - the matrix with columns representing the vectors to be
            %                multiplied by the phi functions.
            %
            % OPTIONAL PARAMETERS:
            %   tol        - the convergence tolerance required.
            %   m_init     - an estimate of the appropriate Krylov size.
            %   mmin, mmax - let the Krylov size vary between mmin and mmax
            %   trimcoef   - If true, return the linear combination above. If set
            %                to false, return the linear combinaison
            %                w = phi_0(tA) u(:, 1) + t phi_1(tA) u(:, 2) + t^2 phi_2(tA) u(:, 3) + ...
            %                This allows for compatibility with other solvers (eg. phipm)
            %
            % RETURNS:
            %   w        - the linear combination of the phi functions
            %              evaluated at tA acting on the vectors from u.
            %   m        - the Krylov size of the last substep.
            %   stats(1) - number of substeps
            %   stats(2) - number of rejected steps
            %   stats(3) - number of Krylov steps
            %   stats(4) - number of matrix exponentials
            
            % n is the size of the original problem
            % p is the highest indice of the phi functions
            [n, ppo] = size(u);
            p = ppo - 1;
            
            if p == 0
                p = 1;
                % Add extra column of zeros
                u = [u, zeros(size(u))];
            end
            
            if ~(isa(A, 'function_handle'))
                A = @(vec) A*vec;
            end
            
            % Krylov parameters
            orth_len = 2;
            
            % Check inputs
            if nargin < 8
                trimcoef = false;
                if nargin < 7
                    mmax = 128;
                    if nargin < 6
                        mmin = 10;
                        if nargin < 5
                            m_init = mmin;
                            if nargin < 4
                                tol = 1.0e-7;
                                if nargin < 3
                                    error('Not enough input arguments.');
                                end
                            end
                        end
                    end
                end
            end
            
            % We only allow m to vary between mmin and mmax
            m = max(mmin, min(m_init, mmax));
            
            % Preallocate matrix
            V = zeros(n + p, mmax + 1);
            H = zeros(mmax + 1, mmax + 1);
            
            step    = 0;
            krystep = 0;
            ireject = 0;
            reject  = 0;
            exps    = 0;
            sgn     = sign(t(end));
            t_now   = 0;
            t_out   = abs(t(end));
            happy   = false;
            j       = 0;
            min_tau = this.min_tau;
            clean_exit = true;
            
            numSteps = size(t, 2);
            
            % Initial condition
            w     = zeros(n, numSteps);
            w_aug = zeros(p, 1);
            w(:, 1) = u(:, 1);
            
            % Normalization factors
            normU = norm(u(:, 2:end),1);
            if ppo > 1 && normU > 0
                ex = ceil(log2(normU));
                nu = 2^(-ex);
                mu = 2^(ex);
            else
                nu = 1;
                mu = 1;
            end
            
            % Flip the rest of the u matrix
            u_flip = nu*fliplr(u(:, 2:end));
            
            % Compute and initial starting approximation for the timestep
            tau = t_out;
            
            % Setting the safety factors and tolerance requirements
            if t_out > 1
                gamma = 0.2;
                gamma_mmax = 0.1;
            else
                gamma = 0.9;
                gamma_mmax = 0.6;
            end
            delta = 1.4;
            
            % Used in the adaptive selection
            oldm = NaN; oldtau = NaN; omega = NaN;
            orderold = true; kestold = true;
            
            l=1;
            % Iterate until we reach the final time t
            while t_now < t_out
                
                % Compute necessary starting information
                if j == 0
                    
                    % Update the last part of w
                    for k=1:p-1
                        i = p - k;
                        w_aug(k) = (t_now^i)/factorial(i) * mu;
                    end
                    w_aug(p) = mu;
                    
                    % Initialize the matrices V and H
                    H(:, :) = 0;
                    
                    % Normalize initial vector (this norm is nonzero)
                    beta = sqrt( w(:,l)' * w(:,l) + w_aug' * w_aug );
                    
                    % The first Krylov basis vector
                    V(1:n, 1)     = (1/beta) * w(:,l);
                    V(n+1:n+p, 1) = (1/beta) * w_aug;
                    
                end
                
                % Incomplete orthogonalization process
                while j < m
                    
                    j = j + 1;
                    
                    % Augmented matrix - vector product
                    V(1:n      , j + 1) = A( V(1:n, j) ) + u_flip * V(n+1:n+p, j);
                    V(n+1:n+p-1, j + 1) = V(n+2:n+p, j);
                    V(end      , j + 1) = 0;
                    
                    % Modified Gram-Schmidt
                    for i = max(1,j-orth_len+1):j
                        H(i, j) = V(:, i)' * V(:, j + 1);
                        V(:, j + 1) = V(:, j + 1) - H(i, j) * V(:, i);
                    end
                    
                    nrm = norm(V(:, j + 1));
                    
                    % Happy breakdown
                    if nrm < tol
                        happy = true;
                        break;
                    end
                    
                    H(j + 1, j) = nrm;
                    V(:, j + 1) = (1/nrm) * V(:, j + 1);
                    
                    krystep = krystep + 1;
                    
                end
                
                % To obtain the phi_1 function which is needed for error estimate
                H(1, j + 1) = 1;
                
                % Save h_j+1,j and remove it temporarily to compute the exponential of H
                nrm = H(j + 1, j);
                H(j + 1, j) = 0;
                
                % -- NaN emergency exit condition ---
                subH = H(1:j + 1, 1:j + 1);
                if(any(isnan(subH(:))) || any(isinf(subH(:))))
                    warning('KIOPS: H matrix contains NaN or Inf.')
                    w = NaN;
                    clean_exit = false;
                    break;
                end
                
                % Compute the exponential of the augmented matrix
                F = this.expm_raw(sgn * tau * subH);
                exps = exps + 1;
                
                % Restore the value of H_{m+1,m}
                H(j + 1, j) = nrm;
                
                if happy
                    
                    % Happy breakdown; wrap up
                    omega   = 0;
                    happy   = false;
                    m_new   = m;
                    tau_new = min(t_out - (t_now + tau), tau);
                    
                else
                    
                    % Local truncation error estimation
                    err = abs(beta * nrm * F(j, j + 1));
                    
                    % Error for this step
                    oldomega = omega;
                    omega = t_out * err / (tau * tol);
                    
                    % Estimate order
                    if m == oldm && tau ~= oldtau && ireject >= 1
                        order = max(1, log(omega/oldomega) / log(tau/oldtau));
                        orderold = false;
                    elseif orderold || ireject == 0
                        orderold = true;
                        order = j/4;
                    else
                        orderold = true;
                    end
                    % Estimate k
                    if m ~= oldm && tau == oldtau && ireject >= 1
                        kest = max(1.1, (omega/oldomega)^(1/(oldm-m)));
                        kestold = false;
                    elseif kestold || ireject == 0
                        kestold = true;
                        kest = 2;
                    else
                        kestold = true;
                    end
                    
                    if omega > delta
                        remaining_time = t_out - t_now;
                    else
                        remaining_time = t_out - (t_now + tau);
                    end
                    
                    % Krylov adaptivity
                    
                    same_tau = min(remaining_time, tau);
                    tau_opt  = tau * (gamma / omega)^(1 / order);
                    tau_opt  = min(remaining_time, max(tau/5, min(5*tau, tau_opt)));
                    
                    m_opt = ceil(j + log(omega / gamma) / log(kest));
                    m_opt = max(mmin, min(mmax, max(floor(3/4*m), min(m_opt, ceil(4/3*m)))));
                    
                    if j == mmax
                        if omega > delta
                            m_new = j;
                            tau_new = tau * (gamma_mmax / omega)^(1 / order);
                            tau_new = min(t_out - t_now, max(tau/5, tau_new));
                        else
                            tau_new = tau_opt;
                            m_new = m;
                        end
                    else
                        m_new = m_opt;
                        tau_new = same_tau;
                    end
                    
                end
                
                % Check error against target
                if omega <= delta
                    
                    % Yep, got the required tolerance; update
                    reject = reject + ireject;
                    step = step + 1;
                    
                    % Udate for t in the interval (t_now, t_now + tau)
                    blownTs = 0;
                    nextT = t_now + tau;
                    for k = l:numSteps
                        if abs(t(k)) < abs(nextT)
                            blownTs = blownTs + 1;
                        end
                    end
                    
                    if blownTs ~= 0
                        % Copy current w to w we continue with.
                        w(:,l + blownTs) = w(:,l);
                        
                        for k = 0:blownTs - 1
                            tauPhantom = t(l+k) - t_now;
                            F2 = this.expm_raw(sgn * tauPhantom * H(1:j, 1:j));
                            w(:, l+k) = beta * V(1:n, 1:j) * F2(1:j, 1);
                        end
                        
                        % Advance l.
                        l = l + blownTs;
                    end
                    
                    % Using the standard scheme
                    w(:, l) = beta * V(1:n, 1:j) * F(1:j, 1);
                    
                    % Update t
                    t_now = t_now + tau; % time update
                    
                    j = 0;
                    ireject = 0;
                    
                else
                    
                    % Nope, try again
                    ireject = ireject + 1;
                    
                    % Restore the original matrix
                    H(1, j + 1) = 0;
                end
                
                oldtau = tau;
                tau    = tau_new;
                
                oldm = m;
                m    = m_new;
                
                % -- check if tau is above min tolerance ---------------------------------------------------------------
                %tau
                %min_tau
                %(t_now - t_out)
                %[tau, t_now, t_out min_tau]
                if((oldm == m) && (tau < min_tau) && (t_out - t_now) > min_tau)
                    warning('KIOPS: tau below minimum tolerance; Returning NaN.')
                    w = NaN;
                    clean_exit = false;
                    break;
                end
                
            end
                        
            if ( (t(1) ~= 1) && (trimcoef) && (clean_exit) ) % if trim is true, attempt to remove factors of t
                active_phi_n = find(any(u)) - 1;
                if(isempty(active_phi_n))
                    active_phi_n = 0; % ensure that active_phi is a number
                end

                if(length(active_phi_n) > 1)
                    warning('coefficients cannot be trimmed (more than on phi function present)!');
                else
                    w = bsxfun(@times, w, t.^(-active_phi_n));
                end
            end
            
            
            
            
            m_ret=m;
            
            stats = [step, reject, krystep, exps];
            
        end
        
        function F = expm_raw(this, A)
            %   Compute the matrix exponential of A using a scaling and squaring
            %   algorithm with a Pade approximation.
            %
            %   Reference:
            %   N. J. Higham, The scaling and squaring method for the matrix
            %   exponential revisited. SIAM J. Matrix Anal. Appl.,
            %   26(4) (2005), pp. 1179-1193.
            %
            switch class(A)
                case 'double'
                    m_vals = [3 5 7 9 13];
                    % theta_m for m=1:13.
                    theta = [%3.650024139523051e-008
                        %5.317232856892575e-004
                        1.495585217958292e-002  % m_vals = 3
                        %8.536352760102745e-002
                        2.539398330063230e-001  % m_vals = 5
                        %5.414660951208968e-001
                        9.504178996162932e-001  % m_vals = 7
                        %1.473163964234804e+000
                        2.097847961257068e+000  % m_vals = 9
                        %2.811644121620263e+000
                        %3.602330066265032e+000
                        %4.458935413036850e+000
                        5.371920351148152e+000];% m_vals = 13
                case 'single'
                    m_vals = [3 5 7];
                    % theta_m for m=1:7.
                    theta = [%8.457278879935396e-004
                        %8.093024012430565e-002
                        4.258730016922831e-001  % m_vals = 3
                        %1.049003250386875e+000
                        1.880152677804762e+000  % m_vals = 5
                        %2.854332750593825e+000
                        3.925724783138660e+000];% m_vals = 7
                otherwise
                    error(message('MATLAB:expm:inputType'))
            end
            
            normA = norm(A,1);
            if normA <= theta(end)
                
                % no scaling and squaring is required.
                for i = 1:length(m_vals)
                    if normA <= theta(i)
                        F = this.epade_raw(m_vals(i),A);
                        break;
                    end
                end
                
            else
                
                % Scale the matrix A so that a Padé approximant of degree 13 will be accurate
                [t, s] = log2(normA/theta(end));
                s = s - (t == 0.5); % adjust s if normA/theta(end) is a power of 2.
                
                A  = A*(1/2^s);    % Scaling
                
                F = this.epade_raw(m_vals(end),A);
                
                for i = 1:s
                    F = F*F;  % Squaring
                end
                
            end
            
        end
        
        function F = epade_raw(~, m,A)
            %   Pade' approximant to exponential.
            %   Compute the degree M diagonal Pade' approximant to EXP(A),
            %   where M = 3, 5, 7, 9 or 13.
            %   Series are evaluated in decreasing order of powers, which is
            %   in approx. increasing order of maximum norms of the terms.
            n = length(A);
            I = eye(n,class(A));
            A2 = A*A;
            % Evaluate Pade approximant.
            switch m
                case 3
                    c = [120, 60, 12, 1];
                    U = A * (c(4)*A2 + c(2)*I);
                    V = c(3)*A2 + c(1)*I;
                case 5
                    c = [30240, 15120, 3360, 420, 30, 1];
                    A4 = A2*A2;
                    U = A * (c(6)*A4 + c(4)*A2 + c(2)*I);
                    V = c(5)*A4 + c(3)*A2 + c(1)*I;
                case 7
                    c = [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1];
                    A4 = A2*A2;
                    A6 = A2*A4;
                    U = A * (c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*I);
                    V = c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*I;
                case 9
                    c = [17643225600, 8821612800, 2075673600, 302702400, 30270240, ...
                        2162160, 110880, 3960, 90, 1];
                    A4 = A2*A2;
                    A6 = A2*A4;
                    A8 = A4*A4;
                    U = A * (c(10)*A8 + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*I);
                    V = c(9)*A8 + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*I;
                case 13
                    c = [64764752532480000, 32382376266240000, 7771770303897600, ...
                        1187353796428800,  129060195264000,   10559470521600,   ...
                        670442572800,      33522128640,       1323241920,       ...
                        40840800,          960960,            16380,  182,  1];
                    A4 = A2*A2;
                    A6 = A2*A4;
                    U = A * (A6*(c(14)*A6 + c(12)*A4 + c(10)*A2) + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*I);
                    V = A6*(c(13)*A6 + c(11)*A4 + c(9)*A2) + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*I;
            end
            F = (V-U)\(2*U) + I;
        end
        
    end
    
end