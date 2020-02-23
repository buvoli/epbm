classdef epiSDCConst < ExponentialIntegratorConst
    
    properties
        graph_line_style = {};
        use_prev_b = false;
    end
    
    properties(SetAccess = protected)
        name = '';
        description = '';
        order = [];
        starting_times = [0]; 
        nodes      = {};
        num_nodes  = [];
        num_sweeps = [];
        FDweights  = {};
    end
    
    methods
        function this = epiSDCConst(options)
            if(nargin == 0)
                options = struct();
            end
            default_field_value_pairs = {
                {'nodes', (chebpts(3) + 1) / 2,} ...
                {'sweeps', 3} ...
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            this = this@ExponentialIntegratorConst(options);
            this.setNodes(options.nodes, options.sweeps);
        end
        
        function setNodes(this, nodes, sweeps)
            if(nargin == 2)
                sweeps = 0;
            end
            sweeps = round(sweeps); % ensure integer value ----------------
            this.nodes = nodes;
            this.num_nodes = length(nodes);
            this.num_sweeps = sweeps;
            this.order = min(sweeps + 1, length(nodes));
            this.name = ['EPISDC_', num2str(length(nodes)), '^', num2str(sweeps)];
            % -- compute weights ------------------------------------------
            n   = length(nodes);
            this.FDweights = cell(n - 1, 1);
            for i=1:n-1
                q = (nodes - nodes(i));                                     % translated quadrature points
                this.FDweights{i} = weights(0,q,n-1);                       % finite difference matrix
            end
        end
        
        function setNumSweeps(this, num_sweeps)
            this.setNodes(this.nodes, num_sweeps);
        end
    end
    
    methods (Access = protected)
        
        function [step_struct, y_in] = initStepStruct(this, t_in, y_in, problem)
            ode_dim = size(y_in, 1);
            step_struct = struct(                           ...
                'y',    zeros(ode_dim, this.num_nodes),     ...             % struct for storing solution at substep
                'hR',   zeros(ode_dim, this.num_nodes),     ...             % struct for storing remainder at substeps
                'hF',   zeros(ode_dim, this.num_nodes),     ...             % struct for storing F evaluated at substeps
                'b',    zeros(ode_dim, this.num_nodes + 1), ...             % struct for storing vectors that will be multiplied by phi functions  
                'eta',  diff(this.nodes)                    ...             % normalized substeps
            );
            if(this.use_prev_b)
                step_struct.prev_b = cell(this.num_nodes - 1);
                step_struct.prev_projections = cell(this.num_nodes - 1);
            end
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, final_step)
            
            step_start_time = tic;
            
            % -- store parameters into local variables --------------------
            h   = this.h;
            wts = this.FDweights;
            num_nodes  = this.num_nodes;
            num_sweeps = this.num_sweeps;
            
            % -- Define local remainder Function --------------------------
            y_n  = y_in;
            hF_n = h * problem.RHS(y_n);
            hJ_n = h * problem.J(y_n);                            % calculate Jacobian
            hR   = @(y, hFy) hFy - hF_n - hJ_n*(y - y_n);         % remainder function
            

            
            if(this.use_prev_b) % compute phi functions applied to delta vector
            
                % -- Exponential Euler Step -----------------------------------
                step_struct.y(:, 1)  = y_n;
                step_struct.hF(:, 1) = hF_n;
                step_struct.b(:)     = 0;
                for j = 1 : num_nodes - 1
                    step_struct.b(:, 2) = step_struct.hF(:, j);   % form b vector. Note b(:, 1) = 0.
                    % -- save projection vector and result ---
                    step_struct.prev_b{j} = step_struct.b; % save projection vector
                    step_struct.prev_projections{j} = this.phi_evaluator.solve(step_struct.eta(j), hJ_n, step_struct.b(:,1:2));
                    step_struct.y(:, j + 1)  = step_struct.y(:,j) + step_struct.prev_projections{j};
                    step_struct.hF(:, j + 1) = h * problem.RHS(step_struct.y(:, j + 1));     % store hF at each substep.
                    step_struct.hR(:, j + 1) = hR(step_struct.y(:, j + 1), step_struct.hF(:, j + 1)); % store hR at each substep. Note: hR(:, 1) = 0.
                end
                
                % -- Exponential Correction Sweeps -------------------------------------------------------------------------
                for k = 1 : num_sweeps
                    % -- correction sweep ----------------------------------------------------------------------------------
                    for j = 1 : num_nodes - 1
                        % -- form b vectors --------------------------------------------------------------------------------
                        step_struct.b(:, 2 : end) = step_struct.hR * wts{j};
                        step_struct.b(:, 2) = step_struct.b(:, 2) + hF_n + hJ_n * (step_struct.y(:,j) - y_n);
                        if(j > 1) % -- add contribution from h \varphi_1(h J) (F_{n,j}^{[k+1]} - F_{n,j}^{[k]}) --- NOTE: for j = 1 term is zero since F_{n,j}^{[k+1]} = F_{n,j}^{[k]}
                            hF_new = h * problem.RHS(step_struct.y(:,j));
                            step_struct.b(:, 2)  = step_struct.b(:, 2) + (hF_new - step_struct.hF(:,j));
                            step_struct.hF(:, j) = hF_new;
                        end
                        % -- step solution ---------------------------------------------------------------------------------
                        delta = step_struct.b - step_struct.prev_b{j};
                        step_struct.prev_b{j} = step_struct.b;
                        step_struct.prev_projections{j} = step_struct.prev_projections{j} + this.phi_evaluator.solve(step_struct.eta(j), hJ_n, delta);
                        step_struct.y(:,j+1) = step_struct.y(:,j) + step_struct.prev_projections{j};
                    end
                    % -- update final f evaluation and hR values at each substep ---------------------------------------
                    step_struct.hF(:, num_nodes) = h * problem.RHS(step_struct.y(:,num_nodes));
                    for j = 2 : num_nodes
                        step_struct.hR(:,j) = hR(step_struct.y(:, j), step_struct.hF(:, j));
                    end
                end
                
            else
            
                % -- Exponential Euler Step ----------------------------------------------------------------------------
                step_struct.y(:, 1)  = y_n;
                step_struct.hF(:, 1) = hF_n;
                for j = 1 : num_nodes - 1
                    step_struct.b(:, 2)      = step_struct.hF(:, j);   % form b vector. Note b(:, 1) = 0.
                    [projection, clean_exit] = this.phi_evaluator.solve(step_struct.eta(j), hJ_n, step_struct.b(:, 1:2));
                    % -- Check for emergency exit ----------------------------------------------------------------------
                    step_struct.y(:, j + 1)  = step_struct.y(:, j) + projection;
                    emergency_exit = ( (~clean_exit) || any(isinf(step_struct.y(:,j+1))) || any(isnan(step_struct.y(:,j+1))) );
                    if(emergency_exit)
                        break;
                    end
                    % -- Comptute F and R ------------------------------------------------------------------------------
                    step_struct.hF(:, j + 1) = h * problem.RHS(step_struct.y(:, j + 1));     % store hF at each substep.
                    step_struct.hR(:, j + 1) = hR(step_struct.y(:, j + 1), step_struct.hF(:, j + 1)); % store hR at each substep. Note: hR(:, 1) = 0.
                end
                
                if(~emergency_exit)
                    % -- Exponential Correction Sweeps ---------------------------------------------------------------------
                    for k = 1 : num_sweeps

                        % -- correction sweep ------------------------------------------------------------------------------
                        for j = 1 : num_nodes - 1
                            % -- form b vectors ----------------------------------------------------------------------------
                            step_struct.b(:, 2 : end) = step_struct.hR * wts{j};
                            step_struct.b(:, 2) = step_struct.b(:, 2) + hF_n + hJ_n * (step_struct.y(:,j) - y_n);
                            if(j > 1) % -- add contribution from h \varphi_1(h J) (F_{n,j}^{[k+1]} - F_{n,j}^{[k]}) --- NOTE: for j = 1 term is zero since F_{n,j}^{[k+1]} = F_{n,j}^{[k]}
                                hF_new = h * problem.RHS(step_struct.y(:,j));
                                step_struct.b(:, 2) = step_struct.b(:, 2) + (hF_new - step_struct.hF(:,j));
                                step_struct.hF(:, j) = hF_new;
                            end
                            % -- step solution -----------------------------------------------------------------------------
                            [projection, clean_exit] = this.phi_evaluator.solve(step_struct.eta(j), hJ_n, step_struct.b);
                            step_struct.y(:,j+1) = step_struct.y(:,j) + projection;
                            % -- check for emergency exit ------------------------------------------------------------------
                            emergency_exit = ( (~clean_exit) || any(isinf(step_struct.y(:,j+1))) || any(isnan(step_struct.y(:,j+1))) );
                            if(emergency_exit)
                                break;
                            end
                        end
                        if(emergency_exit)
                            break;
                        end
                        % -- update final f evaluation and hR values at each substep ---------------------------------------
                        step_struct.hF(:, num_nodes) = h * problem.RHS(step_struct.y(:,num_nodes));
                        for j = 2 : num_nodes
                            step_struct.hR(:,j) = hR(step_struct.y(:, j), step_struct.hF(:, j));
                        end
                    end
                end
            
            end
            
            % -- Set Solution at Next Timestep -------------------------------------------------------------------------
            if(emergency_exit)
                t_out = NaN;
                y_out = NaN;
            else
                y_out = step_struct.y(:, end);
                t_out = t_in + h;
            end
            
            % Update counters.
            this.step_stats.recordStep(toc(step_start_time));
            
        end
        
    end
    
end