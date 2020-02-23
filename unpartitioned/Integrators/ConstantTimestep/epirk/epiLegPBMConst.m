classdef epiLegPBMConst < ExponentialIntegratorConst
    
    properties
        graph_line_style = {};
        m     = 0;
        alpha = 2;
    end
    
    properties(SetAccess = protected)
        name = '';
        description = '';
        order = [];
        starting_times   = []; 
        wts   = [];
    end
     
    properties(Access = private)
        q
        z
    end
    
    methods
        function this = epiLegPBMConst(options)
            if(nargin == 0)
                options = struct();
            end
            default_field_value_pairs = {
                {'q', 3}
                {'m', 0}
                {'alpha', 2}
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            this = this@ExponentialIntegratorConst(options);
            this.m = options.m;
            this.alpha = options.alpha;
            this.setQ(options.q);
        end
        
        function setQ(this, q)
            leg = legpts(q - 1, [0 2]) / this.alpha;
            if(q > 3)
                this.order = q + 1;
            else
                this.order = q;
            end
            this.q               = q;
            this.name            = ['EPIPBM', num2str(q)];
            this.wts             = weights(0, leg, q - 2);
            this.starting_times  = 0;
            this.z               = [0, leg(:).'];           
        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct, y_in] = initStepStruct(this, t_in, y_in, problem)
            
            % -- define remainder function -----------------------------------------------------------------------------
            h    = this.h;
            w    = this.wts;
            hF_n = h * problem.RHS(y_in);                            % calculate right hand side
            hJ_n = h * problem.J(y_in);                              % calculate Jacobian
            hR   = @(y) h * problem.RHS(y) - hF_n - hJ_n*(y - y_in); % remainder function
            
            step_struct = struct(...
                'b', zeros(size(y_in,1), this.q), ...
                'r', zeros(size(y_in,1), this.q - 1), ...
                'hF_n', hF_n, ...
                'hJ_n', hJ_n ...
            );
        
            % -- iterate for initial conditions ------------------------------------------------------------------------
            for i = 1 : this.q
                % -- form vectors b ------------------------------------------------------------------------------------
                step_struct.b(:, 2 : end) = step_struct.r * w;
                step_struct.b(:, 2) = step_struct.b(:,2) + hF_n;
                % -- advance solution ----------------------------------------------------------------------------------
                y_out = y_in + this.phi_evaluator.solve(this.z(2:end), hJ_n, step_struct.b);               
                for j = 1 : this.q - 1
                    step_struct.r(:, j) = hR(y_out(:,j));
                end                
            end
        
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, ~)
            
            step_start_time = tic;
            
            h    = this.h;
            % -- form vectors b ----------------------------------------------------------------------------------------
            step_struct.b(:, 2:end) = step_struct.r * this.wts;
            step_struct.b(:, 2) = step_struct.b(:,2) + step_struct.hF_n;
            % -- advance solution --------------------------------------------------------------------------------------
            [y_out, exit_flag] = this.phi_evaluator.solve(this.z + 1, step_struct.hJ_n, step_struct.b);
            if(exit_flag == false)
                t_out = NaN; y_out = NaN; return;
            end
            y_out = y_in + y_out;
            t_out = t_in + h;
            % -- define remainder function -----------------------------------------------------------------------------
            y_n  = y_out(:, 1);
            hF_n = h * problem.RHS(y_n);                            % calculate right hand side
            hJ_n = h * problem.J(y_n);                              % calculate Jacobian
            hR   = @(y) h * problem.RHS(y) - hF_n - hJ_n*(y - y_n); % remainder function
            % -- evaluate remainder terms ------------------------------------------------------------------------------
            for j = 2 : this.q
                step_struct.r(:, j-1) = hR(y_out(:,j));
            end
            step_struct.hF_n = hF_n;
            step_struct.hJ_n = hJ_n;
            y_out = y_out(:,1);            
            
            % -- optional iteration ------------------------------------------------------------------------------------
            for i = 1 : this.m
                % -- form vectors b ------------------------------------------------------------------------------------
                step_struct.b(:, 2 : end) = step_struct.r * this.wts;
                step_struct.b(:, 2) = step_struct.b(:,2) + hF_n;
                % -- advance solution ----------------------------------------------------------------------------------
                ys_out = y_out + this.phi_evaluator.solve(this.z(2:end), hJ_n, step_struct.b);               
                for j = 1 : this.q - 1
                    step_struct.r(:, j) = hR(ys_out(:,j));
                end                
            end
            
            % Update counters.
            this.step_stats.recordStep(toc(step_start_time));
            
        end
        
    end
    
end