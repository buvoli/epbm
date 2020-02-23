classdef epiABConst < ExponentialIntegratorConst
    
    properties
        graph_line_style = {};
        
    end
    
    properties(SetAccess = protected)
        name = '';
        description = '';
        order = [];
        starting_times   = []; 
        step_backwards_for_ic = false;   % if true, initial conditions will be computed behind first timestep
        wts   = [];
    end
     
    properties(Access = private)
        q
        z
    end
    
    methods
        function this = epiABConst(options)
            if(nargin == 0)
                options = struct();
            end
            default_field_value_pairs = {{'order', 2}, {'step_backwards_for_ic', false}};
            options = setDefaultOptions(options, default_field_value_pairs);
            this = this@ExponentialIntegratorConst(options);
            this.step_backwards_for_ic = options.step_backwards_for_ic;
            this.setOrder(options.order);
        end
        
        function setOrder(this, order)
            this.order           = order;
            this.name            = ['EPIAB', num2str(order)];
            if(this.step_backwards_for_ic)
                this.starting_times = linspace(-order + 1, 0, order);
            else
                this.starting_times = linspace(0, order - 1, order);
            end
            this.wts = weights(order - 1, linspace(0, order - 1, order), order - 1);
        end
        
    end
    
    methods (Access = protected)
        
        function setStepsize(this, problem)
            if(this.step_backwards_for_ic)
                this.h = (problem.tspan(end) - problem.tspan(1)) / (this.num_timesteps);
            else
                this.h = (problem.tspan(end) - problem.tspan(1)) / (this.num_timesteps + this.order - 1);
            end
        end
        
        function [step_struct, y_in] = initStepStruct(this, t_in, y_in, problem)
            
            % -- define remainder function -----------------------------------------------------------------------------
            step_struct = struct(...
                'b', zeros(size(y_in,1), this.order + 1), ...
                'r', zeros(size(y_in,1), this.order) ...
            );
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, final_step)
            
            step_start_time = tic;
            
            h    = this.h;
            y_n  = y_in(:, end);
            hF_n = h * problem.RHS(y_n);                            % calculate right hand side
            hJ_n = h * problem.J(y_n);                              % calculate Jacobian
            hR   = @(y) h * problem.RHS(y) - hF_n - hJ_n*(y - y_n); % remainder function
            for j = 1 : this.order
                step_struct.r(:, j) = hR(y_in(:,j));
            end
            % -- form vectors b ----------------------------------------------------------------------------------------
            step_struct.b(:, 2:end) = step_struct.r * this.wts;
            step_struct.b(:, 2) = step_struct.b(:,2) + hF_n;
            % -- advance solution --------------------------------------------------------------------------------------
            y_new = y_n + this.phi_evaluator.solve(1, hJ_n, step_struct.b);
            y_out = [y_in(:,2:end), y_new]; 
            t_out = t_in + h;
            
            
            if(final_step)
                y_out = y_out(:, end);
            end
            
            % Update counters.
            this.step_stats.recordStep(toc(step_start_time));
            
        end
        
    end
    
end