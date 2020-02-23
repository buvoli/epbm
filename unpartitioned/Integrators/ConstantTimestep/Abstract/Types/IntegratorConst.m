% NOTE: There should be two kinds of integrators. FixedTimestepIntegrator and VariableTimeIntegratior

% NOTE: There should be ImplicitIntegrator, ExpoenntialIntegrator, ImplicitExponentialIntegrator so IMEX Integrator, ...
% and each type of integrator should have its own combination of linear_solver, nonlinear_solver, phi_initializer, etc

classdef IntegratorConst < handle
    
    properties (Abstract = true)
        graph_line_style
    end
    
    properties (Abstract = true, SetAccess = protected)
        name
        description
        starting_times      % array of required starting times
        order               % method order
    end
    
    properties (Access = protected)
        h
    end
    
    properties
        % -- integrator properties -------------------------------------------------------------------------------------
        num_timesteps = 1;
        starting_integrator;
        starting_integrator_timeout = 600; % max time in seconds for computing each initial condition
        % -- stats objects ---------------------------------------------------------------------------------------------
        step_stats
        rhs_stats
    end
    
    methods (Abstract = true, Access = protected)
        step(this, t_in, y_in, step_struct, problem, final_step);
        initStepStruct(this, t_in, y_in, problem);
    end
    
    methods
        
        function this = IntegratorConst(options)
            
            default_field_value_pairs = { ...
                {'starting_integrator', []} ...
                {'graph_line_style',    {'k.--'}} ...
                {'record_stats',        false} ...
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            
            % -- set props ---------------------------------------------------------------------------------------------
            props = {'starting_integrator', 'graph_line_style'};
            this.setEmptyClassProps(props, options);
            
            % -- initialize stats object -------------------------------------------------------------------------------
            this.step_stats = TimestepStats();
                this.step_stats.record = options.record_stats;            
            this.rhs_stats = RHSStats();
                this.rhs_stats.record = options.record_stats;
        end
        
        function set.num_timesteps(this, num_timesteps)
            %SET.NUM_TIMESTEP ensure that num_timesteps is valid
            input_is_int = floor(num_timesteps) == num_timesteps;
            intput_is_positive = num_timesteps > 0;
            if(input_is_int && intput_is_positive)
                this.num_timesteps = num_timesteps;
            else
                warning('Integrator: num_timesteps must be strictly positive integer');
            end
        end
        
        function [t_out, y_out] = solve(this, problem)
            if(isempty(this.num_timesteps))
                error('Integrator: number of timesteps has not been set.');
            end
            % -- set stepsize ------------------------------------------------------------------------------------------
            this.setStepsize(problem);
            % -- obtain initial conditions -----------------------------------------------------------------------------
            [t_in, y_in] = this.initialConditions(problem);
            % -- solve problem -----------------------------------------------------------------------------------------
            [struct_in, y_in] = this.initStepStruct(t_in, y_in, problem); % note: code was modified on Dec 4, 2018. Must return y_in as well
            num_steps = this.num_timesteps;
            for i = 1 : num_steps
                final_step_flag = (i == num_steps);              
                [t_in, y_in, struct_in] = this.step(t_in, y_in, struct_in, problem, final_step_flag);
                % -- emergency exit conditions ---------------------------------------------------------------------
                if(any(isinf(y_in(:))) || any(isnan(y_in(:))))
                    t_in = NaN;
                    y_in = NaN;
                    warning('Integrator: Emergency exit at step %i', i);
                    break;
                end
            end
            t_out = t_in; y_out = y_in;
            % -- output solution if only one output argument -----------------------------------------------------------
            if(nargout <= 1)
                t_out = y_out;
            end
        end
        
        function resetStats(this)
            this.step_stats.reset();
            this.rhs_stats.reset();
            if(isprop(this, 'nonlinear_solver'))
                this.nonlinear_solver.stats.reset();
                this.nonlinear_solver.linear_solver.stats.reset();
            end
            if(isprop(this, 'linear_solver'))
                this.linear_solver.stats.reset();
            end
            if(isprop(this, 'phi_evaluator'))
                this.phi_evaluator.stats.reset();
            end
        end
        
        function stats_struct = getStats(this)
            stats_struct = struct( ...
                'timestep', copy(this.step_stats), ...
                'rhs', copy(this.rhs_stats) ...
                );
            if(isprop(this, 'nonlinear_solver'))
            	stats_struct.nonlinear = {copy(this.nonlinear_solver.stats), copy(this.nonlinear_solver.linear_solver.stats)};
            end
            if(isprop(this, 'linear_solver'))
                stats_struct.linear = copy(this.linear_solver.stats);
            end
            if(isprop(this, 'phi_evaluator'))
                stats_struct.phi = copy(this.phi_evaluator.stats);
            end
        end
        
    end
    
    methods(Access = protected)
        
        function setStepsize(this, problem)
            this.h = (problem.tspan(end) - problem.tspan(1))/this.num_timesteps;
        end
        
        % == START FUNCTIONS FOR COMPUTING INITIAL INPUTS ==============================================================
        
        function [t_initial, y_initial] = initialConditions(this, problem, scaling_factor)
            %INITIALCONDITIONS Summary of this function goes here
            %   Detailed explanation goes here
            if(nargin == 2)
                scaling_factor = this.h;
            end
            
            % define local timeout function for ode15s
            function status = timeout(start_time, max_duration)
                if(toc(start_time) > max_duration)
                   error('timeout:Interupt', 'Computation exceeded allowed time interval!') 
                end
                status = 0; % function must return 0 for integrator to progress
            end
            
            t_initial   = problem.tspan(1) + scaling_factor * [0; this.starting_times(:)]; % prepent initial condition \tau = 0
            num_outputs = length(t_initial) - 1;
            ode_dim     = length(problem.initial_condition);
            y_initial   = zeros(ode_dim, num_outputs+1); % prenent initial condition
            y_initial(:,1) = problem.initial_condition;
            
            integration_order = startingValueOrder(this, this.starting_times);
            
            for i=1:num_outputs
                % -- start and end integration times -----------
                start_index = integration_order(1,i) + 1;
                end_index   = integration_order(2,i) + 1;
                t_start     = t_initial(start_index);
                y_start     = y_initial(:, start_index);
                t_end       = t_initial(end_index);
                % -- wrap into mutable problem ---------------------
                problem_i = this.initCondProblemWrapper(problem, t_start, t_end, y_start);
                if(diff(problem_i.tspan) == 0)
                    y_end = y_start;
                else
                    if(~isempty(this.starting_integrator))
                        % -- set starting integrator num_timesteps ---
                        est_error = max(1e-16, (this.h ^ this.order)/100); % assume local error of integrator at each step is proprotional to h^order / 100
                        this.starting_integrator.num_timesteps = max(10, ceil(diff(problem_i.tspan)/(est_error ^ (1/this.starting_integrator.order)))); % increased min to 10 steps for ADR2D
                        y_end = this.starting_integrator.solve(problem_i);
                    else % USE ODE15
                        % -- set up timeout function -------------------------------------------------------------------
                        start_time = tic();
                        outputFun = @(t, y, flag) timeout(start_time, this.starting_integrator_timeout);
                        % -- compute initial condition -----------------------------------------------------------------
                        try
                            ode_options = odeset('RelTol', 3e-14, 'AbsTol', 3e-14, 'Jacobian', @(varargin) problem_i.J(varargin{2:end}), 'OutputFcn', outputFun);
                            [~, sol] = ode15s(@(varargin) problem_i.RHS(varargin{2:end}), problem_i.tspan, problem_i.initial_condition, ode_options);
                            y_end = sol(end,:);
                        catch E
                            if strcmp(E.identifier,'timeout:Interupt')
                                warning('Initial condition computation timed out');
                                y_initial(:) = NaN;
                                break;
                            else
                            	rethrow(E);
                            end
                        end
                    end
                end
                y_initial(:, end_index) = y_end;
            end
            % remove \tau = 0 initial condition
            y_initial = y_initial(:,2:end);
            t_initial = t_initial(2:end);
        end
        
        function integration_indices = startingValueOrder(~, times)
            % OPTIMALSTARTINGINTEGRATIONINDICES  computes optimal order to compute output_times. Used by initialConditions
            % PARAMETERS
            %   times (vector) : vector containing output times in local coordinates.
            % RETURNS
            % integration_indices (array) : returns array of size 2 x length(output_times) indicating optimal way to
            %                               compute values times: times(integration_indices(i,2)) should be computed
            %                               using times(integration_indices(i,1)) as initial condition value for i = 1, ...,
            %                               length(output_times). If integration_indices(i,1) = 0, then t = 0 should be
            %                               initial condition.
            
            num_outputs = length(times);
            
            % -- ic solve order --
            z_known = zeros(1, num_outputs + 1); % known z values
            z_known_index = zeros(1, num_outputs + 1); % index of known z values
            
            z_unknown = double(times); %
            z_unknown_index = 1:num_outputs;
            integration_indices = zeros(2, num_outputs);
            
            for i=1:num_outputs
                % -- find closest z_known to any z_unknown ---------
                min_dist = Inf;
                min_start_zindex = -1;
                min_end_zindex = -1;
                min_end_index = -1;
                for j=1:i
                    [min_dist_cand, index_cand] = min(abs(z_known(j) - z_unknown)); % closest z_unknown to z_known(j)
                    if(min_dist_cand < min_dist)
                        min_dist         = min_dist_cand;
                        min_start_zindex = z_known_index(j);
                        min_end_zindex   = z_unknown_index(index_cand);
                        min_end_index    = index_cand;
                    end
                end
                % -- move closest z_unknown to z_known ---
                z_known(i + 1) = times(min_end_zindex);
                z_known_index(i + 1) = min_end_zindex;
                % -- remove unknown point ----------------
                z_unknown(min_end_index) = [];
                z_unknown_index(min_end_index) = [];
                % -- store start and end index -----------
                integration_indices(:,i) = [min_start_zindex, min_end_zindex];
            end
            
        end
        
        function [problem_ic] = initCondProblemWrapper(~, problem, t_start, t_end, y_start)
            %COMPLEXPROBLEMWRAPPER converts a problem with complex time integration bounds into a initial value problem with real
            % integation bounds. This is accomplished by wrapping the original problem into a MutableProblem object.
            
            % -- compute real integration path -------------------------------------------------------------------------
            z           = t_end - t_start;
            r           = abs(z);
            exp_i_theta = z/abs(z);
            tspan       = [0 r];
            % -- wrap problem ------------------------------------------------------------------------------------------
            problem_ic = MutableProblem(struct('tspan', tspan, 'dimension', length(y_start), 'initial_condition', y_start, 'real_valued', problem.real_valued));
            if(problem.has('RHS'))
                problem_ic.RHS = @(varargin) exp_i_theta * problem.RHS(varargin{:});
            end
            if(problem.has('L'))
                problem_ic.L = @(varargin) exp_i_theta * problem.L(varargin{:});
            end
            if(problem.has('N'))
                problem_ic.N = @(varargin) exp_i_theta * problem.N(varargin{:});
            end
            if(problem.has('J'))
                problem_ic.J = @(varargin) exp_i_theta * problem.J(varargin{:});
            end
            if(problem.has('Jxy'))
                problem_ic.Jxy =  @(varargin) exp_i_theta * problem.Jxy(varargin{:});
            end
        end
        
        function setEmptyClassProps(this, props, values)
            for i = 1 : length(props)
                if(isempty(this.(props{i})))
                    this.(props{i}) = values.(props{i});
                end
            end
        end
        
    end
    
end