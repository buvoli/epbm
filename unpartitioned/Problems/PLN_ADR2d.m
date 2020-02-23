% ======================================================================================================================
% Partitioned ADR Equation - Partitioning with linear and nonlinear terms so that
%
%           F^{1} : u_x + u_xx
%           F^{2} : u(u - 1/2)(1 - u)
%
% ======================================================================================================================

classdef PLN_ADR2d < Problem

    properties(SetObservable)
        tspan     = [0, 0.1];
        params    = struct(     ...
            'N',        200,     ...
            'epsilon',  1/100,  ...
            'alpha',    -10,    ...
            'gamma',    100     ...      
        );        
    end
    
	properties(SetAccess = protected)
    	name = 'pADR2d';
        dimension;
        initial_condition;
        real_valued = true;
        description;
    end
    
    properties(Access = private)
    	linear_operator         % stores stores epsilon * (u_xx + u_yy) + delta * (u_x + u_y)
        UxUy_operator           % stores delta * (u_x + u_y)
        UxxUyy_operator         % stores epsilon * (u_xx + u_yy)
    end
        
    methods
        
        function desc = get.description(this)
            desc = sprintf('2D ADR, N = %d^2', this.params.N);
        end
        
        function up = RHS(this, u, part)
            if(nargin == 2)
                part = 0;
            end
            
            switch part
                case 0 % -- full RHS -----------------------------------------------------------------------------------
                    N  = this.params.gamma * u .* (u - 1/2) .* (1 - u);
                    up = this.linear_operator * u + N;
                case 1 % -- linear component (ADVECTION & DIFFUSION TERM) ----------------------------------------------
                    up = this.linear_operator * u;
                case 2 % -- second component (NONLINEAR TERM) ----------------------------------------------------------
                    up  = this.params.gamma * u .* (u - 1/2) .* (1 - u);
            end
        end

        function J = J(this, u, part)
            if(nargin == 2)
                part = 0;
            end
            
            switch part
                case 0 % -- full Jacobian ------------------------------------------------------------------------------
                    temp = this.params.gamma * (3*u - 1/2 - 3*u.^2);
                    J = this.linear_operator + spdiags(temp, 0, this.dimension, this.dimension);
                case 1 % -- linear component (ADVECTION & DIFFUSION TERM) ----------------------------------------------
                    J = this.linear_operator;                    
                case 2 % -- second component (NONLINEAR TERM) ----------------------------------------------------------
                    temp = this.params.gamma * (3*u - 1/2 - 3*u.^2);
                    J = spdiags(temp, 0, this.dimension, this.dimension);
            end
        end

        function N = N(this, u)
            N = this.params.gamma * u .* (u - 1/2) .* (1 - u);
        end
        
        function L = L(this)
            L = this.linear_operator;
        end       
        
    end
    
    methods(Access = protected)
        
        function reset(this, varargin)
            this.setInitialCondition();
            this.setLinearOperator();
        end
        
        function setDimension(this)
            this.dimension = this.params.N ^ 2;
        end
        
        function setInitialCondition(this)
            N = this.params.N;
            h = 1 / (N - 1);
            size = N * N;

            y0 = zeros(size,1);
            for j = 1:N
                y = (j-1) * h;
                for i = 1:N
                    x = (i-1) * h;
                    index = (j-1)*(N) + i;
                    y0(index) = 256 * ((1 - x) * x * (1 - y) * y)^2 + 0.3;
                end
            end
            
            this.initial_condition = y0;            
        end
        
        function setLinearOperator(this)
            epsilon = this.params.epsilon;
            alpha = this.params.alpha;
            N = this.params.N;
            h = 1 / (N - 1);
            
            this.UxxUyy_operator =  epsilon / h^2 * FD_OP.UXX_UYY_N(N);
            this.UxUy_operator   =  alpha / h * FD_OP.UX_UY_N(N);
            this.linear_operator =  this.UxxUyy_operator + this.UxUy_operator;
        end
        
    end
    
end