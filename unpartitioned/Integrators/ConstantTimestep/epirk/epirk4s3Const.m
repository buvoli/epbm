classdef epirk4s3Const < ExponentialIntegratorConst
    properties
        graph_line_style = {};
    end
    
    properties (SetAccess = protected)
        order = 4;
        name = 'epirk4s3';
        description = 'epirk4s3'
        starting_times = 0;
    end
    
    % REMOVE ME ONCE KRYLOV PHI IS BETTER
    properties
        max_m = 100;
    end
    
    methods
        function this = epirk4s3Const(options)
            if(nargin == 0)
                options = struct();
            end
            this = this@ExponentialIntegratorConst(options);
        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct, y_in] = initStepStruct(this, t_in, y_in, problem)
            step_struct = struct();
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, final_step)
            
            %  Coefficients for
            %  U_2     = u_n + alpha21*(p211 Phi_1(g21 z)+p212 Phi_2(g21 z)+p213 Phi_3(g21 z))hF
            %  U_3     = u_n + alpha31*(p311 Phi_1(g31 z)+p312 Phi_2(g31 z)+p313 Phi_3(g31 z))hF + alpha32*Psi(g32 z) hR(U2)
            %  u_{n+1} = u_n + Phi_1(z)hF + (b2p3 Phi_3(z) + b2p4 Phi_4(z)) hR(U2) + (b3p3 Phi_3(z) + b3p4 Phi_4(z)) hR(U3)
            
            g21=1/8;
            g31=1/9;
            g32=0;
            
            alpha21=1/8;
            alpha31=1/9;
            alpha32=0;
            
            p211=1;
            p212=0;
            p213=0;
            
            p311=1;
            p312=0;
            p313=0;
            
            b2p3=-1024;
            b2p4=27648;
            
            b3p3=1458;
            b3p4=-34992;
            
            % Setup initial state.
            h = this.h;
            
            % Correctly order the g-coefficients for adaptive-Krylov
            N = length(y_in);
            zeroVec = zeros(N,1);
            if g21 > g31
                gCoeffVec=[g31 g21];
                KryIndex=[2 1];
            else
                gCoeffVec=[g21 g31];
                KryIndex=[1 2];
            end
            
            tol = 1e-12;
            
            
            step_start_time = tic;
            hF = h * problem.RHS(y_in);  % calculate right hand side
            hA = h * problem.J(y_in);  % calculate Jacobian
            
            % -- Stage 1 -------------------------------------------------------------------------------------------
            %temp1 = AdaptiveKrylov_phipm(gCoeffVec, hA, [zeroVec, hF], tol, false, 1, this.max_m);   % 1st Krylov projection :  Used for first and second stage
            [temp1, clean_exit] = this.phi_evaluator.phibs(gCoeffVec, hA, [zeroVec, hF]);
            if(~clean_exit)
                t_out = NaN;
                y_out = NaN;
                return;
            end
            U2 = y_in + alpha21*temp1(:,KryIndex(1));
            hb1 = h * problem.RHS(U2) - hF - hA*(U2 - y_in);     % Calculate residual r(U2)
            
            % -- Stage 2 -------------------------------------------------------------------------------------------
            U3 = y_in + alpha31*temp1(:,KryIndex(2));
            hb2 = h * problem.RHS(U3) - hF - hA*(U3 - y_in);    % Calculate residual r(U3)
            
            % -- Stage 3 -------------------------------------------------------------------------------------------
            %temp3 = AdaptiveKrylov_phipm(1, hA, [zeroVec, hF, zeroVec,b2p3*hb1+b3p3*hb2, b2p4*hb1+b3p4*hb2], tol, false,1, this.max_m);  % 2nd Krylov projection: Computes last stage
            [temp3, clean_exit] = this.phi_evaluator.phibs(1, hA, [zeroVec, hF, zeroVec,b2p3*hb1+b3p3*hb2, b2p4*hb1+b3p4*hb2]);
            if(~clean_exit)
                t_out = NaN;
                y_out = NaN;
                return;
            end
            
            y = y_in + temp3(:,1);
            
            % Advance simulation state.
            t_out = t_in + h;
            y_out = y;
            
            % Update counters.
            this.step_stats.recordStep(toc(step_start_time));
            
        end
        
    end
    
end