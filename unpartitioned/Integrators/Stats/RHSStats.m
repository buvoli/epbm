classdef RHSStats < StatObject
    
    properties
        num_evaluations % total right hand side evaluations
        total_seconds   % total second for evaluating rhs
        seconds         % array of rhs eval times
    end
    
    properties(Access = protected)
        cell_props = {'seconds'};
        double_props = {'num_evaluations', 'total_seconds'};
    end
    
    methods
        
        function this = RHSStats(num_indices)
            if(nargin == 0)
                num_indices = 1;
            end            
            this = this@StatObject(num_indices);
        end
        
        function recordRHSEval(this, seconds)
            if(nargin == 1)
                seconds = toc(this.start_time);
            end            
            this.increment({'num_evaluations', 'total_seconds'}, {1, sum(seconds)});
            this.add({'seconds'}, {seconds});
        end
        
    end
end