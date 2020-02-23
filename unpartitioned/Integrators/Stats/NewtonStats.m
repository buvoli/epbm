classdef NewtonStats < StatObject
    
    properties
        num_solves      % total number of newton solves
        total_seconds   % total computation time
        iterations      % A cell array of length num_solves containing iterations for each projection
        residuals       % A cell array of residuals
        deltas          % A cell array of deltas
        seconds         % array of solve times
    end
    
    properties(Access = protected)
        cell_props = {'iterations', 'residuals',  'deltas', 'seconds'};
        double_props = {'num_solves', 'total_seconds'};
    end
    
    methods
       
        function this = NewtonStats(num_indices)
            if(nargin == 0)
                num_indices = 1;
            end            
            this = this@StatObject(num_indices);
        end
        
        function addSolve(this, iterations, residuals, deltas, seconds)
            this.increment({'num_solves', 'total_seconds'}, {1, sum(seconds)});
            this.add({'iterations', 'residuals', 'deltas', 'seconds'}, {iterations, residuals, deltas, seconds});
        end
        
    end
end