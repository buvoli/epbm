classdef KrylovStats < StatObject
    
    properties
        num_projections % total number of projections
        total_seconds   % total computation time
        iterations      % An array of dimension numProjections containing iterations for each projection
        residuals       % An array of residuals
        seconds         % array of solve times
    end
    
    properties(Access = protected)
        cell_props = {'iterations', 'residuals', 'seconds'};
        double_props = {'num_projections', 'total_seconds'};
    end
    
    methods
       
        function this = KrylovStats(num_indices)
            if(nargin == 0)
                num_indices = 1;
            end            
            this = this@StatObject(num_indices);
        end
        
        function recordSolve(this, iterations, residuals, seconds)
            this.increment({'num_projections', 'total_seconds'}, {1, sum(seconds)});
            this.add({'iterations', 'residuals', 'seconds'}, {iterations, residuals, seconds});
        end
        
    end
end