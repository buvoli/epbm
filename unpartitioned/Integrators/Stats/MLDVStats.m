classdef MLDVStats < StatObject
    
    properties
        num_solves      % total number of newton solves
        total_seconds   % total computation time
        seconds         % array of solve times
    end
    
    properties(Access = protected)
        cell_props = {'seconds'};
        double_props = {'num_solves', 'total_seconds'};
    end
    
    methods
        function this = MLDVStats(num_indices)
            if(nargin == 0)
                num_indices = 1;
            end
            this = this@StatObject(num_indices);
        end
        
        function recordSolve(this, seconds)
            this.increment({'num_solves', 'total_seconds'}, {1, sum(seconds)});
            this.add({'seconds'}, {seconds});
        end
        
    end
end