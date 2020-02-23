classdef BufferedTimestepStats < StatObject
    
    properties
        num_steps      % total number of newton solves
        total_seconds   % total computation time
        seconds         % array of solve times
    end
    
    properties(Access = protected)
        cell_props = {'seconds'};
        double_props = {'num_steps', 'total_seconds'};
    end
    
    methods
        
        function this = BufferedTimestepStats(num_indices)
            if(nargin == 0)
                num_indices = 1;
            end            
            this = this@StatObject(num_indices);
        end
        
        function recordStep(this, seconds)
            this.increment({'num_steps', 'total_seconds'}, {1, sum(seconds)});
            this.add({'seconds'}, {seconds});
        end
        
    end
end