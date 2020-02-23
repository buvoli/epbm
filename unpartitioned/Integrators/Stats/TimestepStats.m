classdef TimestepStats < StatObject
    
    properties
        num_steps       % total number of newton solves
        total_seconds   % total computation time
        seconds         % array of solve times
    end
    
    properties(Access = protected)
        cell_props = {'seconds'};
        double_props = {'num_steps', 'total_seconds'};
    end
    
    methods
        
        function this = TimestepStats(num_indices)
            if(nargin == 0)
                num_indices = 1;
            end
            this = this@StatObject(num_indices);
        end
        
        function recordStep(this, seconds) % records new timestep with duration seconds
            if(nargin == 1)
                seconds = toc(this.start_time);
            end
            if(this.index == 1)
                this.increment({'num_steps', 'total_seconds'}, {1, sum(seconds)});
            else
                this.increment({'total_seconds'}, {sum(seconds)});
            end
            this.add({'seconds'}, {seconds});
        end
        
        function recordSubstep(this, seconds) % appends seconds to current timestep
            if(nargin == 1)
                seconds = toc(this.start_time);
            end
            this.increment({'total_seconds'}, {sum(seconds)});
            this.incrementLastCell({'seconds'}, {seconds});
        end
        
    end
end