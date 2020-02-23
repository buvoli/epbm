classdef MutableProblem < Problem
    
    properties (SetObservable)
        tspan;
        params;        
    end
    
    properties
    	% -- Functions --
        RHS
        Jxy
        J
        L
        N       
    end
    
    properties(SetAccess = protected)
    	dimension;
        name = 'MutableProblem';
        description = 'Mutable Problem';
        initial_condition;
        real_valued;
    end
    
    methods
            function this = MutableProblem(options)
            if(nargin < 1)
                options = struct();
            end
            this.setProps(options);
            this.setDimension();
            end
            
            function setProps(this, prop_struct)                
                fields = fieldnames(prop_struct);
                values = struct2cell(prop_struct);
                
                for i=1:length(fields)
                    if(isprop(this, fields{i}))
                        this.(fields{i}) = values{i};
                    end
                end
            end
    end
    
    methods(Access = protected)
        
        function setDimension(this, varargin)
            this.dimension = length(this.initial_condition);
        end
        
        function reset(varargin)
        end
        
    end

end