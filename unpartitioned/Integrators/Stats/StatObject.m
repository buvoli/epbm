classdef StatObject < matlab.mixin.Copyable % (See https://www.mathworks.com/help/matlab/ref/matlab.mixin.copyable-class.html)
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        record = false;
    end    
    
    properties(Access = protected)
        index = 1;
        start_time = 0;
        num_indices = 1;
    end
    
    properties(Abstract = true, Access = protected)
    	cell_props;
        double_props;
    end
    
    methods
        function this = StatObject(num_indices)
            if(nargin < 1)
                num_indices = 1;
            end            
            this.reset(num_indices);
        end
        
        function setIndex(this, index)
            this.index = min(max(0,index),this.num_indices);
        end
        
        function startTimer(this, index)
            if(nargin > 1)
                this.index = min(max(0,index),this.num_indices);
            end
            this.start_time = tic;
        end
        
        function reset(this, num_indices)
            if(nargin == 2)
                this.num_indices = num_indices;
            end            
            this.resetCellProps(this.cell_props);  
            this.resetDoubleProps(this.double_props);  
        end
    end
    
    methods(Access = protected)
        function add(this, fields, values)
            if(this.record)
                if(this.num_indices == 1)
                    for i=1:length(fields)
                        this.(fields{i}){end+1} = values{i};
                    end   
                else
                    for i=1:length(fields)
                        this.(fields{i}){this.index}{end+1} = values{i};
                    end 
                end            
            end
        end
        
        function incrementLastCell(this, fields, values)
            if(this.record)
                if(this.num_indices == 1)
                    for i=1:length(fields)
                        this.(fields{i}){end} = this.(fields{i}){end} + values{i};
                    end   
                else
                    for i=1:length(fields)
                        this.(fields{i}){this.index}{end} = this.(fields{i}){this.index}{end} + values{i};
                    end 
                end            
            end
        end
            
        function increment(this, fields, values)
            if(this.record)
                for i=1:length(fields)
                	this.(fields{i}) = this.(fields{i}) + values{i};
                end           
            end
        end
        
        function resetCellProps(this, props)
        	if(this.num_indices == 1)
            	for i=1:length(props)
                    this.(props{i}) = {};
                end
            else
            	for i=1:length(props)
                    this.(props{i}) = cell(this.num_indices, 1);
                end
            end
        end
        
        function resetDoubleProps(this, props)
        	for i=1:length(props)
            	this.(props{i}) = 0;
            end
        end
        
    end

end

