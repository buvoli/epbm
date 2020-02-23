classdef ExponentialIntegratorConst < IntegratorConst
    
    properties
        phi_evaluator
    end
    
    methods
        
        function this = ExponentialIntegratorConst(options)
            % -- set basic integrator properties -----------------------------------------------------------------------
            this@IntegratorConst(options);
            % -- set local properties ----------------------------------------------------------------------------------
            default_field_value_pairs = {
                {'phi_evaluator', KIOPS()}
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            % -- set props ---------------------------------------------------------------------------------------------
            props = {'phi_evaluator'};
            this.setEmptyClassProps(props, options);
        end
        
    end
    
end