
%

%   Copyright 2013-2014 The MathWorks, Inc.

classdef MexFunction < handle
    properties(SetAccess = private)
        name = '';
    end
    
    methods
        function reactToBreakpointChange(obj, bp)
            obj.updateBreakpointBitVector(bp.ownerUdd, {});
        end
        
        function reactToBreakpointRemoval(obj, bp)
            obj.updateBreakpointBitVector(bp.ownerUdd, {bp});            
        end
        
        function obj = MexFunction(name)
            obj.name = [name '_sfun'];
        end
    end
    
    methods(Access = private)
        function updateBreakpointBitVector(obj, udd, listOfBpsToIgnore)
            assert(~isempty(udd));
            
            if ~ishandle(udd) || shouldIgnoreBecauseObjectCommented(udd) 
                % do nothing as we should not honor any breakpoints set on
                % commented objects, or closed models             
                return;
            end
            
            machineId = udd.Machine.Id;
            chartId = sfprivate('getChartOf', udd.Id);
            objectType = sf('get', udd.Id, '.isa');
            
            breakpointBitVector = Stateflow.Debug.get_breakpoints_bit_vector_for(udd, listOfBpsToIgnore);     
                        
            objectNumber = getObjectNumberFromUdd(udd);
            
            feval(obj.name,'sf_debug_api','set_breakpoints',...
                machineId,chartId,objectType,objectNumber,breakpointBitVector);
        end
    end
end

function objectNumber = getObjectNumberFromUdd(udd)
    SUB_TYPE = 1; % ugly
    SUPER = 2;
    if isa(udd, 'Stateflow.Transition') && sf('get', udd.Id, '.type') == SUB_TYPE
        % In the case of subwires, we need to get the corresponding
        % superwire because that is where breakpoints ought to live
        superTransId = sf('get', udd.Id, '.subLink.parent');
        assert(sf('get', superTransId, '.type') == SUPER);
        
        objectNumber = sf('get', superTransId, '.number');
    else
        objectNumber = sf('get', udd.Id, '.number');
    end
end

function result = shouldIgnoreBecauseObjectCommented(udd)    
    switch class(udd)
        case {'Stateflow.Chart', 'Stateflow.StateTransitionTableChart'}
            % NOTE: We do not need to worry about breakpoints set on
            % charts/STTs because if the entire chart is commented out
            % then it will not partake in debug initialization
            % (sf_debug_initialize_chart) which eventually calls into
            % here.
            result = false;
        case 'Stateflow.Event'
            % We cannot comment out events
            result = false;
        otherwise            
            result = udd.isCommented;
    end
end
