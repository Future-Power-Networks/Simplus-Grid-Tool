% Author(s): Yitong Li

%% Class

classdef ModelAdvanceNetwork < SimplusGT.Class.ModelBase ...
                     	& matlab.system.mixin.Nondirect ...
                       	& matlab.system.mixin.Propagates
   
%%
% =================================================
% Properties
% =================================================
% ### Public properties 
% Public can be set to []
properties
    Ts;          % Sampling period (s)            
end

% ### Nontunable publica properties
% Will be read first during construction   
properties(Nontunable)
    
end

% ### Discrete state
% CAN be modefied by methods
% MUST be numeric or logical variables or fi objects (not string)
% CAN NOT be set with default values
% MUST be initialized before doing simulation
properties(DiscreteState)
    % Notes: It is better to put only x here, and put other states (such as
    % x[k+1], x[k-1], u[k+1]) to other types of properties, because the
    % size of states characteristics also needs to be defined. Putting
    % other states here would make it confused.
    x;          % It is a column vector generally
end

% ### Protected properties
properties(Access = protected)
        
    % Steady-state operating points
   	x_e;        % State
 	u_e;        % Input  
    
    % Timer
    Timer = 0;    
    
end

properties(GetAccess = protected, Constant)

end

%%
% =================================================
% Methods
% =================================================

methods
    % constructor
    function obj = ModelAdvance(varargin)
        
        % Support name-value pair arguments when constructing object
        setProperties(obj,nargin,varargin{:});

    end
end

% ### Protected default methods provided by matlab
% Notes: The following methods are used for simulink model.
methods(Access = protected)

    % Perform one-time calculations, such as computing constants
    function setupImpl(obj)
         
      	% Initialize x_e, u_e
        [obj.x_e,obj.u_e] = obj.Equilibrium(obj);
        
        % Initialize Timer
        obj.Timer = 0;
           
    end
    
    % Initialize / reset discrete-state properties
    function resetImpl(obj)
        % Notes: x should be a column vector
        obj.x = obj.x_e;
    end
    
    % ### Update discreate states
    function updateImpl(obj, u)
        
       	% ### Forward Euler 
    	% s -> (z-1)/Ts
    	% which leads to
     	% x[k+1]-x[k] = Ts * f(x[k],u[k])
     	% y[k+1] = g(x[k],u[k]);
        delta_x = obj.Ts * obj.StateSpaceEqu(obj, obj.x, u, 1);
     	obj.x = delta_x + obj.x;

        obj.Timer = obj.Timer + obj.Ts;
    end
        
    % ### Calculate output y
	function y = outputImpl(obj,u)
        
     	% ### Forward Euler
     	y = obj.StateSpaceEqu(obj,obj.x,u,2);
    end
    
  	% Set direct or nondirect feedthrough status of input
    function flag = isInputDirectFeedthroughImpl(obj,u)
        flag = 0;
    end
    
    % Release resources, such as file handles
    function releaseImpl(obj)
    end

    % Define total number of inputs for system
    function num = getNumInputsImpl(obj)
        num = 1;
    end

    % Define total number of outputs for system
    function num = getNumOutputsImpl(obj)
        num = 1;
    end
    
    % Set the size of output
    function [size] = getOutputSizeImpl(obj)
        [~,~,Output] = obj.SignalList(obj);
        size = length(Output);
    end
        
    % Set the characteristics of state
    function [size,dataType,complexity] = getDiscreteStateSpecificationImpl(obj, x)        
        [State,~,~] = obj.SignalList(obj);
        size = length(State);
        dataType = 'double';
        complexity = false;
    end
    
end

end     % End class definition