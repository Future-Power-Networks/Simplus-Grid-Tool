% This class shows a template for models based on "ModelAdvance". 

% Author(s): Yitong Li

%% Notes
%
% This is just a template, for more detailed examples, please see
% "Inductor.m", which is a single-phase inductor, and see
% "SynchronousMachine.m", which is a three-phase synchronous machine.
%
% Double check the index consistency between strings and equations.
%
% Double check the index consistency when getting inputs, states,
% paramters, etc.

%% Class

% Please change "ModelTemplate" here to your customized name and save this
% file as your customized name as well.
classdef ModelTemplate < SimplexPS.Class.ModelAdvance
    
  	methods
        % Constructor
        % Please change "ModelTemplate" to your customized name.
        function obj = ModelTemplate(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        
        % Set the strings of state, input, output
        function [State,Input,Output] = SignalList(obj)
        	StateString  = {'x1','x2'};        % x
            InputString  = {'v'};        % u
            OutputString = {'i'};        % y
        end
        
        % Calculate the equilibrium
        function [x_e,u_e,xi] = Equilibrium(obj)
         	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            
            % Get parameters
            obj.Para(1);
            
            % Calculate equilibrium
            
            % Set equilibrium
            x_e = [];
            u_e = [];
            xi  = [];
        end
        
    	% State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
          	% Get parameter
            obj.Para(1);
            
        	% Get state
            x(1);
            
            % Get input
            u(1);
            
            % State space equations
          	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1
                % ### State equation
                
                f_xu = [];
                Output = f_xu;
            elseif CallFlag == 2
                % ### Output equation
                
                g_xu = [];
                Output = g_xu;
            end
        end
        
    end
end