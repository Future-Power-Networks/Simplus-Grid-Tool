% This class shows a template for models based on "ModelAdvance"

% Author(s): Yitong Li

%% Notes
%
% Double check the index consistency between strings and equations.
%
% Double check the index consistency when getting inputs, states,
% paramters, etc.

%% Class

classdef ModelTemplate < SimplexPS.Class.ModelAdvance
    
    methods(Static)
        
        % Set the strings of state, input, output
        function SignalList(obj)
        	obj.StateString  = {''};        % x
            obj.InputString  = {''};        % u
            obj.OutputString = {''};        % y
        end
        
        % Calculate the equilibrium
        function Equilibrium(obj)
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
            obj.x_e = [];
            obj.u_e = [];
            obj.xi  = [];
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
          	% Get parameter
            obj.Para(1);
            
        	% Get state
            x(1);
            
            % Get input
            u(1);
            
            % State space equation
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