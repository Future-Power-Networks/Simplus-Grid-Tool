% This class gives a template for models based on "Class_Model_Advance"

% Author(s): Yitong Li

%% Notes
%
% Double check the index consistency between strings and equations.
%
% Double check the index consistency when getting inputs, states,
% paramters, etc.

%% Class

classdef Class_Model_Template < Class_Model_Advance
    methods(Static)
        
        % Set the strings of state, input, output
        function SetString(obj)
        	obj.StateString  = {};        % x
            obj.InputString  = {};        % u
            obj.OutputString = {};        % y
        end
        
        % Calculate the equilibrium
        function Equilibrium(obj)
         	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get input
            u(1);
            
            % Get state
            x(1);
            
            % Get parameter
            obj.Para(1);
            
            % State space equation
            if CallFlag == 1
                % ### State equation
                Output = f_xu;
            elseif CallFlag == 2
                % ### Output equation
                Output = g_xu;
            end
        end
        
    end
end