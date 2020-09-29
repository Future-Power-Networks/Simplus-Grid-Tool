% This class defines advanced properties and methods for state space models
% used in script and simulink.

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% This class makes the model fit in both simulation (simulink) and
% theoratical analysis (script).
%
% Subclass of this class contains the specific models of different devices.
% The models should satisfy these conditions:
% - First two inputs are "v_d" and "v_q"
% - First two outputs are "i_d" and "i_q"
% - Second output is "w"
% - Final state is "theta"
% - The final state "theta" SHOULD NOT appear in other state equations
% - The D matrix for system should be 0.
% 
% Available device type:
% 00:   Synchronous generator (SG)
% 10:   PLL-controlled voltage source inverter (VSI)
% 20:   Droop-controlled voltage source inverter (VSI)
% 90:   Single-phase inductor, for test

%% References


%% Class

classdef Class_Model_Advance < Class_Model_Linearization ...
                             & matlab.system.mixin.Nondirect ...
                             & matlab.system.mixin.Propagates
   
%%
% =================================================
% Properties
% =================================================
% ### Public properties
properties
    DeviceType = [];  	% Device type
    Para = [];          % Device parameters
    PowerFlow = [];     % Power flow parameters
    Ts = [];            % Sampling period (s)
    x0 = [];            % Initial state
end

% ### Discrete state
% CAN be modefied by methods
% MUST be numeric or logical variables or fi objects (not string)
% CAN NOT be set with default values
% MUST be initialized before doing simulation
properties(DiscreteState)
    % Notes: It is better to put only x here, and put other states (such as
    % x[k+1], x[k-1], u[k+1]) to other types of properties, because the
    % size of states characteristics are also defined below.
    x;          % It is a column vector generally
end

% ### Protected properties
properties(Access = protected)
 	% Steady-state operating points
   	x_e;
 	u_e;
    xi;         % Angle difference
    
    % Used for Trapezoidal method
    Wk;
    fk;
    xk;
    uk;      
    
end

properties(GetAccess = protected, Constant)
	% Discretization methods: 
    % 1-Forward Euler, 2-Trapezoidal, 3-Virtual Damping
    DiscretizationMethod = 2;
    
    % Linearization times:
    % 1-Initial step, 2-Every step
    LinearizationTimes = 1;
    
    % Damping outside:
    % 0-no damping resistor,1-with damping resistor
    DampingResistor = 1;
end



%%
% =================================================
% Methods
% =================================================
% ### Static methods
methods(Static)
  	function [read1,read2,read3,read4] = ReadEquilibrium(obj)
        y_e = obj.StateSpaceEqu(obj,obj.x_e,obj.u_e,2);
        read1 = obj.x_e;
        read2 = obj.u_e;
        read3 = y_e;
        read4 = obj.xi;
    end
end

% ### Protected default methods provided by matlab
% Notes: The following methods are used for simulink model.
methods(Access = protected)

    % Perform one-time calculations, such as computing constants
    function setupImpl(obj)
        obj.SetString(obj);
        
        % Initialize x_e, u_e
        obj.Equilibrium(obj);
        
        % Initialize A, B, C, D
        obj.Linearization(obj,obj.x_e,obj.u_e);
        
        % Initialize u[k] and x[k]
        obj.uk = obj.u_e;
        obj.xk = obj.x_e;
        
        % Initialize W[k]
        obj.Wk = inv(eye(length(obj.A)) - obj.Ts/2*obj.A);
    end

  	% Update states and calculate output in the same function
    % Notes: This function is replaced by "UpdateImpl" and "outputImpl"
    % function y = stepImpl(obj,u)
    % end
    
    % Update discreate states
    function updateImpl(obj, u)
        
        switch obj.DiscretizationMethod
            
            % ### Forward Euler 
            % s -> Ts/(z-1)
            % => x[k+1] - x[k] = Ts * f(x[k],u[k])
          	case 1
                f_xu = obj.StateSpaceEqu(obj, obj.x, u, 1);
                obj.x = f_xu * obj.Ts + obj.x;
                
            % ### Trapezoidal
            % s -> Ts/2*(z+1)/(z-1)
            % => x[k+1] - x[k] = Ts/2 * (f(x[k+1],u(k+1)) + f(x[k],u[k]))
            % => (x[k+1] - x[k])/Ts =
            % f((x[k+1]+x[k])/2,(u[k+1]+u[k])/2) = f(x[k],u[k]) + Ak*(x[k+1]-x[k])/2 + Bk*(u[k+1]-u[k])/2
            case 2
                
                % Update x[k]
                obj.xk = obj.x;
                
                % Linear Trapezoidal
                if obj.LinearizationTimes == 2
                    obj.Linearization(obj,obj.xk,obj.uk);
                    obj.Wk = inv(eye(length(obj.A)) - obj.Ts/2*obj.A);
                end
                x_k1_Trapez = obj.Wk * (obj.Ts*(obj.StateSpaceEqu(obj,obj.xk,obj.uk,1) + obj.B*(u - obj.uk)/2) ) + obj.x;
                
                % Forward Euler
                x_k1_Euler = obj.StateSpaceEqu(obj, obj.xk, obj.uk, 1)*obj.Ts + obj.xk;
                
                % Split the states
                lx = length(obj.xk);
            	x_k1_linear = x_k1_Trapez(1:(lx-1));
              	x_k1_others = x_k1_Euler(lx:end);   % Theta
         
                % Update x[k+1] and u[k]
                obj.x = [x_k1_linear;
                         x_k1_others];
                obj.uk = u;
                
            % ###  Virtual damping: 
            % Euler -> Trapezoidal
            % s -> s/(1+s*Ts/2)
            % which makes the new state x' replace the old state x with
            % x = x' + Ts/2*dx'/dt
            % Old state space system
            % dx/dt = f(x,u)
            % y     = g(x,u)
            % =>
            % New state space system
            % dx'/dt = f(x'+Ts/2*dx'/dt,u)
            % y      = g(x'+Ts/2*dx'/dt,u)
            case 3
                
                % Linearization
                obj.xk = obj.x;
                if obj.LinearizationTimes == 2
                    obj.Linearization(obj,obj.xk,obj.uk);
                    obj.Wk = inv(eye(length(obj.A)) - obj.Ts/2*obj.A);
                end
                x_k1_VD  = obj.Wk * obj.Ts * obj.StateSpaceEqu(obj,obj.xk,obj.uk,1) + obj.xk;

                % Forward Euler
                x_k1_Euler = obj.Ts * obj.StateSpaceEqu(obj, obj.xk, obj.uk, 1) + obj.xk;
                
                % Split the states
                lx = length(obj.x);
                x_k1_linear = x_k1_VD(1:(lx-1));
                x_k1_others = x_k1_Euler(lx:end);   % Theta
                
                obj.x = [x_k1_linear;
                         x_k1_others];
            otherwise
        end
    end
        
    % Calculate output y
	function y = outputImpl(obj,u)
        switch obj.DiscretizationMethod
          	case 1
                y = obj.StateSpaceEqu(obj, obj.x, u, 2);
            case 2
                if obj.DampingResistor == 1
                    Ck1 = obj.C;
                    Ck2 = obj.C;
                    Ck1(1:2,1:2) = zeros(2,2);
                    Ck2 = Ck2-Ck1;
                    Bk1 = obj.B;
                    Bk2 = obj.B;
                    Bk1(1:2,1:2) = zeros(2,2);
                    Bk2 = Bk2 - Bk1;
                    y = obj.StateSpaceEqu(obj, obj.x, u, 2) ...
                        + obj.Ts/2*(Ck1*obj.Wk*Bk1*(u - obj.uk) - Ck2*obj.Wk*Bk2*obj.uk) ...
                        + obj.Ts*obj.C*obj.Wk*obj.StateSpaceEqu(obj,obj.x,obj.uk,1);
                else
                    y = obj.StateSpaceEqu(obj, obj.x, u, 2) ...
                        + obj.Ts/2*obj.C*obj.Wk*obj.B*(u - obj.uk) ...
                        + obj.Ts*obj.C*obj.Wk*obj.StateSpaceEqu(obj,obj.x,obj.uk,1);
                end
            case 3
                xold_k1_VD = obj.x + obj.Ts/2*obj.Wk * obj.StateSpaceEqu(obj,obj.x,u,1);
                
                % Split the states
                lx = length(obj.x);
                xold_linear = xold_k1_VD(1:(lx-1));
                xold_others = obj.x(lx);
                
                xold = [xold_linear;
                      	xold_others];
                
                y = obj.StateSpaceEqu(obj, xold, u, 2);
            otherwise
        end
    end
    
  	% Set direct or nondirect feedthrough status of input
    function flag = isInputDirectFeedthroughImpl(~)
       flag = false;
    end
    
    % Initialize / reset discrete-state properties
    function resetImpl(obj)
        % Notes: x should be a column vector
        obj.x = obj.x0;
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
        obj.SetString(obj);
        size = [length(obj.OutputString)];
    end
        
    % Set the characteristics of state
    function [size,dataType,complexity] = getDiscreteStateSpecificationImpl(obj, x)
        obj.SetString(obj);
        size = [length(obj.StateString)];
        dataType = 'double';
        complexity = false;
    end
    
end

end     % End class definition