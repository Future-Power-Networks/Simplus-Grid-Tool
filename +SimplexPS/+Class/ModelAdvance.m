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
% - Third output is "w"
% - Final state is "theta"
% - The final state "theta" SHOULD NOT appear in other state equations
% - The D matrix for system should be 0.

%% References


%% Class

classdef ModelAdvance < SimplexPS.Class.ModelBase ...
                     	& matlab.system.mixin.Nondirect ...
                       	& matlab.system.mixin.Propagates
   
%%
% =================================================
% Properties
% =================================================
% ### Public properties
properties
    Para = [];          % Device parameters
    PowerFlow = [];     % Power flow parameters
    Ts = [];            % Sampling period (s)
    x0 = [];            % Initial state
    
  	% Discretization methods
    % 1-Forward Euler, 2-Hybrid Euler-Trapezoidal, 3-General virtual
    % dissipation.
    DiscreMethod = 1;
    
  	% Damping flag
    % 0-no damping resistor, 1-with damping resistor
    DiscreDampingFlag = 0;
    
    % Linearization times
    % 1-Initial step, 2-Every step
    LinearizationTimes = 1;
end

properties(Nontunable)
    DeviceType;         % Device type
    DirectFeedthrough;
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
    xi;         % Angle difference
    
    % Used for Trapezoidal method
    Wk;         % (I-Ts/2*A)
    xk;
    uk;     
    
    % Timer
    Timer = 0;
    
end

properties(GetAccess = protected, Constant)

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
        
        % Initialize Timer
        obj.Timer = 0;
    end

  	% Update states and calculate output in the same function
    % function y = stepImpl(obj,u); end
 	% Notes: This function is replaced by two following functions:
 	% "UpdateImpl" and "outputImpl", and hence is commented out.
    
    % ### Update discreate states
    function updateImpl(obj, u)
        
        obj.Timer = obj.Timer + obj.Ts;
        
        switch obj.DiscreMethod
            
            % ### Case 1: Forward Euler 
            % s -> Ts/(z-1)
            % which leads to
            % x[k+1]-x[k] = Ts * f(x[k],u[k])
            % y[k+1] = g(x[k],u[k]);
          	case 1
                f_xu = obj.StateSpaceEqu(obj, obj.x, u, 1);
                obj.x = obj.Ts * f_xu  + obj.x;
                
            % ### Case 2: Hybrid Euler-Trapezoidal (Yunjie's Method)
            % s -> Ts/2*(z+1)/(z-1)
            % which leads to
            % state equations:
            % (x[k+1] - x[k])/Ts 
            % = (f(x[k+1],u(k+1)) + f(x[k],u[k]))/2
            % or
            % = f((x[k+1]+x[k])/2, (u[k+1]+u[k])/2) 
            % -> 
            % f(x[k],u[k]) + Ak*(x[k+1]-x[k])/2 + Bk*(u[k+1]-u[k])/2
            % output equations:
            % y[k+1] - y[k] = Ts*Ck*Wk* (Bk*(u[k+1]-u[k])/2 + fk)
            case 2
                
                % Update x[k]
                obj.xk = obj.x;
                
                % Linear Trapezoidal
                if obj.LinearizationTimes == 2
                    % Update the linearized system every step during the
                    % simulation.
                    obj.Linearization(obj,obj.xk,obj.uk);
                    obj.Wk = inv(eye(length(obj.A)) - obj.Ts/2*obj.A);
                end
                x_k1_Trapez = obj.Wk * (obj.Ts*(obj.StateSpaceEqu(obj,obj.xk,obj.uk,1) + obj.B*(u - obj.uk)/2) ) + obj.x;
                
                % Forward Euler
                x_k1_Euler = obj.StateSpaceEqu(obj, obj.xk, obj.uk, 1)*obj.Ts + obj.xk;
                
                % Split the states
                lx = length(obj.xk);
            	x_k1_linearized = x_k1_Trapez(1:(lx-1));    % Trapez    
              	x_k1_others     = x_k1_Euler(lx:end);       % Euler for theta only
         
                % Update x[k+1] and u[k]
                obj.x = [x_k1_linearized;
                         x_k1_others];
                obj.uk = u;
                
            % ###  Case 3: General virtual dissipation (Yitong's Method)
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
            % => 
            % Approximation of the new state space system
            % dx'/dt = f(x',u) + Ts/2*A*dx'/dt
            % y      = g(x',u) + Ts/2*C*dx'/dt
            % =>
            % dx'/dt = W*f(x',u)
            % y      = g(x',u) + Ts/2*C*W*f(x',u)
            % =>
            % (x'[k+1]-x'[k])/Ts = Wk*(f(x'[k],u[k]))
            % y'[k+1] = g(x'[k+1],u[k+1]) + Ts/2*C*W*f(x'[k+1],u[k+1]) 
            case 3
                % Update x[k]
                obj.xk = obj.x;
                
                % Linearization
                if obj.LinearizationTimes == 2
                    obj.Linearization(obj,obj.xk,obj.uk);
                    obj.Wk = inv(eye(length(obj.A)) - obj.Ts/2*obj.A);
                end
                x_k1_VD  = obj.Wk * obj.Ts * obj.StateSpaceEqu(obj,obj.xk,obj.uk,1) + obj.xk;

                % Forward Euler
                x_k1_Euler = obj.Ts * obj.StateSpaceEqu(obj, obj.xk, obj.uk, 1) + obj.xk;
                
                % Split the states
                lx = length(obj.x);
                x_k1_linearized = x_k1_VD(1:(lx-1));
                x_k1_others = x_k1_Euler(lx:end);   % Theta
                
                obj.x = [x_k1_linearized;
                         x_k1_others];
            otherwise
        end
    end
        
    % ### Calculate output y
	function y = outputImpl(obj,u)
     	y_Euler = obj.StateSpaceEqu(obj,obj.x,u,2);
    	ly = length(y_Euler);
        switch obj.DiscreMethod
            % ### Case 1: Forward Euler
          	case 1
                y = y_Euler;
                
            % ### Case 2: Hybrid Euler-Trapezoidal (Yunjie's Method)
            case 2
                if obj.DiscreDampingFlag == 0
                    % Consider the virtual dissipation here
                	y_Trapez = y_Euler ...
                        + obj.Ts/2*obj.C*obj.Wk*obj.B*(u - obj.uk) ...
                        + obj.Ts*obj.C*obj.Wk*obj.StateSpaceEqu(obj,obj.xk,obj.uk,1);
                else
                    % Consider the virtual dissipation by using a
                    % resistor in simulink model. Exclude the dissipation
                    % here.
                    Ck1 = obj.C;
                    Ck2 = obj.C;
                    Ck1(1:2,1:2) = zeros(2,2);  % Set corresponding Ck to 0
                    Ck2 = Ck2 - Ck1;
                    Bk1 = obj.B;
                    Bk2 = obj.B;
                    Bk1(1:2,1:2) = zeros(2,2);  % Set corresponding Bk to 0
                    Bk2 = Bk2 - Bk1;
                    y_Trapez = y_Euler ...
                        + obj.Ts/2*(Ck1*obj.Wk*Bk1*(u - obj.uk) - Ck2*obj.Wk*Bk2*obj.uk) ...
                        + obj.Ts*obj.C*obj.Wk*obj.StateSpaceEqu(obj,obj.x,obj.uk,1);
                end
                y = [y_Trapez(1:(ly-1));
                     y_Euler(ly)];
                
            % ### Case 3: General virtual dissipation (Yitong's Method)
            case 3
                if obj.DiscreDampingFlag == 0
                    % No linearization of g(x,u)
                    % xold_VD = obj.x + obj.Ts/2*obj.Wk * obj.StateSpaceEqu(obj,obj.x,u,1); 
                    % xold_VD = obj.x + obj.Ts/2*(obj.x - obj.xk);                          
                    % y_VD = g(xold_VD,u);
                    
                    % Linearization of g(x,u)
                    dx = obj.Wk * obj.StateSpaceEqu(obj,obj.x,u,1);   	% Has algebraic loop
                    % dx = obj.Wk*(obj.A*obj.x + obj.B*u);
                    % dx = obj.x - obj.xk;                            	% No algebraic loop
                    y_VD = y_Euler + obj.Ts/2*obj.C*dx;
                else
                    % Move Ts/2*Ck*Wk*Bk*u[k+1] to the outside of the
                    % system as a damping resistor in the following
                    % equation:
                    % y_VD = y_Euler + obj.Ts/2*obj.C*obj.Wk * (obj.A*obj.x + obj.B*u);
                    Ck1 = obj.C;
                    Ck1(1:2,1:2) = zeros(2,2);
                    Bk1 = obj.B;
                    Bk1(1:2,1:2) = zeros(2,2);
                    y_VD = y_Euler ...
                           + obj.Ts/2*(obj.C*obj.Wk*obj.A*obj.x + Ck1*obj.Wk*Bk1*u);
                end
                y = [y_VD(1:(ly-1));
                     y_Euler(ly)];
            otherwise
                error(['Error: discretization method.'])
        end
    end
    
  	% Set direct or nondirect feedthrough status of input
    function flag = isInputDirectFeedthroughImpl(obj)
        if obj.DirectFeedthrough
            flag = true;
        else
            flag = false;
        end
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