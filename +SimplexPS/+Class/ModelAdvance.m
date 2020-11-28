% This class defines advanced properties and methods for state space models
% used in script and simulink.

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% This class makes the model fit in both simulation (simulink) and
% theoratical analysis (script).
%
% Subclass of this class contains the specific models of different devices.

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
    Ts = 1e-4;          % Sampling period (s)
    Para = [];          % Device parameters   
end

% Nontunable properties can be set to []
properties(Nontunable)
    
    PowerFlow = [];     % Power flow parameters
    x0 = [];            % Initial state
    
  	% Discretization methods
    % 1-Forward Euler, 2-Hybrid Euler-Trapezoidal, 2'-General virtual
    % dissipation.
    DiscreMethod = 1;
    
    % Jacobian calculation
    % 1-Initial step, 2-Every step
    LinearizationTimes = 1;
    
    % Void Jacobian
    % The corresponding Jacobian is to set as null (A,C).
    % Equivalently, the marked states are updated by Euler's method.
    % The index is counted backward, i.e. k -> end-k, for compactablity 
    % with Yitong's code
    VoidJacobStates = [0];
    
    % Void feedthrough
    % void the corresponding feedthrough matrix to break algebraic loop
    % the index is counted backward, i.e. k -> end-k, for compactablity 
    % with Yitong's code
    VoidFeedthrough = [0];
    
  	% Takeout feedthrough
    % 0: do not take out feedthrough as resistor (may have algebraic loop)
    % 1: take out feedthrough as resistor (eliminate algebraic loop)
    DiscreDampingFlag = 1;
    
    % Electrical ports
    % The electical port IOs with direct feedthrough (due to Trapezoidal
    % method, for example) can take out the feedthrough parts as virtual
    % resistors to eliminate algebraic loops. This property is used with 
    % DiscreDampingFlag to control which ports to take out.
    ElecPortIOs = [1,2];
    
    % Initialize to steady-state  
    EquiInitial = 0;
    
    % Device type
    DeviceType = [];
    
    % Direct Feedthrough
    DirectFeedthrough = 0;
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
    xi;         % Angle difference calculated by power flow analysis
    
    % Used for Trapezoidal method
    Ak;
    Bk;
    Ck;
    Dk;
    Wk;     % State update gain
    Qk;     % Feedthrough gain
    Gk;     % Virtual resistor gain
    Fk;     % Feedthrough gain by taking out the virtual resistor 
    
    % Store previous input. Used for eliminate algebraic loop.
    uk;     
    
    % Timer
    Timer = 0;    
    
    % Start Flag
    Start = 0;
end

properties(GetAccess = protected, Constant)

end

%%
% =================================================
% Methods
% =================================================


% ### Static methods
methods(Static)
    % Read Equilibrium
  	function [read1,read2,read3,read4] = ReadEquilibrium(obj)
        if (isempty(obj.x_e) || isempty(obj.u_e) || isempty(obj.xi))
            error(['Error: The equilibrium is null']);
        else
            y_e = obj.StateSpaceEqu(obj,obj.x_e,obj.u_e,2);
            read1 = obj.x_e;
            read2 = obj.u_e;
            read3 = y_e;
            read4 = obj.xi;
        end
    end
    
    % Preparation for Hybrid Euler-Trapezoidal update
    function PrepareHybridUpdate(obj) 
        
        obj.Ak = obj.A;
        if ~isempty(obj.VoidJacobStates) % null the corresponding column of A
            %obj.Ak(end-obj.VoidJacobStates,:) = zeros(length(obj.VoidJacobStates),length(obj.Ak));
            obj.Ak(:,end-obj.VoidJacobStates) = zeros(length(obj.Ak),length(obj.VoidJacobStates));   
        end
        
        obj.Bk = obj.B;
        if ~isempty(obj.VoidJacobStates) % null the corresponding row of B
           obj.Bk(end-obj.VoidJacobStates,:) = zeros(length(obj.VoidJacobStates),length(obj.Bk(1,:)));   
        end
        
        obj.Ck = obj.C;
        if ~isempty(obj.VoidJacobStates) % null the corresponding column of C
            obj.Ck(:,end-obj.VoidJacobStates) = zeros(length(obj.Ck(:,1)),length(obj.VoidJacobStates));   
        end
        
        obj.Dk = obj.D;
        
        obj.Wk = (eye(length(obj.Ak))/obj.Ts - 1/2*obj.Ak)^(-1);    % state update
        obj.Qk = obj.Dk + 1/2*obj.Ck*obj.Wk*obj.Bk;                 % feedthrough
        if ~isempty(obj.VoidFeedthrough) % null the corresponding row of Q ???
            obj.Qk(end-obj.VoidFeedthrough,:) = zeros(size(obj.Qk(end-obj.VoidFeedthrough,:)));   
        end
                
        % Split obj.Qk = obj.Gk + obj.Fk
        % where Gk is the virtual resistor gain, and Fk is the feedthrough
        % gain by taking out the virtual resistor
        obj.Gk = zeros(size(obj.Qk));
        obj.Fk = obj.Qk;
        if ~isempty(obj.ElecPortIOs)                    
            Gk_ = obj.Qk(obj.ElecPortIOs,obj.ElecPortIOs);
            Gk_ = diag(Gk_);                        % Get the diagonal elements
            Gk_ = mean(Gk_);                        % Calculate the average
            Gk_ = Gk_*eye(length(obj.ElecPortIOs)); % Form the resistor matrix
            obj.Gk(obj.ElecPortIOs,obj.ElecPortIOs) = Gk_;
            obj.Fk(obj.ElecPortIOs,obj.ElecPortIOs) = zeros(size(Gk_));
        end
    end
    
    % Caculate the virtual resistor
    function Rv = GetVirtualResistor(obj)
        if ~isempty(obj.Gk)
            Gk_ = obj.Gk(obj.ElecPortIOs,obj.ElecPortIOs);
            Rv = (Gk_(1,1))^(-1);
        else
            error(['Error: Gk is empty.'])
        end
    end
    
    % State space equation for the system
    function rtn = StateSpaceEqu(obj,x,u,flag)
        rtn = [];
        % should be overrided in the subclass in the following format
        % rtn = dx/dt = f(x,u), if flag == 1
        % rtn =  y    = g(x,u), if flag == 2
    end
    
    % Get equilibrium x_e and u_e from power flow
    function Equilibrium(obj)
        % should be overrided in the subclass in the following format
        % set x_e,u_e and xi according to PowerFlow
        % x_e and u_e are column vectors with the same dimention as x and u
    end
end

% ### Protected default methods provided by matlab
% Notes: The following methods are used for simulink model.
methods(Access = protected)

    % Perform one-time calculations, such as computing constants
    function setupImpl(obj)
        % Set string vectors for x, u, y
        obj.SetString(obj);
        
        % Initialize x_e, u_e
        obj.Equilibrium(obj);
                
        % Initialize A, B, C, D
        obj.Linearization(obj,obj.x_e,obj.u_e);              
        obj.PrepareHybridUpdate(obj);
        
        % Initialize uk
        obj.uk = obj.u_e;
                
        % Initialize Timer
        obj.Timer = 0;
        
        % Initialize start flag
        obj.Start = 0;       
    end
    
    % Initialize / reset discrete-state properties
    function resetImpl(obj)
        % Notes: x should be a column vector
        if obj.EquiInitial
            obj.x = obj.x_e;
        else
            obj.x = obj.x0;
        end
    end

  	% Update states and calculate output in the same function
    % function y = stepImpl(obj,u); end
 	% Notes: This function is replaced by two following functions:
 	% "UpdateImpl" and "outputImpl", and hence is commented out, to avoid
 	% repeated update in algebraic loops
    
    % ### Update discreate states
    function updateImpl(obj, u)
        
        switch obj.DiscreMethod
            
            % ### Case 1: Forward Euler 
            % s -> (z-1)/Ts
            % which leads to
            % x[k+1]-x[k] = Ts * f(x[k],u[k])
            % y[k+1] = g(x[k],u[k]);
          	case 1
                delta_x = obj.Ts * obj.StateSpaceEqu(obj, obj.x, u, 1);
                obj.x = delta_x + obj.x;
                
            % ### Case 2 : Hybrid Euler-Trapezoidal (Yunjie's Method)
            % s -> 2/Ts*(z-1)/(z+1)
            % which leads to
            %
            % state equations:
            % dx[k]/Ts = (x[k+1] - x[k])/Ts = f(x[k+1/2], u[k+1/2])
            %          = f(x[k], u[k+1/2]) + 1/2*A[k]*dx[k] ->
            % dx[k] = W[k]*f(x[k],u[k+1/2]), W[k]=[I/Ts - 1/2*A[k]]^(-1)
            % where x[k+1/2] = (x[k]+x[k+1])/2
            %
            % output equations:
            % y[k+1/2] = g(x[k+1/2],u[k+1/2]) 
            %          = g(x[k],u[k+1/2]) + 1/2*C[k]*dx[k]
            %          = g(x[k],u[k+1/2]) + 1/2*C[k]*W[k]*f(x[k],u[k+1/2])
            % eliminate algebraic loops:
            % y[k+1/2] = g(x[k],u[k-1/2]) + 1/2*C[k]*W[k]*f(x[k],u[k-1/2])
            %            +(D[k] + 1/2*C[k]*W[k]*B[k])*du[k-1/2]
            %          = ... + Q[k]*du[k-1/2], Q[k]=D[k]+1/2*C[k]*W[k]*B[k]
            % 1) Set Q[k]=0, eliminate algebraic loops with delay
            % 2) Q[k]*du[k-1/2] = Q[k]*u[k+1/2] - Q[k]*u[k-1/2])
            %                     ------------- take out as resistor
            %
            % ###  Case 2': General virtual dissipation (Yitong's Method)
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
            % dx'/dt = W[k]*f(x',u)/Ts, W[k] same as case 2
            % y      = g(x',u) + 1/2*C*W[k]*f(x',u)/Ts
            % =>
            % dx'[k] = (x'[k+1]-x'[k]) = W[k]*(f(x'[k],u[k]))
            % y[k]   = g(x'[k],u[k]) + 1/2*C[k]*W[k]*f(x'[k],u[k]) 
            %
            % case 2' turns out to be equivalent to case 2 and therefore 
            % is combined with case 2  
            case 2               
                if obj.LinearizationTimes == 2 
                    % update Jacobian in real time    
                    obj.Linearization(obj,obj.x,u);
                    obj.PrepareHybridUpdate(obj);
                end
                delta_x = obj.Wk * obj.StateSpaceEqu(obj,obj.x,u,1);               
                obj.x = obj.x + delta_x;    % update x             
                
            otherwise
                error(['Error: discretization method invalid.']);
        end
        
        obj.Timer = obj.Timer + obj.Ts;
    end
        
    % ### Calculate output y
	function y = outputImpl(obj,u)
        
        switch obj.DiscreMethod
            
            % ### Case 1: Forward Euler
          	case 1
                y = obj.StateSpaceEqu(obj,obj.x,u,2);
                
            % ### Case 2 : Hybrid Euler-Trapezoidal (Yunjie's Method)
            % ### Case 2': General virtual dissipation (Yitong's Method)
            % case 2' turns out to be equivalent to case 2 and therefore 
            % is combined with case 2
            case 2
                if obj.DiscreDampingFlag
                    y = obj.StateSpaceEqu(obj,obj.x,obj.uk,2) + 1/2*obj.Ck*obj.Wk*obj.StateSpaceEqu(obj,obj.x,obj.uk,1) - obj.Gk*obj.uk + obj.Fk*(u-obj.uk);
                else
                    y = obj.StateSpaceEqu(obj,obj.x,obj.uk,2) + 1/2*obj.Ck*obj.Wk*obj.StateSpaceEqu(obj,obj.x,obj.uk,1) + obj.Qk*(u-obj.uk);
                end
                
            otherwise
                error(['Error: discretization method invalid.']);
        end
        
        obj.uk = u;                 % store the current u=u[k+1/2]
        obj.Start = 1;        
    end
    
  	% Set direct or nondirect feedthrough status of input
    function flag = isInputDirectFeedthroughImpl(obj)
        if obj.DirectFeedthrough
            flag = true;
        else
            flag = false;
        end
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