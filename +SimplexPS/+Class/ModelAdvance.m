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
% Public can be set to []
properties
    Ts = 1e-4;          % Sampling period (s)
    Para = [];          % Device parameters            
end

% ### Nontunable publica properties
% Will be read first during construction   
properties(Nontunable)
 
    DeviceType = 0;    % Device type
    PowerFlow = [];     % Power flow parameters
    x0 = [];            % Initial state
    
  	% Discretization methods
    % 1-Forward Euler, 2-Hybrid Euler-Trapezoidal, 2'-General virtual
    % dissipation.
    DiscreMethod = 2;
    
    % Jacobian calculation
    % 1-Initial step, 2-Every step
    LinearizationTimes = 1;
    
    % Direct Feedthrough
    DirectFeedthrough = false;
    
  	% Virtual Resistor
    % 0: do not take out feedthrough as resistor
    % 1: take out feedthrough as resistor
    VirtualResistor = true;
    
    % Initialize to steady-state  
    EquiInitial = true;
    
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
        
    % Void Jacobian
    % The corresponding Jacobian is to set as null (A,C).
    % Equivalently, the marked states are updated by Euler's method.
    % The index is counted backward, i.e. k -> end-k, for compactablity 
    % with Yitong's code
    % This property is only effective when LinearizationTimes == 1
    VoidJacobStates = [0]; %#ok<*NBRAK>
    
    % Void feedthrough
    % void the corresponding feedthrough matrix to break algebraic loop
    % the index is counted backward, i.e. k -> end-k, for compactablity 
    % with Yitong's code
    % This property is only effective when DirectFeedthrough == 0
    VoidFeedthrough = [0];
    
    % Electrical ports
    % The electical port IOs with direct feedthrough (due to Trapezoidal
    % method, for example) can take out the feedthrough parts as virtual
    % resistors to eliminate algebraic loops. 
    % This property is only effective when DirectFeedthrough == 0 and
    % VirtualResistor == 1.
    ElecPortIOs = [1,2];
    
    % Steady-state operating points
   	x_e;        % State
 	u_e;        % Input
    xi;         % Angle difference calculated by power flow analysis
    
    % Dynamic SS Model
    Ak;
    Bk;
    Ck;
    Dk;
    Wk;         % State update gain
    Qk;         % Feedthrough gain
    Gk;         % Virtual resistor gain
    Fk;         % Feedthrough gain by taking out the virtual resistor
    
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

methods
    % constructor
    function obj = ModelAdvance(varargin)
        
        % Support name-value pair arguments when constructing object
        setProperties(obj,nargin,varargin{:});

    end
end

% ### Static methods
methods(Static)

    %% Equilibrium
  	function SetEquilibrium(obj)
        [obj.x_e,obj.u_e,obj.xi] = obj.Equilibrium(obj);
    end
    
  	function [read1,read2,read3,read4] = GetEquilibrium(obj)
        if (isempty(obj.u_e) || isempty(obj.xi))
            % Only check u_e and xi, cause x_e can be empty for certain
            % device such as floating bus.
            error(['Error: The equilibrium is null']);
        else
            y_e = obj.StateSpaceEqu(obj,obj.x_e,obj.u_e,2);
            read1 = obj.x_e;
            read2 = obj.u_e;
            read3 = y_e;
            read4 = obj.xi;
        end
    end
    
    % calc equilibrium x_e and u_e from power flow
    function [x_e, u_e, xi] = Equilibrium(obj)
        error('The Equilibrium method should be overloaded in subclasses.');
        % set x_e,u_e and xi according to PowerFlow
        % x_e and u_e are column vectors with the same dimention as x and u
    end
        
    %% Dynamic SS
    function SetDynamicSS(obj,xk,uk)
        [Ak,Bk,Ck,Dk,Wk,Qk,Gk,Fk] = obj.CalcDynamicSS(obj,xk,uk);
        
        obj.Ak = Ak;
        obj.Bk = Bk;
        obj.Ck = Ck;
        obj.Dk = Dk;
        obj.Wk = Wk;
        obj.Qk = Qk;
        obj.Gk = Gk;
        obj.Fk = Fk;           
    end
    
    function [Ak,Bk,Ck,Dk,Wk,Qk,Gk,Fk] = CalcDynamicSS(obj,xk,uk) 
        
        [Ak,Bk,Ck,Dk] = obj.Linearization(obj,xk,uk);
        
        if obj.LinearizationTimes == 1
            if ~isempty(obj.VoidJacobStates) % null the corresponding column of A
                Ak(:,end-obj.VoidJacobStates) = zeros(length(Ak),length(obj.VoidJacobStates));   
            end

            if ~isempty(obj.VoidJacobStates) % null the corresponding row of B
                Bk(end-obj.VoidJacobStates,:) = zeros(length(obj.VoidJacobStates),length(Bk(1,:)));   
            end

            if ~isempty(obj.VoidJacobStates) % null the corresponding column of C
                Ck(:,end-obj.VoidJacobStates) = zeros(length(Ck(:,1)),length(obj.VoidJacobStates));   
            end
        end
        
        Wk = (eye(length(Ak))/obj.Ts - 1/2*Ak)^(-1);    % state update
        Qk = Dk + 1/2*Ck*Wk*Bk;                         % feedthrough
        
        switch obj.DiscreMethod          
            case 1
                Fk = Dk;                
            case 2                
                Fk = Qk;
            otherwise
                error('Invalid discretization method.');
        end
        
        Gk = zeros(size(Fk));
                
        if ~obj.DirectFeedthrough
            if obj.VirtualResistor
                if ~isempty(obj.ElecPortIOs)                    
                    Gk_ = Fk(obj.ElecPortIOs,obj.ElecPortIOs);
                    Gk_ = diag(Gk_);                            % Get the diagonal elements
                    Gk_ = mean(Gk_);                            % Calculate the average
                    Gk_ = Gk_*eye(length(obj.ElecPortIOs));     % Form the resistor matrix
                    Gk(obj.ElecPortIOs,obj.ElecPortIOs) = Gk_;
                    Fk(obj.ElecPortIOs,obj.ElecPortIOs) = zeros(size(Gk_));
                end                       
            end

            if ~isempty(obj.VoidFeedthrough)                    
                if obj.VoidFeedthrough(1) == -1                 % null all of F
                     Fk = zeros(size(Fk));
                else                                            % null the corresponding row of F
                     Fk(end-obj.VoidFeedthrough,:) = zeros(size(Fk(end-obj.VoidFeedthrough,:)));
                end                     
            end
        end
    end
    
    %% Virtual resistor
    % get the virtual resistor from property obj.Gk
    function Rv = GetVirtualResistor(obj)
        Rv = 1/obj.Gk(1,1)/1.002;
    end
    
    % calculate the virtual resistor from parameter
    function Rv = CalcVirtualResistor(obj)
        [x_e,u_e,~] = obj.Equilibrium(obj);
        [~,~,~,~,~,~,Gk,~] = CalcDynamicSS(obj,x_e,u_e);        
        Rv = 1/Gk(1,1)/1.002;
    end
    
end

% ### Protected default methods provided by matlab
% Notes: The following methods are used for simulink model.
methods(Access = protected)

    % Perform one-time calculations, such as computing constants
    function setupImpl(obj)
        
        % Initialize x_e, u_e
        obj.SetEquilibrium(obj);
        
        % Initialize uk
        if obj.EquiInitial
            obj.uk = obj.u_e;
        else
            obj.uk = 0*obj.u_e;
        end  
        
        % Initialize A, B, C, D
        if obj.LinearizationTimes == 1
            obj.SetDynamicSS(obj,obj.x_e,obj.u_e);            
        else
            if obj.EquiInitial 
                obj.SetDynamicSS(obj,obj.x_e,obj.u_e);
            else
                obj.SetDynamicSS(obj,obj.x0,0*obj.u_e);
            end                          
        end
                
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
                    obj.SetDynamicSS(obj,obj.x,obj.u);
                end
                delta_x = obj.Wk * obj.StateSpaceEqu(obj,obj.x,u,1);               
                obj.x = obj.x + delta_x;             
        end
        
        obj.Timer = obj.Timer + obj.Ts;
    end
        
    % ### Calculate output y
	function y = outputImpl(obj,u)
        
        switch obj.DiscreMethod
            
            % ### Case 1: Forward Euler
            % Notes: Can we also seperate the virtual resistor for forward
            % Euler method?
          	case 1
                if obj.DirectFeedthrough
                    y = obj.StateSpaceEqu(obj,obj.x,u,2);
                else
                    if obj.VirtualResistor
                        y = obj.StateSpaceEqu(obj,obj.x,obj.uk,2) - obj.Gk*obj.uk + obj.Fk*(u-obj.uk);
                    else
                        y = obj.StateSpaceEqu(obj,obj.x,obj.uk,2) + obj.Fk*(u-obj.uk);
                    end
                end
                
            % ### Case 2 : Hybrid Euler-Trapezoidal (Yunjie's Method)
            % ### Case 2': General virtual dissipation (Yitong's Method)
            % case 2' turns out to be equivalent to case 2 and therefore 
            % is combined with case 2
            case 2
                
                if obj.DirectFeedthrough
                    y = obj.StateSpaceEqu(obj,obj.x,u,2) + 1/2*obj.Ck*obj.Wk*obj.StateSpaceEqu(obj,obj.x,u,1);
                else
                    if obj.VirtualResistor
                        y = obj.StateSpaceEqu(obj,obj.x,obj.uk,2) + 1/2*obj.Ck*obj.Wk*obj.StateSpaceEqu(obj,obj.x,obj.uk,1) - obj.Gk*obj.uk + obj.Fk*(u-obj.uk);
                    else
                        y = obj.StateSpaceEqu(obj,obj.x,obj.uk,2) + 1/2*obj.Ck*obj.Wk*obj.StateSpaceEqu(obj,obj.x,obj.uk,1) + obj.Fk*(u-obj.uk);
                    end
                end
        end
        
        obj.uk = u;                 % store the current u=u[k+1/2]
        obj.Start = 1;        
    end
    
  	% Set direct or nondirect feedthrough status of input
    function flag = isInputDirectFeedthroughImpl(obj)
        flag = obj.DirectFeedthrough;
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