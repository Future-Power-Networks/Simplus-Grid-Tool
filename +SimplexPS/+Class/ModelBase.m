% Base class for state space model

% Author(s): Yitong Li

%% Notes
%
% This class defines the basic properties used for subclass models. These
% properties are defined with "protected" attribute so that these
% properties CAN NOT be accessed by Simulink model but CAN be accessed by
% script (by using methods).
%
% Format of the linear state space model
% dx/dt = A*x + B*u
% y     = C*x + D*u
% Format of the linear descriptor state space model
% E*dx/dt = A*x + B*u
% y       = C*x + D*u
%
% Tips:
% - Use methods to be the ports with external script or simulink model
% - Use "obj" to transfer paramerters between methods
% - Capitalize the first character of customized methods
% - "shift+enter" works only within a function, rather than in the whole
% file.

%%
classdef ModelBase < matlab.System
    
%%   
% =================================================
% Properties
% =================================================

% ### Protected properties
properties(Access = protected)
    % State space matrices
  	A;
    B;
    C;
    D;
    E;
    
    % Strings of state, input, output
    StateString;    % StateString = {'x1','x2',...}
    InputString;
    OutputString;
    
    % Model type
    ModelBaseType;  % 1-ss; 2-dss
end

%%
% =================================================
% Methods
% =================================================
% Please use method to output private variables, i.e., "converting private/protected to public".
% CAN invoke/call class-related functions.
% CAN invoke/call external matlab functions which is saved in the visiable matlab path.

% Only the functions begin with "Set" can set/write properties.
% Only the functions begin with "Get" can get/read properties.

methods
    % constructor
    function obj = ModelBase(varargin)
        
        % Support name-value pair arguments when constructing object
        setProperties(obj,nargin,varargin{:});
        
    end
end

% Static methods
% Can be called/invoked without defining a specific instance.
methods(Static)    
  	%% Set/write properties
    function SetString(obj,varargin)
        if nargin > 1
            obj.StateString = varargin{1};
            obj.InputString = varargin{2};
            obj.OutputString = varargin{3};
        else
            % The function "SignalList" is defined in subclass
            [State,Input,Output] = obj.SignalList(obj);
            obj.StateString = State;
            obj.InputString = Input;
            obj.OutputString = Output;
        end
    end
    
    % Set the ss model
  	function SetSS(obj,G)
        % Get the date from G
        A = G.A;
        B = G.B;
        C = G.C;
        D = G.D;
        
        % Check if G is in descriptor state space form
        try E = G.E;
            if ( ~isempty(E) )
                error('Error: the model is in dss rather than ss form.');
            end
        catch
            E = [];
        end
        
        % Set properties
        obj.A = A;
        obj.B = B;
        obj.C = C;
        obj.D = D;
        obj.E = E;
        obj.ModelBaseType = 1;
    end
    
    % Set the dss model
  	function SetDSS(obj,G)
        % Get the date from G
        A = G.A;
        B = G.B;
        C = G.C;
        D = G.D;
        E = G.E;
        
        % Check if G is in descriptor state space form
        if ( isempty(E) && (~isempty(A)) )
            error('Error: the model is in ss rather than dss form.');
        end
        
        % Set properties
        obj.A = A;
        obj.B = B;
        obj.C = C;
        obj.D = D;
        obj.E = E;
        obj.ModelBaseType = 2;
    end
    
    % Set the ss model from linearization function.
    function SetSSLinearized(obj,x_e,u_e)
        % Calculate linearized state space model
        [A,B,C,D] = obj.Linearization(obj,x_e,u_e);
        
        % Set properties
        obj.A = A;
        obj.B = B;
        obj.C = C;
        obj.D = D;
        obj.E = [];
        obj.ModelBaseType = 1;
    end
    
   	%% Get properties  
    function Value = GetProperty(obj,Name)
        Value = obj.(Name);
    end
    
    function [State,Input,Output] = GetString(obj)
        State = obj.StateString;
        Input = obj.InputString;
        Output = obj.OutputString;
    end
    
	function [Read1,Read2] = GetSS(obj)
        if obj.ModelBaseType == 1
            MatrixSS = {obj.A,obj.B,obj.C,obj.D};
            ModelSS = ss(obj.A,obj.B,obj.C,obj.D);
            Read1 = MatrixSS;
            Read2 = ModelSS;
        else
            error('Error: The model is in dss rather than ss form.');
        end
    end
    
    function [Read1,Read2] = GetDSS(obj)
        if obj.ModelBaseType == 2
          	MatrixDSS = {obj.A,obj.B,obj.C,obj.D,obj.E};
            ModelDSS = dss(obj.A,obj.B,obj.C,obj.D,obj.E);
            Read1 = MatrixDSS;
            Read2 = ModelDSS;
        else
            error('Error: The model is in ss rather than dss form.');
        end
    end
    
    %% Linearization
 	% Linearize state and output equations to get the linearized state
    % space matrices at a given steady-state operating point
    function [A,B,C,D] = Linearization(obj,x_e,u_e)

        % Calculate equilibrium of dx_e and y_e
        % The function "StateSpaceEqu" is defined in the subclass
        dx_e = obj.StateSpaceEqu(obj, x_e, u_e, 1);
        y_e  = obj.StateSpaceEqu(obj, x_e, u_e, 2);

        % Calculate length
        lx = length(x_e);
        lu = length(u_e);
        ly = length(y_e);

        % Initialize A,B,C,D
        A = zeros(lx,lx);
        B = zeros(lx,lu);
        C = zeros(ly,lx);
        D = zeros(ly,lu);

        % Get the perturb size
        perturb_factor = 1e-6;

        % Perturb x to calculate Ass and Css
        for i = 1:lx
            x_p = x_e;                      % Reset xp
            perturb_x = perturb_factor * abs(1+abs(x_e(i))); 	% Perturb size
            x_p(i) = x_e(i) + perturb_x;                        % Positive perturb on the ith element of xp
            dx_p = obj.StateSpaceEqu(obj, x_p, u_e, 1);
            y_p  = obj.StateSpaceEqu(obj, x_p, u_e, 2);
            A(:,i) = (dx_p - dx_e)/(x_p(i) - x_e(i));
            C(:,i) = (y_p - y_e)/(x_p(i) - x_e(i));
        end

        % Perturb u to calculate Bss and Dss
        for i = 1:lu
            up = u_e;                       % Reset up
            perturb_u = perturb_factor * abs(1+abs(u_e(i)));    % Perturb size
            up(i) = u_e(i) + perturb_u;                         % Positve perturb on the ith element of up
            dx_p = obj.StateSpaceEqu(obj, x_e, up, 1);
            y_p  = obj.StateSpaceEqu(obj, x_e, up, 2);
            B(:,i) = (dx_p - dx_e)/(up(i) - u_e(i));
            D(:,i) = (y_p - y_e)/(up(i) - u_e(i));
        end 
    end   
    
    %% virtual functions to be overloaded in subclass
    % State space equation for the system
    function rtn = StateSpaceEqu(obj, x, u, flag)
        error('The StateSpaceEqu method should be overloaded in subclasses.');
        % rtn = dx/dt = f(x,u), if flag == 1
        % rtn =  y    = g(x,u), if flag == 2
    end
    
    % Signal list in the state space model
    function [State,Input,Output] = SignalList(obj)
        error('The SignalList method should be overloaded in subclasses.');
        % set signal list by string arrays in the format below:
        % State  = {'x1','x2', ...};
        % Input  = {'u1','u2', ...};
        % Output = {'y1','y2', ...};
    end
end

end % End class definition

%% Class-related functions
% CAN NOT be accessed by script.
% CAN be called by methods.