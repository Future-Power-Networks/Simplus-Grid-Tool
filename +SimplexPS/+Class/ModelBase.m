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

methods(Static)
   	% Read properties
  	function [Read1,Read2,Read3] = ReadString(obj)
        Read1 = obj.StateString;
        Read2 = obj.InputString;
        Read3 = obj.OutputString;
    end
	function [Read1,Read2] = ReadSS(obj)
        if obj.ModelBaseType == 1
            MatrixSS = {obj.A,obj.B,obj.C,obj.D};
            ModelSS = ss(obj.A,obj.B,obj.C,obj.D);
            Read1 = MatrixSS;
            Read2 = ModelSS;
        else
            error('Error: The model is in dss rather than ss form.');
        end
    end
    function [Read1,Read2] = ReadDSS(obj)
        if obj.ModelBaseType == 2
          	MatrixDSS = {obj.A,obj.B,obj.C,obj.D,obj.E};
            ModelDSS = dss(obj.A,obj.B,obj.C,obj.D,obj.E);
            Read1 = MatrixDSS;
            Read2 = ModelDSS;
        else
            error('Error: The model is in ss rather than dss form.');
        end
    end
    
    % Write properties
 	function WriteString(obj,StateString,InputString,OutputString)
        obj.StateString  = StateString;     
        obj.InputString  = InputString;     
        obj.OutputString = OutputString;
    end
  	function LoadSS(obj,G)
        % Get the date from G
        A = G.A;
        B = G.B;
        C = G.C;
        D = G.D;
        E = G.E;
        
        % Check if G is in descriptor state space form
        if ( ~isempty(E) )
            error(['Error: the model is in dss rather than ss form.']);
        end
        
        % Set the properties
        obj.A = A;
        obj.B = B;
        obj.C = C;
        obj.D = D;
        obj.ModelBaseType = 1;
    end
  	function LoadDSS(obj,G)
        % Get the date from G
        A = G.A;
        B = G.B;
        C = G.C;
        D = G.D;
        E = G.E;
        
        % Check if G is in descriptor state space form
        if ( isempty(E) && (~isempty(A)) )
            error(['Error: the model is in ss rather than dss form.']);
        end
        
        % Set the properties
        obj.A = A;
        obj.B = B;
        obj.C = C;
        obj.D = D;
        obj.E = E;
        obj.ModelBaseType = 2;
    end
    
 	% Linearize state and output equations to get the linearized state
    % space matrices at a given steady-state operating point
    function Linearization(obj,x_e,u_e)

        % Calculate equilibrium of dx_e and y_e
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

        obj.A = A; obj.B = B; obj.C = C; obj.D = D;
        obj.ModelBaseType = 1;
    end
    
end

end % End class definition

%% Class-related functions
% CAN NOT be accessed by script.
% CAN be called by methods.