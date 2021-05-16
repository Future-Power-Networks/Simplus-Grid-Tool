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
classdef Class_Model_Base < matlab.System
    
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
    
    % Models
    MatrixSS;   % MatrixSS = {A,B,C,D}
    ModelSS;    % ModelSS = ss(A,B,C,D)
 	MatrixDSS;
    ModelDSS;
end

%%
% =================================================
% Methods
% =================================================
% Please use method to output private variables, i.e., "converting private/protected to public".
% CAN invoke/call class-related functions.
% CAN invoke/call external matlab functions which is saved in the visiable matlab path.

methods(Static)
    % Write properties
 	function WriteString(obj,StateString,InputString,OutputString)
        obj.StateString  = StateString;     
        obj.InputString  = InputString;     
        obj.OutputString = OutputString;
    end
    function WriteSS(obj,MatrixSS,ModelSS)
        obj.MatrixSS = MatrixSS;            
        obj.ModelSS = ModelSS;
    end
    function WriteDSS(obj,MatrixDSS,ModelDSS)
        obj.MatrixDSS = MatrixDSS;
        obj.ModelDSS = ModelDSS;
    end
    
    % Read properties
  	function [Read1,Read2,Read3] = ReadString(obj)
        Read1 = obj.StateString;
        Read2 = obj.InputString;
        Read3 = obj.OutputString;
    end
	function [Read1,Read2] = ReadSS(obj)
        Read1 = obj.MatrixSS;
        Read2 = obj.ModelSS;
    end
    function [Read1,Read2] = ReadDSS(obj)
        Read1 = obj.MatrixDSS;
        Read2 = obj.ModelDSS;
    end
    
    % Construct models
    function ConstructSS(obj)
        obj.ModelSS = ss(obj.A,obj.B,obj.C,obj.D);
    end
    function ConstructDSS(obj)
        obj.ModelDSS = dss(obj.A,obj.B,obj.C,obj.D,obj.E);
    end
end

end % End class definition

%% Class-related functions
% CAN NOT be accessed by script.
% CAN be called by methods.