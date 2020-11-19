% This class inherits the base class, and adds an function for writting the
% model properties.

% Author(s): Yitong Li

%% 
% Notes:
% 
% This class CAN be used in script, but CAN NOT be used in simulations.

%% 
classdef Class_Model_DSS < Class_Model_Base

methods(Static)
    % Write properties
    function LoadModel(obj,G)
        % Get the date from G
        A = G.A;
        B = G.B;
        C = G.C;
        D = G.D;
        E = G.E;
        
        % Check if G is in descriptor state space form
        if ( isempty(E) && (~isempty(A)) )
            error(['Error: the system is not in descriptor-state-space form']);
        end
        
        % Set the properties
        obj.A = A;
        obj.B = B;
        obj.C = C;
        obj.D = D;
        obj.E = E;
        obj.MatrixDSS   = {A,B,C,D,E};
        obj.ModelDSS    = dss(A,B,C,D,E);
    end
end

end     % End for class def
