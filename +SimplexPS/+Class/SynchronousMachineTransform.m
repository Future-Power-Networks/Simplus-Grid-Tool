% This class defines the model of synchronous machine

% Author(s): Yitong Li, Yunjie Gu

%% Notes
%
% The model is in load convention, admittance form.

%% Class

classdef SynchronousMachineTransform < SimplexPS.Class.SynchronousMachine ...
                                     & SimplexPS.Class.ModelAdvanceTransform 
                                     % Diamond shape inheritance is not supported in code generation mode
    
    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
            State	 = {'i_d','i_q','w','theta'};            
            Input	 = {'v_a','v_b','v_c','T_m','v_ex'};
            Output   = {'v_d','v_q','i_a','i_b','i_c','i_d','i_q','w','i_ex','theta'};
        end
        
    end

end     % End class definition