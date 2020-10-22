% This function connects the admittance model of devices and the impedance
% model of network lines.

% Author(s): Yitong Li

function [GsysObj,GsysDSS,Port_v_dq,Port_i_dq,Port_w,Port_T_m,Port_ang_r] = GmZbus_Connect(GmObj,ZbusObj)

% Get the strings
[~,InputStr,OutputStr] = GmObj.ReadString(GmObj);

% Initialize port
Port_v_dq  = [];
Port_i_dq  = [];
Port_w  = [];
Port_T_m = [];
Port_ang_r = [];

% Find vdq for feedin, and other port for later use
CountIn_v = 0;
for n = 1:length(InputStr)
    if ( strcmp(InputStr{n},'v_d') || strcmp(InputStr{n},'v_q') ...
      || strcmp(InputStr{n},'vd') || strcmp(InputStr{n},'vq'))
        CountIn_v = CountIn_v+1;
        Port_v_dq(CountIn_v) = n;      % Index of port v
    end
    % Count_v_In/2 can be seen as the index of bus number
    IndexIn_Bus = CountIn_v/2;
    if  ( strcmp(InputStr{n},'T_m') || strcmp(InputStr{n},'Tm') )
        Port_T_m(IndexIn_Bus) = n;   % Index of port Tm
    end
    if strcmp(InputStr{n},'ang_r')
        Port_ang_r(IndexIn_Bus) = n;
    end
end

% Find idq for feedout, and other output for later use
CountOut_i = 0;
for n = 1:length(OutputStr)
    if ( strcmp(OutputStr{n},'i_d') || strcmp(OutputStr{n},'i_q') ...
      || strcmp(OutputStr{n},'id') || strcmp(OutputStr{n},'iq'))
        CountOut_i = CountOut_i+1;
        Port_i_dq(CountOut_i) = n;     % Index of port i
    end
    % Count_i_Out/2 can be seen as the index of bus number
    IndexOut_Bus = CountOut_i/2;
    if ( strcmp(OutputStr{n},'w') || strcmp(OutputStr{n},'omega') )
        Port_w(IndexOut_Bus) = n;   % Index of port omega
    end
end

% Connect
GsysObj = obj_Feedback(GmObj,ZbusObj,Port_v_dq,Port_i_dq);
obj_CheckDim(GsysObj);

% Output auxiliary data
[~,GsysDSS] = GsysObj.ReadDSS(GsysObj);	% Descriptor state space model

end