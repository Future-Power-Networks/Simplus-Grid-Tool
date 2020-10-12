% This function connects the admittance model of devices and the impedance
% model of network lines.

% Author(s): Yitong Li

function [GsysObj,GsysDSS,Port_v,Port_i,Port_w,Port_Tm] = GmZbus_Connect(GmObj,ZbusObj)

% Get the strings
[~,InputStr,OutputStr] = GmObj.ReadString(GmObj);

% Find vdq for feedin, and other port for later use
CountIn_v = 0;
for n = 1:length(InputStr)
    if ( strcmp(InputStr{n},'v_d') || strcmp(InputStr{n},'v_q') ...
      || strcmp(InputStr{n},'vd') || strcmp(InputStr{n},'vq'))
        CountIn_v = CountIn_v+1;
        Port_v(CountIn_v) = n;      % Index of port v
    end
    % Count_v_In/2 can be seen as the index of bus number
    IndexIn_Bus = CountIn_v/2;
    if  ( strcmp(InputStr{n},'T_m') || strcmp(InputStr{n},'Tm') )
        Port_Tm(IndexIn_Bus) = n;   % Index of port Tm
    end
end

% Find idq for feedout, and other output for later use
CountOut_i = 0;
for n = 1:length(OutputStr)
    if ( strcmp(OutputStr{n},'i_d') || strcmp(OutputStr{n},'i_q') ...
      || strcmp(OutputStr{n},'id') || strcmp(OutputStr{n},'iq'))
        CountOut_i = CountOut_i+1;
        Port_i(CountOut_i) = n;     % Index of port i
    end
    % Count_i_Out/2 can be seen as the index of bus number
    IndexOut_Bus = CountOut_i/2;
    if ( strcmp(OutputStr{n},'w') || strcmp(OutputStr{n},'omega') )
        Port_w(IndexOut_Bus) = n;   % Index of port omega
    end
end

% Connect
GsysObj = obj_Feedback(GmObj,ZbusObj,Port_v,Port_i);
obj_CheckDim(GsysObj);

% Output auxiliary data
[~,GsysDSS] = GsysObj.ReadDSS(GsysObj);	% Descriptor state space model

end