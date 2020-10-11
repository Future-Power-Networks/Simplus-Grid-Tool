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
        feedin(CountIn_v) = n;
    end
    if  strcmp(InputStr{n},'T_m')
        % Count_v_In/2 can be seen as the index of bus number
        Port_Tm(CountIn_v/2) = n;
    end
end

% Find idq for feedout, and other output for later use
CountOut = 0;
for n = 1:length(OutputStr)
    if ( strcmp(OutputStr{n},'i_d') || strcmp(OutputStr{n},'i_q') ...
      || strcmp(OutputStr{n},'id') || strcmp(OutputStr{n},'iq'))
        Count_i_Out = Count_i_Out+1;
        feedout(Count_i_Out) = n;
    end
    if strcmp(OutputStr{n},'w')
        % Count_i_Out/2 can be seen as the index of bus number
        Port_w(Count_i_Out/2) = n;
    end
end

% Connect
GsysObj = obj_Feedback(GmObj,ZbusObj,feedin,feedout);
obj_CheckDim(GsysObj);

% Output auxiliary data
[~,GsysDSS] = GsysObj.ReadDSS(GsysObj);	% Descriptor state space model
Port_v = feedin;                         % Index of vdq ports in input vector
Port_i = feedout;                        % Index of idq ports in output vector

end