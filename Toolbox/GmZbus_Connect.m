% This function connects the admittance model of devices and the impedance
% model of network lines.

% Author(s): Yitong Li

function [GsysObj,GsysDSS,vport,iport] = GmZbus_Connect(GmObj,ZbusObj)

% Get the strings
[~,InputStr,OutputStr] = GmObj.ReadString(GmObj);

% Find vdq for feedin
CountIn = 0;
for n = 1:length(InputStr)
    if ( strcmp(InputStr{n},'v_d') || strcmp(InputStr{n},'v_q') ...
      || strcmp(InputStr{n},'vd') || strcmp(InputStr{n},'vq'))
        CountIn = CountIn+1;
        feedin(CountIn) = n;
    end
end

% Find idq for feedout
CountOut = 0;
for n = 1:length(OutputStr)
    if ( strcmp(OutputStr{n},'i_d') || strcmp(OutputStr{n},'i_q') ...
      || strcmp(OutputStr{n},'id') || strcmp(OutputStr{n},'iq'))
        CountOut = CountOut+1;
        feedout(CountOut) = n;
    end
end

% Connect
GsysObj = obj_Feedback(GmObj,ZbusObj,feedin,feedout);
obj_CheckDim(GsysObj);

% Output auxiliary data
[~,GsysDSS] = GsysObj.ReadDSS(GsysObj);	% Descriptor state space model
vport = feedin;                         % Index of vdq ports in input vector
iport = feedout;                        % Index of idq ports in output vector

end