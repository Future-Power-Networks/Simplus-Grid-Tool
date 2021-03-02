% This function connects the admittance model of devices and the impedance
% model of network lines.

% Author(s): Yitong Li

function [GsysObj,GsysDSS,Port_v,Port_i,Port_w,Port_T_m,Port_ang_r,Port_P_dc,Port_v_dc] = ConnectGmZbus(GmObj,ZbusObj,N_Bus)

% Get the strings
[~,GmInStr,GmOutStr] = GmObj.GetString(GmObj);

% Initialize
Port_v = [];
Port_i = [];

Port_w  = zeros(1,N_Bus);
Port_T_m = zeros(1,N_Bus);
Port_ang_r = zeros(1,N_Bus);
Port_P_dc = zeros(1,N_Bus);
Port_v_dc = zeros(1,N_Bus);

% Find ports
for i = 1:N_Bus
    [~,i1] = SimplexPS.CellFind(GmInStr,['v_d',num2str(i)]);
    [~,i2] = SimplexPS.CellFind(GmInStr, ['v',num2str(i)]);
    
    [~,o1] = SimplexPS.CellFind(GmOutStr,['i_d',num2str(i)]);
    [~,o2] = SimplexPS.CellFind(GmOutStr,['i',num2str(i)]);
    if ~isempty(i1)
        Port_v = [Port_v,i1,i1+1];
        Port_i = [Port_i,o1,o1+1];
    elseif ~isempty(i2)
        Port_v = [Port_v,i2];
        Port_i = [Port_i,o2];
    else
        error(['Error']);
    end
    
 	% Other input ports
    [~,p1] = SimplexPS.CellFind(GmInStr,['T_m',num2str(i)]);
    if ~isempty(p1) 
        Port_T_m(i) = p1;
    end
	[~,p2] = SimplexPS.CellFind(GmInStr,['ang_r',num2str(i)]);
    if ~isempty(p2) 
        Port_ang_r(i) = p2;
    end
   	[~,p3] = SimplexPS.CellFind(GmInStr,['P_dc',num2str(i)]);
    if ~isempty(p3) 
        Port_P_dc(i) = p3;
    end
        
    % Output output ports
  	[~,p4] = SimplexPS.CellFind(GmOutStr,['w',num2str(i)]);
    if ~isempty(p4) 
        Port_w(i) = p4;
    end
	[~,p5] = SimplexPS.CellFind(GmOutStr,['v_dc',num2str(i)]);
    if ~isempty(p5) 
        Port_v_dc(i) = p5;
    end
end

% Connect
GsysObj = SimplexPS.ObjFeedback(GmObj,ZbusObj,Port_v,Port_i);
SimplexPS.ObjCheckDim(GsysObj);

% Output auxiliary data
[~,GsysDSS] = GsysObj.GetDSS(GsysObj);	% Descriptor state space model

end