% This function connects the admittance model of apparatuses and the impedance
% model of network lines.

% Author(s): Yitong Li

function [GsysObj,GsysDSS,Port_v_feedin,Port_i_feedout,BusPort_v,BusPort_i] = ConnectGmZbus(GmObj,ZbusObj,N_Bus)

% Get the strings
[~,GmInStr,GmOutStr] = GmObj.GetString(GmObj);

% Initialize
Port_v_feedin = [];
Port_i_feedout = [];

% Port_w  = zeros(1,N_Bus);
% Port_T_m = zeros(1,N_Bus);
% Port_ang_r = zeros(1,N_Bus);
% Port_P_dc = zeros(1,N_Bus);
% Port_v_dc = zeros(1,N_Bus);

% Find ports
for i = 1:N_Bus
    [~,in1] = SimplusGT.CellFind(GmInStr,['v_d',num2str(i)]);
    [~,in2] = SimplusGT.CellFind(GmInStr, ['v',num2str(i)]);
    
    [~,out1] = SimplusGT.CellFind(GmOutStr,['i_d',num2str(i)]);
    [~,out2] = SimplusGT.CellFind(GmOutStr,['i',num2str(i)]);
    if ~isempty(in1)
        Port_v_feedin = [Port_v_feedin,in1,in1+1];
        Port_i_feedout = [Port_i_feedout,out1,out1+1];
        BusPort_v{i} = [in1,in1+1];
        BusPort_i{i} = [out1,out1+1];
    elseif ~isempty(in2)
        Port_v_feedin = [Port_v_feedin,in2];
        Port_i_feedout = [Port_i_feedout,out2];
        BusPort_v{i} = in2;
        BusPort_i{i} = out2;
    else
        error(['Error']);
    end
    
%  	% Other input ports
%     [~,p1] = SimplusGT.CellFind(GmInStr,['T_m',num2str(i)]);
%     if ~isempty(p1) 
%         Port_T_m(i) = p1;
%     end
% 	[~,p2] = SimplusGT.CellFind(GmInStr,['ang_r',num2str(i)]);
%     if ~isempty(p2) 
%         Port_ang_r(i) = p2;
%     end
%    	[~,p3] = SimplusGT.CellFind(GmInStr,['P_dc',num2str(i)]);
%     if ~isempty(p3) 
%         Port_P_dc(i) = p3;
%     end
%         
%     % Output output ports
%   	[~,p4] = SimplusGT.CellFind(GmOutStr,['w',num2str(i)]);
%     if ~isempty(p4) 
%         Port_w(i) = p4;
%     end
% 	[~,p5] = SimplusGT.CellFind(GmOutStr,['v_dc',num2str(i)]);
%     if ~isempty(p5) 
%         Port_v_dc(i) = p5;
%     end
end

% Connect
GsysObj = SimplusGT.ObjFeedback(GmObj,ZbusObj,Port_v_feedin,Port_i_feedout);
SimplusGT.ObjCheckDim(GsysObj);

% Output auxiliary data
[~,GsysDSS] = GsysObj.GetDSS(GsysObj);	% Descriptor state space model

end