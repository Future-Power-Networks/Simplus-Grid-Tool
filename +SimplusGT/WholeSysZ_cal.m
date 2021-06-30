% Use the toolbox results to acquire whole-system impedance mode
% ###ONLY support AC system at present.###
% Method: Zsys = feedback(Zm,Yb)
%    1) get impedance model for each apparatus Zm
%    2) Link apparatus impedance with network admittance Yb
%Author: Yue Zhu


function [ZsysObj,ZsysDSS] = WholeSysZ_cal(GmObj_Cell,YbusObj,N_Apparatus)

% 1) Switch the input v and output i of all apparatus to change into impedance model Zm
for i = 1: N_Apparatus 
    ZmObj_Cell{i} = SimplusGT.ObjSwitchInOut(GmObj_Cell{i},2);
end
ZmObj = SimplusGT.Toolbox.ApparatusModelLink(ZmObj_Cell);

% 2) acquire current-input port and voltage-ouput port number of Zmobj:
[~,ZmInStr,ZmOutStr] = ZmObj.GetString(ZmObj);%get name string
Port_i_feedin = []; %declare
Port_v_feedout = [];
for i = 1:N_Apparatus
    [~,in1] = SimplusGT.CellFind(ZmInStr,['i_d',num2str(i)]);
    [~,in2] = SimplusGT.CellFind(ZmInStr, ['i',num2str(i)]);
    
    [~,out1] = SimplusGT.CellFind(ZmOutStr,['v_d',num2str(i)]);
    [~,out2] = SimplusGT.CellFind(ZmOutStr,['v',num2str(i)]);
    if ~isempty(in1) %ac apparatus
        Port_i_feedin = [Port_i_feedin,in1,in1+1];
        Port_v_feedout = [Port_v_feedout,out1,out1+1];
    elseif ~isempty(in2) % dc: not activated at this moment.
        Port_i_feedin = [Port_i_feedin,in2];
        Port_v_feedout = [Port_v_feedout,out2];
    else
        error(['Error']);
    end
end

% 3) Feedback Zm and Ybus
ZsysObj = SimplusGT.ObjFeedback(ZmObj,YbusObj,Port_i_feedin,Port_v_feedout);
[~,ZsysDSS] = ZsysObj.GetDSS(ZsysObj);

end