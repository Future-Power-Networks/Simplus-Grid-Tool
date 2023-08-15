% Use the toolbox results to acquire whole-system impedance model
% Output: the state-space model of whole-system impedance.
% Method: Zsys = feedback(Zm,Yb)
%    1) Switch i_out and v_in of Gm, lead to Zm.
%    2) Feedback Zm and Ybus.
%Author: Yue Zhu


function Zsys_SS = WholeSysZ_cal(GmObj,YbusObj,Port_i,Port_v)

[~,Gm_dss] = GmObj.GetDSS(GmObj);
[~,YbusDSS] = YbusObj.GetDSS(YbusObj);
Gm_dss_trim = Gm_dss(Port_i,Port_v);
Zm_dss = inv(Gm_dss_trim);
Zsys_dss = feedback(Zm_dss,YbusDSS);
Zsys_SS = SimplusGT.dss2ss(Zsys_dss);

end