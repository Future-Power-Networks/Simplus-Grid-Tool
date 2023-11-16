% 
% vdq2_d=out.vdq2;
% save("Dome/ERA_DQ/vdq2_d","vdq2_d");

vdq2_q=out.vdq2;
save("Dome/ERA_DQ/vdq2_q","vdq2_q");

save("Dome/ERA_DQ/ListPowerFlowNew","ListPowerFlowNew")

Zsys_SS = SimplusGT.WholeSysZ_cal(GmObj,YbusObj,Port_i,Port_v); % theoritical results
save("Dome\ERA_DQ\Zsys_SS","Zsys_SS");