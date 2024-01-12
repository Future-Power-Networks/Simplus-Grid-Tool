
axis_sel=2; % input selection: 1 for d, 2 for q

if axis_sel==1
  vdq2_d=out.vdq2;
  vdq3_d=out.vdq3;
  vdq_all_d=out.vdq_all;
  idq_in_d =out.i_in_dq;
  save("Dome/ERA_DQ/vdq2_d","vdq2_d");
  save("Dome/ERA_DQ/vdq3_d","vdq3_d");
  save("Dome/ERA_DQ/vdq_all_d","vdq_all_d");
  save("Dome/ERA_DQ/idq_in_d","idq_in_d");
elseif axis_sel==2
    vdq2_q=out.vdq2;
    vdq3_q=out.vdq3;
    vdq_all_q=out.vdq_all;
    idq_in_q =out.i_in_dq;
    save("Dome/ERA_DQ/vdq2_q","vdq2_q");
    save("Dome/ERA_DQ/vdq3_q","vdq3_q");
    save("Dome/ERA_DQ/vdq_all_q","vdq_all_q");
    save("Dome/ERA_DQ/idq_in_q","idq_in_q");
end

save("Dome/ERA_DQ/ListPowerFlowNew","ListPowerFlowNew")

Zsys_SS = SimplusGT.WholeSysZ_cal(GmObj,YbusObj,Port_i,Port_v); % theoritical results
save("Dome\ERA_DQ\Zsys_SS","Zsys_SS");
YN=YbusDSS;
save("Dome\ERA_DQ\YN","YN");