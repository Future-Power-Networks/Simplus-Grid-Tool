% State space model of a grid-following inverter
%
% Author(s): Yitong Li

function Tss = StateSpaceEquGfl(kp_pll,ki_pll,w_PLL_LPF,Enable_PLL_LPF)

if Enable_PLL_LPF == 0 
% Without PLL LPF
% dx/dt = [dw_pll_i]/dt = [0 ,0]*[w_pll_i] + [ki]*[W]
%         [dtheta  ]      [1, 0] [theta  ]   [kp]
% y = [omega] = [1, 0]*[w_pll_i] + [kp]*[W]
%     [theta]   [0 ,1] [theta  ]   [0 ]
Ai = [0, 0;
      1, 0];
Bi = [ki_pll;
      kp_pll];
Ci = [1, 0;
      0, 1];
Di = [kp_pll;
      0];
Tss = ss(Ai,Bi,Ci,Di);
Tss = Tss*1/ki_pll;

else
% With PLL LPF
% dx/dt = [dw      ]/dt
%         [dw_pll_i]
%         [dtheta  ]
% y = [omega]
%     [theta]
tau = 1/w_PLL_LPF;
Ai = [-1/tau,1/tau,0;
      0,0,0;
      1,0,0];
Bi = [kp_pll/tau;
      ki_pll;
      0];
Ci = [1,0,0;
      0,0,1];
Di = [0;
      0];
Tss = ss(Ai,Bi,Ci,Di);
Tss = Tss*1/(w_PLL_LPF*kp_pll);
end

end