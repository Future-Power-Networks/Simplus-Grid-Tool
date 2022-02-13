% State space model of a synchronous generator
%
% Author(s): Yitong Li

function [Tss] = StateSpaceEquSg(J,D)

% State space form:
% dx/dt = [domega]/dt = [-D/J,0]*[omega] + [1/J]*[W];
%         [dtheta]      [1   ,0] [theta]   [0  ]
% y = [omega] = [1,0]*[omega] + [0]*[W]
%     [theta]   [0,1] [theta]   [0]
Av = [-D/J, 0;
      1,    0];
Bv = [1/J;
      0];
Cv = [1,0;
      0,1];
Dv = [0;
      0];
Tss = ss(Av,Bv,Cv,Dv);
Tss = Tss*J;    	% T_V_ss actually represents this system: 1/(D/J+s) without Hinv

end