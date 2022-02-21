% This function handles different nodes in the power networks
%
% Author(s): Yitong Li

function [V,I,J,D,kp_pll,ki_pll,Ybus,Ybus_,N_Node] = HandleNode(N_Bus,V,I,ApparatusParaNew,Ybus,Ybus_,NumVbus1st,NumIbus1st,NumFbus1st,ExistVbus,ExistIbus,ExistFbus,Wbase,dW,Enable_VoltageNode_InnerLoop,Enable_CurrentNode_InnerLoop)

% Symbolic Laplace operator s
s = sym('s');

% =============================
% Handle voltage node
% =============================
J = [];
D = [];
if ExistVbus == 0
    fprintf('Warning: The system has no voltage node.\n')
else
fprintf('Handle voltage node...\n')

for i = 1:(NumIbus1st-1)
J{i} = ApparatusParaNew{i}.J*2/Wbase^2;
D{i} = ApparatusParaNew{i}.D/Wbase^2;
J{i} = J{i}*Wbase;
D{i} = D{i}*Wbase;
% Notes: 
% Adding '*Wbase' is because the power swing equation rather than the
% torque swing equation is used, and P=T*w0 if w is not in per unit system.

Lsg = ApparatusParaNew{i}.wL/Wbase;
Rsg = ApparatusParaNew{i}.R;
Zsg = s*Lsg + Rsg;
Ysg{i} = 1/Zsg;

% Convert Ysg to double
Ysg_{i} = double(subs(Ysg{i},'s',1i*(Wbase+dW)));
Ysg{i} = double(subs(Ysg{i},'s',1i*Wbase));
end    

% Doing D-Y conversion
if Enable_VoltageNode_InnerLoop

% Prepare star-delta conversion by adding new buses
Ybus = SimplusGT.Synchron.PrepareConvertDY(Ybus,NumIbus1st,N_Bus,Ysg);
Ybus_ = SimplusGT.Synchron.PrepareConvertDY(Ybus_,NumIbus1st,N_Bus,Ysg_);
    
% Doing the star-delta conversion.
% Notes: Assume old voltage bus as zero current bus, and then switch the
% current and voltage for these buses so that current becomes input, and
% finally remove corresponding blocks because the input current is zero.
Ybus = SimplusGT.Synchron.HybridMatrixYZ(Ybus,N_Bus+1);
Ybus_ = SimplusGT.Synchron.HybridMatrixYZ(Ybus_,N_Bus+1);

% Eliminate the old voltage bus, i.e., zero current bus
Ybus = Ybus(1:N_Bus,1:N_Bus);
Ybus_ = Ybus_(1:N_Bus,1:N_Bus);

% Update V and I
% Notes: The steady-state voltage at voltage buses are changed if we split
% the inductor outside the apparatus. The "for loop" voltage calculation is
% equivalent to the matrix form "V = inv(Ybus)*I". But this matrix form
% would lead to wrong results in some cases when Ybus is almost
% non-invertible.
I = I(1:N_Bus,:);
for i = 1:(NumIbus1st-1)
    V(i) = V(i) + I(i)/Ysg{i};
end

else
    fprintf('Warning: The voltage-node inner loop has been disabled.\n')
end

end

% =============================
% Handle current node
% =============================
if ExistIbus == 0
    fprintf('Warning: The system has no current node.\n')
else

fprintf('Handle current node...\n')

% Get the inner loop parameters
kp_pll = [];
ki_pll = [];
for i = NumIbus1st:(NumFbus1st-1)

% Notes: 
% All inverters have same current controllers
Rf   = ApparatusParaNew{i}.R;
Lf   = ApparatusParaNew{i}.wLf/Wbase;
kp_i = ApparatusParaNew{i}.f_i_dq*2*pi*Lf;  
ki_i = (ApparatusParaNew{i}.f_i_dq*2*pi)^2*Lf/4;

% Notes:
% The bandwidth of vq-PLL and Q-PLL would be different, because Q is
% proportional to id*vq. If we want to ensure the vq-PLL bandwidth of
% different inverters to be same, there Q-PLL bandwidth would probably be
% different because of different id output or active power output. We
% ensure the vq-PLL bandwidth here.
kp_pll{i} = ApparatusParaNew{i}.f_pll*2*pi; 
ki_pll{i} = (ApparatusParaNew{i}.f_pll*2*pi)^2/4;
PI_pll{i} = kp_pll{i} + ki_pll{i}/s;

end

% alpha/beta
Z_PIi = kp_i + ki_i/(s-1i*Wbase);        
Z_Lf = s*Lf+Rf;
Yinv = (s-1i*Wbase)/((kp_i + s*Lf+Rf)*(s-1i*Wbase) + ki_i);

Yinv_ = double(subs(Yinv,'s',1i*(Wbase+dW)));
Yinv = double(subs(Yinv,'s',1i*Wbase));

% Notes:
% When current controller is very fast, Y_inv -> 0 and can be ignored.

% Add Y_inv to nodal admittance matrix
if Enable_CurrentNode_InnerLoop 
    
for i = NumIbus1st:(NumFbus1st - 1)
    % Self branch
    Ybus(i,i) = Ybus(i,i) + Yinv;
    Ybus_(i,i) = Ybus_(i,i) + Yinv_;
end

% Update I
I = Ybus*V;

else
    fprintf('Warning: The current-node inner loop has been disabled.\n')
end

end


% =============================
% Handle floating node
% =============================
% Notes:
% The floating bus (i.e., no device bus) is assumed as zero-current bus,
% and eliminated here after converting the Y matrix to Y-Z hybrid matrix.
N_Node = NumFbus1st-1;              % N_Node contains Vbus and Ibus only
if ExistFbus == 0
    fprintf('Warning: The system has no floating node.\n')
else
    
fprintf('Eliminate floating node...\n')

Ybus = SimplusGT.Synchron.HybridMatrixYZ(Ybus,NumFbus1st);
Ybus_ = SimplusGT.Synchron.HybridMatrixYZ(Ybus_,NumFbus1st);

Ybus = Ybus(1:N_Node,1:N_Node);
Ybus_ = Ybus_(1:N_Node,1:N_Node);

V = V(1:N_Node,:);
I = I(1:N_Node,:);

end

end