% =============================
% Handle voltage node
% =============================
if Exist_Vbus == 0
    fprintf('Warning: The system has no voltage node.\n')
else
fprintf('Handle voltage node...\n')

for i = 1:(n_Ibus_1st-1)
J{i} = ApparatusParaNew{i}.J*2/W0^2;
D{i} = ApparatusParaNew{i}.D/W0^2;
J{i} = J{i}*Wbase;
D{i} = D{i}*Wbase;
% Notes: 
% Adding '*Wbase' is because the power swing equation rather than the
% torque swing equation is used, and P=T*w0 if w is not in per unit system.

Lsg = ApparatusParaNew{i}.wL/W0;
Rsg = ApparatusParaNew{i}.R;
Zsg = s*Lsg + Rsg;
Y_sg{i} = 1/Zsg;

% Convert Ysg to double
Y_sg_{i} = double(subs(Y_sg{i},'s',1i*(W0+dW)));
Y_sg{i} = double(subs(Y_sg{i},'s',1i*W0));
end    

% Doing D-Y conversion
if Enable_VoltageNode_InnerLoop

% Prepare star-delta conversion by adding new buses
Ybus = SimplusGT.Communication.PrepareConvertDY(Ybus,n_Ibus_1st,N_Bus,Y_sg);
Ybus_ = SimplusGT.Communication.PrepareConvertDY(Ybus_,n_Ibus_1st,N_Bus,Y_sg_);
    
% Doing the star-delta conversion.
% Notes: Assume old voltage bus as zero current bus, and then switch the
% current and voltage for these buses so that current becomes input, and
% finally remove corresponding blocks because the input current is zero.
Ybus = SimplusGT.Communication.HybridMatrixYZ(Ybus,N_Bus+1);
Ybus_ = SimplusGT.Communication.HybridMatrixYZ(Ybus_,N_Bus+1);

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
for i = 1:(n_Ibus_1st-1)
    V(i) = V(i) + I(i)/Y_sg{i};
end

else
    fprintf('Warning: The voltage-node inner loop has been disabled.\n')
end

end

% =============================
% Handle current node
% =============================
if Exist_Ibus == 0
    fprintf('Warning: The system has no current node.\n')
else

fprintf('Handle current node...\n')

% Get the inner loop parameters

for i = n_Ibus_1st:(n_Fbus_1st-1)

% Notes: 
% All inverters have same current controllers
Rf   = ApparatusParaNew{i}.R;
Lf   = ApparatusParaNew{i}.wLf/W0;
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
Z_PIi = kp_i + ki_i/(s-1i*W0);        
Z_Lf = s*Lf+Rf;
Y_inv = (s-1i*W0)/((kp_i + s*Lf+Rf)*(s-1i*W0) + ki_i);

Y_inv_ = double(subs(Y_inv,'s',1i*(W0+dW)));
Y_inv = double(subs(Y_inv,'s',1i*W0));

% Notes:
% When current controller is very fast, Y_inv -> 0 and can be ignored.

% Add Y_inv to nodal admittance matrix
if Enable_CurrentNode_InnerLoop                                                 % ??? 
    
for i = n_Ibus_1st:(n_Fbus_1st - 1)
    % Self branch
    Ybus(i,i) = Ybus(i,i) + Y_inv;
    Ybus_(i,i) = Ybus_(i,i) + Y_inv_;
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
if Exist_Fbus == 0
    fprintf('Warning: The system has no floating node.\n')
    YbusVIF = Ybus;
    YbusVI = Ybus;
else
    
fprintf('Eliminate floating node...\n')

YbusVIF = Ybus;

Ybus = SimplusGT.Communication.HybridMatrixYZ(Ybus,n_Fbus_1st);
Ybus_ = SimplusGT.Communication.HybridMatrixYZ(Ybus_,n_Fbus_1st);

N_Bus = n_Fbus_1st-1;
Ybus = Ybus(1:N_Bus,1:N_Bus);
Ybus_ = Ybus_(1:N_Bus,1:N_Bus);
YbusVI = Ybus;
V = V(1:N_Bus,:);
I = I(1:N_Bus,:);

end