%%
fprintf('Calculating the linearized network matrix: K and FreqShift...\n')
% Get K matrix
for m = 1:N_Bus
    for n = 1:N_Bus
        if n ~= m
            ang_K(m,n) = mu(m) - angle(S(m,n)) - angle(Gbus(m,n));
            K(m,n) = abs(S(m,n))*abs(Gbus(m,n))*sin(ang_K(m,n));
        end
    end
end
for m = 1:N_Bus
    K_temp = 0;
    for n = 1:N_Bus
        if n~=m
            K_temp = K_temp - K(m,n);
        end
    end
    K(m,m) = K_temp;
end
K = double(K);
K = -K;                             % For negative feedback

% Get FreqShift matrix
for m = 1:N_Bus
    for n = 1:N_Bus
        ang_FreqShift(m,n) = mu(m) - angle(S(m,n)) - angle(Gbus_prime(m,n));
        FreqShift(m,n) = abs(S(m,n))*abs(Gbus_prime(m,n))*sin(ang_FreqShift(m,n));
    end
end
FreqShift = double(FreqShift);
FreqShift = -FreqShift;             % For negative feedback

% Notes:
% K can be interpreted as the synchronizing torque coefficient, as dS is
% proportional to K*dtheta. FreqShift can be interpreted as the damping
% torque coefficient, as dS is proportional to FreqShift*dw. K and
% FreqShift are all linearized results.

%% K analysis
fprintf('Analyze K...\n')
KH = Hinv*K;
[KH11,KH12,KH21,KH22] = SimplusGT.PartitionMatrix(KH,n_Ibus_1st-1,n_Ibus_1st-1);
[phi,xi,psi] = eig(KH);                 % phi is the right eigenvector matrix, psi is the left eigenvector matrix, xi is the eigenvalue matrix. Noting that phi^(-1) can also be regarded as a left eigenvector matrix.
phi_inv = phi^(-1);
xi_diag = diag(xi);
xi_diag = vpa(xi_diag,5)
[xi_min,xi_min_index] = min( real(xi_diag) );
fprintf(['KH stability: ']);
if xi_min < -1e-5
    fprintf(['unstable xi.\n']);
    xi_min;
    xi_min_index;
else
    fprintf(['stable xi.\n']);
  	xi_min;
    xi_min_index;
end

%% Apparatus Matrix: T
fprintf('Calculate the apparatus state space model...\n')

% Notes
% The representation of Hinv is very important, especially when considering
% the whole system KH. When seperately considerring KH_V and KH_I, H_V >>
% H_I, or equivalently, Hinv_V << Hinv_I, should be valid so that the
% voltage source is much slower than the current source.

% Choose the reference node
Select_Ibus_Ref = 1;
if Select_Ibus_Ref == 1
    n_i_ref = n_Ibus_1st;  	% Select the first current node as the reference
elseif Select_Ibus_Ref == 2
    n_i_ref = N_Bus;        % Select the final current node as the reference
else
    error(['Error;']);
end
n_v_ref = 1;                % Select the first voltage node as the reference

% ================================
% Voltage node
% ================================
if Exist_Vbus == 1

% State space form:
% dx/dt = [domega]/dt = [-D/J,0]*[omega] + [1/J]*[W];
%         [dtheta]      [1   ,0] [theta]   [0  ]
% y = [omega] = [1,0]*[omega] + [0]*[W]
%     [theta]   [0,1] [theta]   [0]
Av = [-D{n_v_ref}/J{n_v_ref}, 0;
      1,                      0];
Bv = [1/J{n_v_ref};
      0];
Cv = [1,0;
      0,1];
Dv = [0;
      0];
T_V_ss = ss(Av,Bv,Cv,Dv);
T_V_ss = T_V_ss*J{n_v_ref};         % T_V_ss actually represents this system: 1/(D/J+s) without Hinv

end

% ================================
% Current node
% ================================
if Exist_Ibus == 1

% State space form
if Enable_PLL_LPF == 0                                                        	% ???
% Without PLL LPF
% dx/dt = [dw_pll_i]/dt = [0 ,0]*[w_pll_i] + [ki]*[W]
%         [dtheta  ]      [1, 0] [theta  ]   [kp]
% y = [omega] = [1, 0]*[w_pll_i] + [kp]*[W]
%     [theta]   [0 ,1] [theta  ]   [0 ]
Ai = [0, 0;
      1, 0];
Bi = [ki_pll{n_i_ref};
      kp_pll{n_i_ref}];
Ci = [1, 0;
      0, 1];
Di = [kp_pll{n_i_ref};
      0];
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
Bi = [kp_pll{n_i_ref}/tau;
      ki_pll{n_i_ref};
      0];
Ci = [1,0,0;
      0,0,1];
Di = [0;
      0];
end
T_I_ss = ss(Ai,Bi,Ci,Di);
T_I_ss = T_I_ss/Hinv(n_i_ref,n_i_ref);

end

%% Calculate the state space representation
fprintf('Calculate the whole system state space model...\n')
% Get whole system Tss
Tss = [[],[],[],[]];
for i = 1:N_Bus
    if ApparatusSourceType(i) == 1
        Tss = append(Tss,T_V_ss);
    elseif ApparatusSourceType(i) == 2
        Tss = append(Tss,T_I_ss);
    else
        error(['Error']);
    end
end

% Calculate the whole-system closed loop state space model
feedin = [1:N_Bus];
feedout_L1 = [1:N_Bus]*2;           % theta port
feedout_L2 = [1:N_Bus]*2-1;         % omega port
T1cl = feedback(Tss*Hinv,K,feedin,feedout_L1);
T12cl = feedback(T1cl,FreqShift,feedin,feedout_L2);
if ~isempty(T1cl.E) || ~isempty(T12cl.E)
    error(['Error: T1cl or T2cl is a dss system.']);
end
% T1cl = minreal(T1cl);
% T12cl = minreal(T12cl);

% Calculate the whole system pole
if 1
    [~,pole_T1cl] = eig(T1cl.A);
    pole_T1cl = diag(pole_T1cl)/2/pi;
    [~,pole_T12cl] = eig(T12cl.A);
    pole_T12cl = diag(pole_T12cl)/2/pi;
else
    pole_T1cl = pole(T1cl)/2/pi;
    pole_T12cl = pole(T12cl)/2/pi;
end

%% Plot
fprintf('Plot...\n')

% Plot: poles of state space system
if Enable_Plot_Eigenvalue
Fig_N = Fig_N+1;
figure(Fig_N)
scatter(real(pole_T1cl),imag(pole_T1cl),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_T12cl),imag(pole_T12cl),'x','LineWidth',1.5); hold on; grid on;

legend('Loop1','Loop12')
end

if Enable_ComparisonToolbox
Fig_N = Fig_N+1;
figure(Fig_N)
scatter(real(pole_T1cl),imag(pole_T1cl),'x','LineWidth',1.5); hold on; grid on;
scatter(real(pole_T12cl),imag(pole_T12cl),'x','LineWidth',1.5); hold on; grid on;
pole_toolbox = load('pole_sys').pole_sys;
index = find(abs(imag(pole_toolbox))<35);
pole_toolbox = pole_toolbox(index);
index = find(real(pole_toolbox)>-1e3);
pole_toolbox = pole_toolbox(index);
scatter(real(pole_toolbox),imag(pole_toolbox),'x','LineWidth',1.5); hold on; grid on;
legend('Wihtout Freq Shift','With Freq Shift','Toolbox')
end

%% Check stability
fprintf('Check the stability by poles:\n')
UnstablePoleIndex     = find(real(pole_T12cl)>1e-9);
UnstablePoleIndexRisk = find(real(pole_T12cl)>0);
if isempty(UnstablePoleIndex)
    fprintf('Stable!\n')
else
    fprintf('Unstable!\n')
end
% UnstablePole = pole_T12cl(UnstablePoleIndex)
RiskPole = pole_T12cl(UnstablePoleIndexRisk)