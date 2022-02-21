%%
fprintf('Calculate the linearized network matrix: K and FreqShift...\n')
% Get K matrix
for m = 1:N_Node
    for n = 1:N_Node
        if n ~= m
            ang_K(m,n) = mu(m) - angle(S(m,n)) - angle(Gbus(m,n));
            K(m,n) = abs(S(m,n))*abs(Gbus(m,n))*sin(ang_K(m,n));
        end
    end
end
for m = 1:N_Node
    K_temp = 0;
    for n = 1:N_Node
        if n~=m
            K_temp = K_temp - K(m,n);
        end
    end
    K(m,m) = K_temp;
end
K = double(K);
K = -K;                             % For negative feedback

% Get FreqShift matrix
for m = 1:N_Node
    for n = 1:N_Node
        ang_FreqShift(m,n) = mu(m) - angle(S(m,n)) - angle(GbusPrime(m,n));
        FreqShift(m,n) = abs(S(m,n))*abs(GbusPrime(m,n))*sin(ang_FreqShift(m,n));
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
% [KH11,KH12,KH21,KH22] = SimplusGT.PartitionMatrix(KH,NumIbus1st-1,NumIbus1st-1);
[phi,xi,~] = eig(KH);        	% phi is the right eigenvector matrix, psi is the left eigenvector matrix, xi is the eigenvalue matrix. Noting that phi^(-1) can also be regarded as a left eigenvector matrix.
% phi_inv = phi^(-1);
xi_diag = diag(xi);
[XiMin,XiMinIndex] = min( real(xi_diag) );
fprintf(['KH stability: ']);
if XiMin < -1e-5
    fprintf(['unstable xi.\n']);
else
    fprintf(['stable xi.\n']);
end
XiMin;
XiMinIndex;

%% Apparatus Matrix: T
fprintf('Calculate the apparatus state space model...\n')

if ExistVbus == 1
% Select the first voltage node as the reference
TssV = SimplusGT.Synchron.StateSpaceEquSg(J{1},D{1});
end

if ExistIbus == 1
% Select the first current node as the reference     
TssI = SimplusGT.Synchron.StateSpaceEquGfl(kp_pll{NumIbus1st},ki_pll{NumIbus1st},w_PLL_LPF,Enable_PLL_LPF);
end

% Notes:
% This needs J is proportional to D, and kp_pll is proportional to ki_pll.

%% Calculate the state space representation
fprintf('Calculate the whole system state space model...\n')
% Get whole system Tss
Tss = [[],[],[],[]];
for i = 1:N_Node
    if ApparatusSourceType(i) == 1
        Tss = append(Tss,TssV);
    elseif ApparatusSourceType(i) == 2
        Tss = append(Tss,TssI);
    else
        error(['Error']);
    end
end

% Calculate the whole-system closed loop state space model
feedin = [1:N_Node];
feedoutL1 = [1:N_Node]*2;           % theta port
feedoutL2 = [1:N_Node]*2-1;         % omega port
T1cl = feedback(Tss*Hinv,K,feedin,feedoutL1);
T12cl = feedback(T1cl,FreqShift,feedin,feedoutL2);
if ~isempty(T1cl.E) || ~isempty(T12cl.E)
    error(['Error: T1cl or T2cl is a dss system.']);
end
% T1cl = minreal(T1cl);
% T12cl = minreal(T12cl);

fprintf('Calculate the poles:\n')
[~,PoleT1cl] = eig(T1cl.A);
PoleT1cl = diag(PoleT1cl)/2/pi;
[~,PoleT12cl] = eig(T12cl.A);
PoleT12cl = diag(PoleT12cl)/2/pi;
% pole_T1cl = pole(T1cl)/2/pi;
% pole_T12cl = pole(T12cl)/2/pi;

%% Check stability
fprintf('Check stability: ')
IndexUnstablePole   = find(real(PoleT12cl)>1e-9);
IndexRiskPole       = find(real(PoleT12cl)>0);
if isempty(IndexUnstablePole)
    fprintf('stable poles!\n')
else
    fprintf('unstable poles!\n')
end
% UnstablePole = pole_T12cl(IndexUnstablePole)
RiskPole = PoleT12cl(IndexRiskPole)

%% Plot
fprintf('Plot...\n')

% Plot: poles of state space system
FigN = FigN+1;
figure(FigN)
scatter(real(PoleT1cl),imag(PoleT1cl),'x','LineWidth',1.5); hold on; grid on;
scatter(real(PoleT12cl),imag(PoleT12cl),'x','LineWidth',1.5); hold on; grid on;
legend('Loop1','Loop12')
SimplusGT.mtit('Pole map');

%if Enable_ComparisonToolbox
%Fig_N = Fig_N+1;
%figure(Fig_N)
%scatter(real(pole_T1cl),imag(pole_T1cl),'x','LineWidth',1.5); hold on; grid on;
%scatter(real(pole_T12cl),imag(pole_T12cl),'x','LineWidth',1.5); hold on; grid on;
%pole_toolbox = load('pole_sys').pole_sys;
%index = find(abs(imag(pole_toolbox))<35);
%pole_toolbox = pole_toolbox(index);
%index = find(real(pole_toolbox)>-1e3);
%pole_toolbox = pole_toolbox(index);
%scatter(real(pole_toolbox),imag(pole_toolbox),'x','LineWidth',1.5); hold on; grid on;
%legend('Wihtout Freq Shift','With Freq Shift','Toolbox')
%end