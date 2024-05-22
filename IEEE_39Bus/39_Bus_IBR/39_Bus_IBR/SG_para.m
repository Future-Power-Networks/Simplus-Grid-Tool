%% Simulation Parameters (for the user - optional)
fault_sim = 0;
Bus_fault = [1:39];
t = 200*ones(1,39);
Ts = 5e-5;

%% Line
line=[...
    1	2	0.0035	0.0411	0.6987	0	100	345
    1	39	0.001	0.025	0.75	0	100	345
    2	3	0.0013	0.0151	0.2572	0	100	345
    2	25	0.007	0.0086	0.146	0	100	345
    2	30	0	    0.00181	0	  1.025	100	22
    3	4	0.0013	0.0213	0.2214	0	100	345
    3	18	0.0011	0.0133	0.2138	0	100	345
    4	5	0.0008	0.0128	0.1342	0	100	345
    4	14	0.0008	0.0129	0.1382	0	100	345
    5	8	0.0008	0.0112	0.1476	0	100	345
    6	5	0.0002	0.0026	0.0434	0	100	345
    6	7	0.0006	0.0092	0.113	0	100	345
    6	11	0.0007	0.0082	0.1389	0	100	345
    7	8	0.0004	0.0046	0.078	0	100	345
    8	9	0.0023	0.0363	0.3804	0	100	345
    9	39	0.001	0.025	1.2	    0	100	345
    10	11	0.0004	0.0043	0.0729	0	100	345
    10	13	0.0004	0.0043	0.0729	0	100	345
    10	32	0	    0.002	0	  1.07	100	22
    12	11	0.0016	0.0435	0	  1.006	100	345
    12	13	0.0016	0.0435	0	  1.006	100	345
    13	14	0.0009	0.0101	0.1723	0	100	345
    14	15	0.0018	0.0217	0.366	0	100	345
    15	16	0.0009	0.0094	0.171	0	100	345
    16	17	0.0007	0.0089	0.1342	0	100	345
    16	19	0.0016	0.0195	0.304	0	100	345
    16	21	0.0008	0.0135	0.2548	0	100	345
    16	24	0.0003	0.0059	0.068	0	100	345
    17	18	0.0007	0.0082	0.1319	0	100	345
    17	27	0.0013	0.0173	0.3216	0	100	345
    19	33	0.0007	0.00142	0	  1.07	100	22
    19	20	0.0007	0.0138	0	  1.06	100	345
    20	34	0.0009	0.0018	0	  1.009	100	22
    21	22	0.0008	0.014	0.2565	0	100	345
    22	23	0.0006	0.0096	0.1846	0	100	345
    22	35	0	    0.00143	0	  1.025	100	22
    23	24	0.0022	0.035	0.361	0	100	345
    23	36	0.0005	0.00272	0	    1	100	22
    25	26	0.0032	0.0323	0.513	0	100	345
    25	37	0.0006	0.00232	0	  1.025	100	22
    26	27	0.0014	0.0147	0.2396	0	100	345
    26	28	0.0043	0.0474	0.7802	0	100	345
    26	29	0.0057	0.0625	1.029	0	100	345
    28	29	0.0014	0.0151	0.249	0	100	345
    29	38	0.0008	0.00156	0	  1.025	100	345
    31	6	0	    0.0025	0	    1	100	22];    %pu based on 100MW

Bus = [ ...
    1   345
    2   345
    3   345
    4   345
    5   345
    6   345
    7   345
    8   345
    9   345
    10   345
    11   345
    12   230
    13   345
    14   345
    15   345
    16   345
    17   345
    18   345
    19   345
    20   345
    21   345
    22   345
    23   345
    24   345
    25   345
    26   345
    27   345
    28   345
    29   345
    30    22
    31    22
    32    22
    33    22
    34    22
    35    22
    36    22
    37    22
    38   345
    39   345];

s=1;
s1=10;

zbase=(line(:,8).^2)./line(:,7);
line(:,3)=line(:,3).*zbase;
line(:,4)=line(:,4).*zbase/(120*pi);
line(:,5)=line(:,5)./zbase/(120*pi);
I_base = 100./Bus(:,2)*1000;

%% Synchronous machine parameters
% Nominal power and frequency
num_gen = 10;
Pn_sm = [1000,1000,1000,1000,1000,1000,1000,1000,1000,1000];              % nominal rated active power [MW]

cosPhi_sm = ones(1,num_gen);              % nominal rated reactive power [MVar]
Un_sm = [345,22,22,22,22,22,22,22,345,22];                % nominal rated voltage [kV]
fn_sm = 60*ones(1,10);                % nominal rated frequency [Hz]

Kd_sm = 6*ones(1,10);                    % damping factor [p.u.]
H_sm = [5,3.03,3.58,2.86,2.6,3.48,2.64,2.43,3.45,4.2]*10;       % inertia coefficient [s]

% Machine setpoints
wn_sm = 1;
wb_sm = 60*2*pi;
Pref_sm = [1000 520.81 650 632 508 650 560 540 830 250]./Pn_sm;%Active Power Generation of PV units;                   % power reference [p.u.]
Vref_sm = 1*ones(1,10);                        % excitation voltage refrence [p.u.]

% Computed parameters
Sn_sm = Pn_sm./cosPhi_sm;           % nominal rated power factor
Qn_sm = sqrt(Sn_sm.^2-Pn_sm.^2);    % nominal rated apparent power [MVA]
In_sm = Sn_sm./Un_sm;               % nominal rated current [kA]
Zb_sm = Un_sm.^2./Sn_sm;
Lb_sm = Zb_sm./(2*pi*fn_sm);

% Per unit reactances
xd = 1.305;                         % synchronous reactance in d-axis [p.u.]
xd_t = 0.296;                       % transient reactance in d-axis [p.u.]
xd_st = 0.252;                      % subtransient reactance in d-axis [p.u.]
xq = 0.474;                         % synchronous reactance in q-axis [p.u.]
xq_t = 0.5;                         % transient reactance in q-axis [p.u.]
xq_st = 0.243;                      % subtransient reactance in q-axis [p.u.]
xl = 0.18;                          % armature leakage reactance [p.u.]
% Time constants
Td_t = 1.01;                        % transient short-circuit time constant [s]
Td_st = 0.053;                      % subtransient short-circuit time constant [s]
Tq_t = 0.6;                         % transient short-circuit time constant [s]
Tq_st = 0.1;                        % subtransient short-circuit time constant [s]

% Other parameters
rs = 2.8544e-3;                     % stator resistance [p.u.]
p = 1;                              % number of pole pairs

%rs = 0.02;
xd = xd_t;

%% Initial conditions
dw0_sm = 0;                         % frequency deviation [%]
theta0_sm = -94.3*0;                  % rotor angle [deg]
ia0_sm = 0.75*0;                      % stator currents in phase a [p.u.]
ib0_sm = 0.75*0;                      % stator currents in phase b [p.u.]
ic0_sm = 0.75*0;                      % stator currents in phase c [p.u.]
pha0_sm = -24.95*0;                   % angle of phase a [deg]
phb0_sm = -144.94*0;                  % angle of phase b [deg]
phc0_sm = 95.05*0;                    % angle of phase c [deg]
vf0_sm = 1.29;                      % excitation voltage [p.u.]

%% Ideal transformer
SnT = Sn_sm;                        % nominal rated apparent power [MVA]
PnT = SnT;                          % nominal rated active power [MW]
UnT1 = Un_sm;                       % nominal rated voltage on the low-voltage side [kV]
UnT2 = 230*ones(size(UnT1));        % nominal rated voltage on the high-voltage side [kV]
fnT = [60,60,60,60,60,60];          % nominal rated transformer frequency [Hz]
rmT = [1e4,1e4,1e4];    % magnetization resistance [p.u.]
lmT = [inf,inf,inf];    % magnetization inductance [p.u.]

%% Transmission line
rlp = 0.03;                         % line resistance [Ohm/km]
xlp = 0.3;                          % line reactance [Ohm/km]
llp = xlp/(2*pi*fn_sm(1));          % line inductance [H/km]
clp = 10e-9;                        % line shunt capacitance [F/km]
Ll = [25,10,110,10,25,110,110,10,25];  % line length [km]
rl = rlp * Ll;
ll = llp * Ll;
cl = clp * Ll;

%% Constant impedance load
UnL = [UnT2,UnT2,UnT2,UnT2,UnT2];   % nominal load voltage [kV]
fnL = [fnT,fnT,fnT,fnT,fnT];        % nomnal load frequency [Hz]
% PL =  [400,567,490,1000,100,1570];    % active load power [MW]
%PL =  [400,567,490,1000,100,900];    % active load power [MW]
PL =  [400,367,390,750,100,400];    % active load power [MW]

QL =  0*[0,100,0,100,0,0];            % reactive load power (inductive) [MVar]
QC =  0*[0,200,0,350,0,0];          % reactive load power (capacitive) [MVar]
dPL = 120;                          % increase in load [MW]

%% Governor
R = 1e-3;       % controller droop gain
T1 = 0.05;      % governor time constant [sec]
T2 = 0.5;%0.5;       % turbine derivative time constant [sec]
T3 = 0.5;% 1.5;       % turbine delay time constant [sec]

Tg = [10,10,10];
Kg = [0.95,0.95,0.95];
Fg = [0.1,0.1,0.1];
Rg = [0.05,0.05,0.05];

Dt = 0;         % frictional losses factor
Vmin = 0;       % minimum gate limit
Vmax = inf;       % maximum gate limit
wref = 1;

%% Automatic Voltage Regulator (AVR)
Ke = 200;       % AVR controller gain [150-300]
Te = 0.05;      % exciter time constant [sec]
Tfa = 3;        % filter derivative time constant [sec]
Tfb = 10;       % filter delay time [sec]
Emin = 0.5;       % controller minimum output
Emax = 1;       % controller maximum output
Kas = 2.16e-4;  % AVR scaling gain

%% Power System Stabilizer (PSS)
Ts1 = 0.18;
Ts2 = 2;
Ts3 = 5;
Ts4 = 10;
Ks = 2;
Ts5 = 15;

%% Transformer
Pn_TF = 1e9;


%% ==================== Solver Setup ==================== %%
hCs = getActiveConfigSet(model.name);
hS = hCs.getComponent('Solver');
hS.StartTime = num2str(simu.T_start);
hS.StopTime = num2str(simu.T_stop);
switch simu.Sim_type
    case 'fixed'
        hS.SolverType = 'Fixed-step';
        solverPrefix = 'FixedStep';
        hS.FixedStep = num2str(simu.T_step);
    case 'variable'
        hS.SolverType = 'Variable-step';
        solverPrefix = 'VariableStep';
end
hS.Solver = [solverPrefix, simu.Solver];

%% Simulation
open_system(model.name)

%% fault simulation
% for i = 1:39
%     t = 200*ones(1,39);
%     t(i) = 20;
%     sim('NE39bus2_PQ.slx');
%     save(['If_' num2str(i)],['If_' num2str(i)])
%     i
% end

%% SCC plot
% Isc_sim = zeros(39,1);
% for i = 1:39
%     load(['If_' num2str(i)])
%     I_sim = eval(['If_' num2str(i)]);
%     Isc_sim(i,1) = -max(abs(I_sim.Data))/(I_base(i)*sqrt(2));
% end
% 
% % Isc_11111 = Isc(1023,:)';
% figure(10)
% plot([1:39],abs(Isc_sim))

%%
% if fault_sim == 1
    for i = 33:38
        if any(i == Bus_fault)
            t = 200*ones(1,39);
            t(i) = 10;
            sim('IEEE_39_IBR_Simplified');
            save(['If_' num2str(i)],['If_' num2str(i)])
            save(['Data_' num2str(i)], ['Data_' num2str(i)])
            i
        end
    end
    % Isc_sim = zeros(39,1);
    % for i = 1:39
    %     if any(i == Bus_fault)
    % 
    %         load(['If_' num2str(i)])
    %         I_sim = eval(['If_' num2str(i)]);
    %         Isc_sim(i,1) = -max(abs(I_sim.Data))/(I_base(i))./sqrt(2);
    %     end
    % end
%     SG_binary = '1111111111';
%     SG_binary = '1111110000';
%     SG_index = bin2dec(SG_binary);
%     Isc_11111 = Isc(length(Isc)/2 + SG_index,:)'./10;
%     Nend = 39;
%     plot([1:Nend],abs(Isc_11111(1:Nend)),'LineWidth',4)
%     hold on
%     plot([1:Nend],abs(Isc_sim(1:Nend)),'LineWidth',4)
%     diff = Isc_11111-Isc_sim
%     %     plot([1:Nend],abs(diff(1:Nend)),'LineWidth',4)
%     legend('Analytical results','Simulation results')%,'diff')
%     grid on;
%     box on;
%     set(gca,'fontsize',45);
%     xlabel('Bus','Interpreter','LaTex','FontSize',45);
%     ylabel('SCC [p.u.]','Interpreter','LaTex','FontSize',45);
%     
% end

% f9 = [downsample(w9.time,40)-28, downsample(w9.Data,40).*60-0.02];
% index = find(f9(:,1)==2);
% f9(1:index,2) = 60;
% plot(f9(:,1),f9(:,2))
% xlim([0,70])
% 
% p9 = [downsample(P9.time,40)-28, downsample(P9.Data,40)];
% index = find(p9(:,1)==2);
% % p9(1:index,2) = 0.843150322344298;
% plot(p9(:,1),p9(:,2))
% xlim([0,70])
% 
% v9 = [downsample(V9.time,40)-28, downsample(V9.Data,40)];
% index = find(v9(:,1)==2);
% plot(v9(:,1),v9(:,2))
% xlim([0,70])

% f_pll = [downsample(fpll.time,100)-48, downsample(fpll.Data,100)+50];
% plot(f_pll(:,1),f_pll(:,2))
% xlim([0,10])



%% Save the current profile and re-plot.
% save(['If_4_Case2'],['If_4'])
% 
% 
% 
figure(10)
plot(If_17/1e3,'linewidth',1.5)
set(gca,'xlim',[50 50.2],'fontsize',18)
xlabel('Time (s)','FontSize',18)
ylabel('SCC (kA)','FontSize',18)
title([],'FontSize',18)
box off
I_sim = eval(['If_2']);
Isc_sim = -max(abs(I_sim.Data))/(I_base(4))./sqrt(2)




%% Power flow (pu)
% PG1 = 0.8568;  VG1 = 1.069;
% PG2 = 0.381;  VG2 = 1.086;
% PG3 = 0.5092;  VG3 = 1.098;
% PG4 = 0.4918;  VG4 = 1.141;
% PG5 = 0.3682;  VG5 = 1.153;
% PG6 = 0.5094;  VG6 = 1.133;
% PG7 = 0.4199;  VG7 = 1.138;
% PG8 = 0.3986;  VG8 = 0.937;
% PG9 = 0.6893;  VG9 = 1.075;
% PG10 = 0.1093; VG10 = 1.026;