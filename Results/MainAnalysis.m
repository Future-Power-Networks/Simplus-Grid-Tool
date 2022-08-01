% Author(s): Yitong Li

%% 
% Clear
clear all
clc
close all

%%
% Load GFM data
Yss_Gfm_P0_Q0 = load('Yss_Gfm_P0_Q0');

Yss_Gfm_Pp05_Q0 = load('Yss_Gfm_Pp05_Q0');
Yss_Gfm_Pp1_Q0 = load('Yss_Gfm_Pp1_Q0');
Yss_Gfm_Pn05_Q0 = load('Yss_Gfm_Pn05_Q0');
Yss_Gfm_Pn1_Q0 = load('Yss_Gfm_Pn1_Q0');

Yss_Gfm_P0_Qp05 = load('Yss_Gfm_P0_Qp05');
Yss_Gfm_P0_Qp1 = load('Yss_Gfm_P0_Qp1');
Yss_Gfm_P0_Qn05 = load('Yss_Gfm_P0_Qn05');
Yss_Gfm_P0_Qn1 = load('Yss_Gfm_P0_Qn1');

Yss_Gfm_Pp1_Qp1 = load('Yss_Gfm_Pp1_Qp1');

%%
% Convert state space model to symbolic transfer function
Ysym_Gfm_P0_Q0 = SimplusGT.ss2sym(Yss_Gfm_P0_Q0.Yss{2});

Ysym_Gfm_Pp05_Q0 = SimplusGT.ss2sym(Yss_Gfm_Pp05_Q0.Yss{2});
Ysym_Gfm_Pp1_Q0 = SimplusGT.ss2sym(Yss_Gfm_Pp1_Q0.Yss{2});
Ysym_Gfm_Pn05_Q0 = SimplusGT.ss2sym(Yss_Gfm_Pn05_Q0.Yss{2});
Ysym_Gfm_Pn1_Q0 = SimplusGT.ss2sym(Yss_Gfm_Pn1_Q0.Yss{2});

Ysym_Gfm_P0_Qp05 = SimplusGT.ss2sym(Yss_Gfm_P0_Qp05.Yss{2});
Ysym_Gfm_P0_Qp1 = SimplusGT.ss2sym(Yss_Gfm_P0_Qp1.Yss{2});
Ysym_Gfm_P0_Qn05 = SimplusGT.ss2sym(Yss_Gfm_P0_Qn05.Yss{2});
Ysym_Gfm_P0_Qn1 = SimplusGT.ss2sym(Yss_Gfm_P0_Qn1.Yss{2});

Ysym_Gfm_Pp1_Qp1 = SimplusGT.ss2sym(Yss_Gfm_Pp1_Qp1.Yss{2});

%%
% Load GFL data
Yss_Gfl_P0_Q0 = load('Yss_Gfl_P0_Q0');
Yss_Gfl_Pp1_Q0 = load('Yss_Gfl_Pp1_Q0');
Yss_Gfl_P0_Qp1 = load('Yss_Gfl_P0_Qp1');

Ysym_Gfl_P0_Q0 = SimplusGT.ss2sym(Yss_Gfl_P0_Q0.Yss{2});
Ysym_Gfl_Pp1_Q0 = SimplusGT.ss2sym(Yss_Gfl_Pp1_Q0.Yss{2});
Ysym_Gfl_P0_Qp1 = SimplusGT.ss2sym(Yss_Gfl_P0_Qp1.Yss{2});

%%
% Plot

% Set omega
omega_p = logspace(-1,4,1e3)*2*pi;
omega_pn = [-flip(omega_p),omega_p];

%% GFM active power analysis
figure(1)
SimplusGT.bode_c(Ysym_Gfm_P0_Q0(1,1),1j*omega_p,'PhaseOn',1);
SimplusGT.bode_c(Ysym_Gfm_Pp05_Q0(1,1),1j*omega_p,'PhaseOn',1); 
SimplusGT.bode_c(Ysym_Gfm_Pp1_Q0(1,1),1j*omega_p,'PhaseOn',1); 

SimplusGT.bode_c(Ysym_Gfm_Pn05_Q0(1,1),1j*omega_p,'PhaseOn',1); 
SimplusGT.bode_c(Ysym_Gfm_Pn1_Q0(1,1),1j*omega_p,'PhaseOn',1); 

legend('P=0,Q=0','P=0.5,Q=0','P=1,Q=0','P=-0.5,Q=0','P=-1,Q=0')

subplot(2,1,1)
ylabel('Admittance')

subplot(2,1,2)
ylabel('Phase')
xlabel('Frequency')

%% GFM reactive power analysis
figure(2)
SimplusGT.bode_c(Ysym_Gfm_P0_Q0(1,1),1j*omega_p,'PhaseOn',1);
SimplusGT.bode_c(Ysym_Gfm_P0_Qp05(1,1),1j*omega_p,'PhaseOn',1); 
SimplusGT.bode_c(Ysym_Gfm_P0_Qp1(1,1),1j*omega_p,'PhaseOn',1); 

SimplusGT.bode_c(Ysym_Gfm_P0_Qn05(1,1),1j*omega_p,'PhaseOn',1); 
SimplusGT.bode_c(Ysym_Gfm_P0_Qn1(1,1),1j*omega_p,'PhaseOn',1); 

legend('P=0,Q=0','P=0,Q=0.5','P=0,Q=1','P=0,Q=-0.5','P=0,Q=-1')

subplot(2,1,1)
ylabel('Admittance')

subplot(2,1,2)
ylabel('Phase')
xlabel('Frequency')

%% GFM combined analysis
figure(3)
SimplusGT.bode_c(Ysym_Gfm_P0_Q0(1,1),1j*omega_p,'PhaseOn',1);
SimplusGT.bode_c(Ysym_Gfm_Pp1_Q0(1,1),1j*omega_p,'PhaseOn',1); 
SimplusGT.bode_c(Ysym_Gfm_P0_Qp1(1,1),1j*omega_p,'PhaseOn',1); 
SimplusGT.bode_c(Ysym_Gfm_Pp1_Qp1(1,1),1j*omega_p,'PhaseOn',1); 

legend('P=0,Q=0','P=1,Q=0','P=0,Q=1','P=1,Q=1')

subplot(2,1,1)
ylabel('Admittance')

subplot(2,1,2)
ylabel('Phase')
xlabel('Frequency')

%% GFL active power analysis
figure(4)
SimplusGT.bode_c(Ysym_Gfl_P0_Q0(1,1),1j*omega_p,'PhaseOn',1);
SimplusGT.bode_c(Ysym_Gfl_Pp1_Q0(1,1),1j*omega_p,'PhaseOn',1);
SimplusGT.bode_c(Ysym_Gfl_P0_Qp1(1,1),1j*omega_p,'PhaseOn',1);

legend('P=0,Q=0','P=1,Q=0','P=0,Q=1')

subplot(2,1,1)
ylabel('Admittance')

subplot(2,1,2)
ylabel('Phase')
xlabel('Frequency')