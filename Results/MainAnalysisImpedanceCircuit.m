% Author(s): Yitong Li

%%
clear all
clc
close all

%%
s = sym('s');       % Laplace operator

Wbase = 2*pi*50;    % (rad/s)
Vbase = 1;
Sbase = 1;
Ibase = Sbase/Vbase;
Zbase = Vbase/Ibase;
Ybase = 1/Zbase;

omega_p = logspace(-2,4,1e3)*2*pi;
omega_pn = [-flip(omega_p),omega_p];

EnablePlot = 1;

%% Droop grid-forming inverter

I(1) = 0;
I(2) = 1+1i*1;
% I(3) = 0+1i*1;

T = [1,1i;1,-1i];

for i = 1:length(I)
    [Zpn_FD{i},Zpn{i},Ypn{i}] = CalcSwingDroop('I',I(i));
    
    if EnablePlot
    figure(101)
    SimplusGT.bode_c(Ypn{i}(1,1),1j*omega_pn,'PhaseOn',1);
    figure(102)
    SimplusGT.bode_c(Ypn{i}(1,2),1j*omega_pn,'PhaseOn',1);
    figure(103)
    SimplusGT.bode_c(Zpn{i}(1,1),1j*omega_pn,'PhaseOn',1); 
    figure(104)
    SimplusGT.bode_c(Zpn{i}(1,2),1j*omega_pn,'PhaseOn',1); 
    end
    
    Zdq{i} = inv(T)*Zpn{i}*T;
    Ydq{i} = inv(T)*Ypn{i}*T;
    
    if EnablePlot
    figure(201)
    SimplusGT.bode_c(Ydq{i}(1,1),1j*omega_p,'PhaseOn',1);
    figure(202)
    SimplusGT.bode_c(Ydq{i}(1,2),1j*omega_p,'PhaseOn',1);
    figure(203)
    SimplusGT.bode_c(Ydq{i}(2,1),1j*omega_p,'PhaseOn',1);
    figure(204)
    SimplusGT.bode_c(Ydq{i}(2,2),1j*omega_p,'PhaseOn',1);
    
    figure(205)
    SimplusGT.bode_c(Zdq{i}(1,1),1j*omega_p,'PhaseOn',1);
    figure(206)
    SimplusGT.bode_c(Zdq{i}(1,2),1j*omega_p,'PhaseOn',1);
    figure(207)
    SimplusGT.bode_c(Zdq{i}(2,1),1j*omega_p,'PhaseOn',1);
    figure(208)
    SimplusGT.bode_c(Zdq{i}(2,2),1j*omega_p,'PhaseOn',1);
    end
    
end

if EnablePlot
figure(101)
legend('P=0,Q=0','P=1,Q=1')
SimplusGT.mtit('Y_{dq+-}(1,1)');
figure(102)
SimplusGT.mtit('Y_{dq+-}(1,2)');
figure(103)
SimplusGT.mtit('Z_{dq+-}(1,1)');
figure(104)
SimplusGT.mtit('Z_{dq+-}(1,2)');

figure(201)
legend('P=0,Q=0','P=1,Q=1')
SimplusGT.mtit('Y_{dq}(1,1)');
figure(202)
SimplusGT.mtit('Y_{dq}(1,2)');
figure(203)
SimplusGT.mtit('Y_{dq}(2,1)');
figure(204)
SimplusGT.mtit('Y_{dq}(2,2)');

figure(205)
legend('P=0,Q=0','P=1,Q=1')
SimplusGT.mtit('Z_{dq}(1,1)');
figure(206)
SimplusGT.mtit('Z_{dq}(1,2)');
figure(207)
SimplusGT.mtit('Z_{dq}(2,1)');
figure(208)
SimplusGT.mtit('Z_{dq}(2,2)');
end

%%
w0 = 2*pi*0.1;
% w0 = 0;
Num_Ypn = subs(Ypn{1},'s',1i*w0);
Num_Zpn = subs(Zpn{1},'s',1i*w0);

Num_Ydq = subs(Ydq{1},'s',1i*w0);
Num_Zdq = subs(Zdq{1},'s',1i*w0);

Num_Zpn_FD = subs(Zpn_FD{1},'s',1i*w0)

[~,S_Ypn,~] = svd(Num_Ypn)
[~,S_Zpn,~] = svd(Num_Zpn)

[~,S_Ypn_Abs,~] = svd(abs(Num_Ypn))
[~,S_Zpn_Abs,~] = svd(abs(Num_Zpn))

[~,S_Zpn_FD,~] = svd(Num_Zpn_FD)

% S_Ypn_tot = S_Ypn(1,1)*S_Ypn(2,2)
% S_Zpn_tot = S_Zpn(1,1)*S_Zpn(2,2)

% [~,S_Ydq,~] = svd(Num_Ydq)
% [~,S_Zdq,~] = svd(Num_Zdq)


% [~,S_Y_,~] = svd(Ypn{1})
% [~,S_Z_,~] = svd(Zpn{1})
