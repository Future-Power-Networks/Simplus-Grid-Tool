% This function analyzes the stability by classfying voltage and current
% nodes respectively.

% ### Voltage node
if Exist_Vbus == 1
if isempty(KH22)
    KH_V = KH;
    [phi_V,xi_V] = eig(KH_V);
    GammaHphi_V = GammaHphi;
    [~,sigma_V,~] = svd(GammaHphi_V);
else
    K_V = K11 - K12*inv(K22)*K21;
    Gamma_V =  Gamma11 - K12*inv(K22)*Gamma21;
    Hinv_V = Hinv(1:(n_Ibus_1st-1),1:(n_Ibus_1st-1));
    KH_V = Hinv_V*K_V;
    [phi_V,xi_V] = eig(KH_V);
    GammaHphi_V = inv(phi_V)*Hinv_V*Gamma_V*phi_V;
    [~,sigma_V,~] = svd(GammaHphi_V);
end
sigma_V_max = max(max(sigma_V));
[zeta_m_V,w_min_V] = CalcZeta(T_V_sym{n_v_ref},diag(xi_V));
if min(min(real(xi_V)))<-1e-4
    fprintf(['Warning: xi_V_min = ' num2str(min(min(real(xi_V)))) ' < 0.']);
end
else
	xi_V = [];
    sigma_V = [];
    zeta_m_V = [];
    w_min_V = [];
end

% ### Current node 
if Exist_Ibus == 1
    % KH_I = K22;
    % Gamma_Hphi_I = GH22;
    K_I = K22;
    Gamma_I = Gamma22;
    Hinv_I = Hinv(n_Ibus_1st:N_Bus,n_Ibus_1st:N_Bus);
    KH_I = Hinv_I*K_I;
    [phi_I,xi_I] = eig(KH_I);
    GammaHphi_I = inv(phi_I)*Hinv_I*Gamma_I*phi_I;
    [~,sigma_I,~] = svd(GammaHphi_I);
    sigma_I_max = max(max(sigma_I));
    [zeta_m_I,w_min_I] = CalcZeta(T_I_sym{n_i_ref},diag(xi_I));
  	if min(min(real(xi_I)))<-1e-4
        fprintf(['Warning: xi_I_min = ' num2str(min(min(real(xi_I)))) ' < 0.\n']);
    end
else
    xi_I = [];
    sigma_I = [];
    zeta_m_I = [];
    w_min_I = [];
end