% Participation analysis of K and xi.

for i = 1:N_Bus
    Hinv_PA{i} = Hinv;
    Hii = 1/Hinv(i,i);
    dH_PA = Hii*1e-3;                           % Find the perturbation value
    Hinv_PA{i}(i,i) = 1/(Hii + dH_PA);         	% Perturb the inertia of the ith device
    KH_PA{i} = Hinv_PA{i}*K;
    [phi_PA{i},xi_PA{i}] = eig(KH_PA{i});
    dxi_PA = xi_PA{i} - xi;                    	% Disturbance on the jth xi
    for j = 1:N_Bus
        ParAnalysisK(j,i) = dxi_PA(j,j)/dH_PA;	% Participation factor
    end
end