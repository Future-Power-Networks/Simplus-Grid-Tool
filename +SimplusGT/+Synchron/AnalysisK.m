% KH analysis
KH_r = KH(1,:);
KH_c = KH(:,1);
KH_Vec = KH_c.*(KH_r');
[KH_min,KH_min_Index] = min(abs(KH_Vec));

% Print xi
xi_min_Hz = xi_min/2/pi
xi_min_index

% Right eigenvector
phi_r_xi = phi(:,xi_min_index);
phi_r_xi = phi_r_xi(Order_New2Old_NoFbus,1);
PhiRightMedian = median(phi_r_xi);
PhiRightMean = mean(phi_r_xi);
PhiRightPositive = find(phi_r_xi>=PhiRightMean);
PhiRightNegative = find(phi_r_xi<PhiRightMean);

% Left eigenvector
phi_l_xi = transpose(phi_inv(xi_min_index,:));
phi_l_xi = phi_l_xi(Order_New2Old_NoFbus,1);
PhiLeftPositive = find(phi_l_xi>=0);
PhiLeftNegative = find(phi_l_xi<0);

% Fiedler vector
FiedlerVec = phi_r_xi.*abs(phi_l_xi);
FiedlerMedian = median(FiedlerVec);
FiedlerAverage = mean(FiedlerVec);
FiedlerPositive = find(FiedlerVec>=FiedlerAverage);
FiedlerNegative = find(FiedlerVec<FiedlerAverage);
[~,FiedlerMinIndex] = min(FiedlerVec);
FiedlerAbsVec = abs(FiedlerVec);
[FiedlerAbsMax,FiedlerAbsMaxIndex] = max(FiedlerAbsVec);

Enable_FiedlerAbs = 1;
if Enable_FiedlerAbs
    [FiedlerMaxN,FiedlerMaxNIndex] = maxk(FiedlerAbsVec,6);
else
    [FiedlerMaxN,FiedlerMaxNIndex] = maxk(FiedlerVec,6);
end

% Deal with the floating bus
for i = 1:length(FiedlerMaxNIndex)
    if FiedlerMaxNIndex(i) >= 50
        FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 17;
    elseif FiedlerMaxNIndex(i) >= 49
        FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 15;
    elseif FiedlerMaxNIndex(i) >= 46
        FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 13;
   	elseif FiedlerMaxNIndex(i) >= 44
        FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 11;
  	elseif FiedlerMaxNIndex(i) >= 34
      	FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 10;
   	elseif FiedlerMaxNIndex(i) >= 30
      	FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 9;
   	elseif FiedlerMaxNIndex(i) >= 29
      	FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 7;
   	elseif FiedlerMaxNIndex(i) >= 28
      	FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 5;
  	elseif FiedlerMaxNIndex(i) >= 21
      	FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 2;
   	elseif FiedlerMaxNIndex(i) >= 19
      	FiedlerMaxNIndex(i) = FiedlerMaxNIndex(i) + 1;
    end
end
FiedlerMaxN;
FiedlerMaxNIndex;

% Notes:
%
% We analyze the current node and voltage node seperately because their
% transfer functions T are different. By asuming the current node is much
% "faster" than the voltage node, we can analyze the current node first
% with assume no voltage node. Then, we can analyze the voltage node by
% adding current node back into K and gamma of the voltage node.
%
% If the power system is in the loop topology. KH and Gamma_Hphi should be
% semi-diagonalized and should have some relation. The values of KH will
% rotate in each row.
%
% If analyzing voltage and current nodes seperately, KH and Gamma_Hphi for
% two subloops should be re-calculated based on K, Hinv, Gamma, and can not
% be seperate directly from the whole system KH and Gamma_Hphi.