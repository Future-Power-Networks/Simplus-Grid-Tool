% Impedance participation analysis of layers 1 and 2.
%
% Author(s): Yue Zhu
% Modified by: Yitong Li, Qipeng Zheng

function [Layer1, Layer2] = MdLayer12(Residue,ZmVal,N_Apparatus,ApparatusBus, ApparatusType,ApparatusSel)

for k = 1:N_Apparatus
    if (ApparatusType{k} <= 89) || (ApparatusType{k} >= 1000 && ApparatusType{k} <= 1089) || (ApparatusType{k} >= 2000 && ApparatusType{k} <= 2089)  % only consider apparatus
        % k
        % Residue{k}
        % ZmVal{k}
        Layer1All(k) = norm(Residue{k},"fro") * norm(ZmVal{k},"fro");
        % conj(sum(dot(A,B'))) = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(2,1)*B(1,2) + A(2,2)*B(2,2)
        Layer2All(k) = -1 * conj(sum(dot(Residue{k},ZmVal{k}')));
   end
end

%%
Count=0;
Layer2RealSum=0;
Layer2ImagSum=0;

for k = 1:N_Apparatus
   if (ismember(k,ApparatusSel)) %if selected 
       Count = Count + 1;
       Layer2RealSum = Layer2RealSum + abs(real(Layer2All(k)));
       Layer2ImagSum = Layer2ImagSum + abs(imag(Layer2All(k)));
   end
end

Count=0;
for k = 1:N_Apparatus
   if (ismember(k,ApparatusSel)) %if selected 
       Count = Count + 1;
       Layer1(Count) = Layer1All(k);%/sum(Layer1All);
       Layer2.real(Count) = real(Layer2All(k));
       Layer2.imag(Count) = imag(Layer2All(k));
       Layer2.real_pu(Count) = real(Layer2All(k))/Layer2RealSum;
       Layer2.imag_pu(Count) = imag(Layer2All(k))/Layer2ImagSum;
       VecLegend{Count} = ['Apparatus',num2str(ApparatusBus{k}(1,1))];
       c(Count) = categorical({['Apparatus',num2str(ApparatusBus{k}(1,1))]});
   end
end

end
