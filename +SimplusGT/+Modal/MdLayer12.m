% Impedance participation analysis of layers 1 and 2.
%
% Author(s): Yue Zhu
% Modified by: Yitong Li, Qipeng Zheng

function [Layer1, Layer2] = MdLayer12(Residue,ZmVal,N_Apparatus,ApparatusBus, ApparatusType,modei,ApparatusSel,FreqSel,ModeSel)

for k = 1:N_Apparatus
    if ApparatusType{k} <= 89  %only consider apparatus
        Layer1All(k) = sqrt( Residue(k).dd*conj(Residue(k).dd) + Residue(k).dq*conj(Residue(k).dq)...
            +Residue(k).qd*conj(Residue(k).qd) +Residue(k).qq*conj(Residue(k).qq) )...
            * sqrt( ZmVal(k).dd*conj(ZmVal(k).dd) + ZmVal(k).dq*conj(ZmVal(k).dq)...
            + ZmVal(k).qd*conj(ZmVal(k).qd) + ZmVal(k).qq*conj(ZmVal(k).qq) );
        
        Layer2All(k) = -1 * ( Residue(k).dd*ZmVal(k).dd + Residue(k).qd*ZmVal(k).dq ...
                    + Residue(k).dq*ZmVal(k).qd + Residue(k).qq*ZmVal(k).qq ) ;       
    elseif ApparatusType{k} >= 1010 && ApparatusType{k} <= 1089 % DC apparatuses
        Layer1All(k) = sqrt( Residue(k).dd*conj(Residue(k).dd));
        Layer2All(k) = -1 * ( Residue(k).dd* ZmVal(k).dd )/sqrt( ZmVal(k).dd*conj(ZmVal(k).dd) );
    elseif ApparatusType{k} >= 2000 && ApparatusType{k} <= 2009 % IC apparatuses
        Layer1All(k) = sqrt( Residue(k).dd*conj(Residue(k).dd) + Residue(k).dq*conj(Residue(k).dq)...
            +Residue(k).qd*conj(Residue(k).qd) +Residue(k).qq*conj(Residue(k).qq) ...
            +Residue(k).d_dc*conj(Residue(k).d_dc) +Residue(k).q_dc*conj(Residue(k).q_dc)...
            +Residue(k).dc_d*conj(Residue(k).dc_d) +Residue(k).dc_q*conj(Residue(k).dc_q)...
            +Residue(k).dc_dc*conj(Residue(k).dc_dc));

        Layer2All(k) = -1 * ( Residue(k).dd*ZmVal(k).dd + Residue(k).qd*ZmVal(k).dq ...
                    + Residue(k).dq*ZmVal(k).qd + Residue(k).qq*ZmVal(k).qq...
                    + Residue(k).d_dc*ZmVal(k).dc_d + Residue(k).q_dc*ZmVal(k).dc_q...
                    + Residue(k).dc_d*ZmVal(k).d_dc + Residue(k).dc_q*ZmVal(k).q_dc...
                    + Residue(k).dc_dc*ZmVal(k).dc_dc)/...
                    sqrt( ZmVal(k).dd*conj(ZmVal(k).dd) + ZmVal(k).dq*conj(ZmVal(k).dq)...
            + ZmVal(k).qd*conj(ZmVal(k).qd) + ZmVal(k).qq*conj(ZmVal(k).qq) ...
            + ZmVal(k).d_dc*conj(ZmVal(k).d_dc) + ZmVal(k).q_dc*conj(ZmVal(k).q_dc) ...
            + ZmVal(k).dc_d*conj(ZmVal(k).dc_d) + ZmVal(k).dc_q*conj(ZmVal(k).dc_q) ...
            + ZmVal(k).dc_dc*conj(ZmVal(k).dc_dc));  

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
       VecLegend{Count} = ['Apparatus',num2str(ApparatusBus{k})];
       c(Count) = categorical({['Apparatus',num2str(ApparatusBus{k})]});
   end
end

end
