function [Layer1, Layer2] = Layer12(DeviceSel,...
     DeviceType, N_Bus, GmDSS_Cell, Mode, FreqSel, modei, Phi, Psi)
%This function gives results of Greybox Layer1 and Layer2, also draw pie 
%charts and bar charts.
%Author: Yue Zhu.

%% Greybox Layer 1&2 calculation
pin = 1;   %pointer to input
pout = 1;  %pointer to output
for k = 1:N_Bus
    if DeviceType{k} <= 89  %apparatus
        %calculate residues.
        Residue(k).dd=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
        Residue(k).dq=C(pout,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1);
        Residue(k).qd=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin);
        Residue(k).qq=C(pout+1,:) * Phi(:,ModeSel) * Psi(ModeSel,:) * B(:,pin+1);
        %calculate device impedance.
        GmSS_Cell{k} = minreal(GmDSS_Cell{k});
        GmTf(k).dd=tf(GmSS_Cell{k}(1,1));
        GmTf(k).dq=tf(GmSS_Cell{k}(1,2));
        GmTf(k).qd=tf(GmSS_Cell{k}(2,1));
        GmTf(k).qq=tf(GmSS_Cell{k}(2,2));
        Gm{k} = [GmTf(k).dd, GmTf(k).dq ; GmTf(k).qd, GmTf(k).qq];
        Zm{k} = inv(Gm{k});
        %calculate the value of impedance
        1i;
        ZmVal(k).dd = evalfr( Zm{k}(1,1), 2*pi*FreqSel*1i);
        ZmVal(k).dq = evalfr( Zm{k}(1,2), 2*pi*FreqSel*1i);
        ZmVal(k).qd = evalfr( Zm{k}(2,1), 2*pi*FreqSel*1i);
        ZmVal(k).qq = evalfr( Zm{k}(2,2), 2*pi*FreqSel*1i);
        %Greybox layer 1
        Layer1All(k) = sqrt( Residue(k).dd*conj(Residue(k).dd) + Residue(k).dq*conj(Residue(k).dq)...
            +Residue(k).qd*conj(Residue(k).qd) +Residue(k).qq*conj(Residue(k).qq) )...
            * sqrt( ZmVal(k).dd*conj(ZmVal(k).dd) + ZmVal(k).dq*conj(ZmVal(k).dq)...
            + ZmVal(k).qd*conj(ZmVal(k).qd) + ZmVal(k).qq*conj(ZmVal(k).qq) );
        %Greybox layer 2
        Layer2All(k) = -1 * ( Residue(k).dd*ZmVal(k).dd + Residue(k).qd*ZmVal(k).dq ...
                    + Residue(k).dq*ZmVal(k).qd + Residue(k).qq*ZmVal(k).qq ) ;        
        pin = pin + 4;    %4 inputs and 5 outputs.
        pout = pout + 5;
    else %passive load.
        pin = pin + 2;
        pout = pout + 2;
   end
end

%%
%%plot draw
close(figure(14+modei));
figure(14+modei)
Count=0;
for k = 1:N_Bus
   if (ismember(k,DeviceSel)) %if selected 
       Count = Count + 1;
       Layer1(Count) = Layer1All(k);
       Layer2.real(Count) = real(Layer2All(k));
       Layer2.imag(Count) = imag(Layer2All(k));
       VecLegend{Count} = ['Device',num2str(k)];
       c(Count) = categorical({['Device',num2str(k)]});
   end
end
clear title
subplot(2,2,[1,3]);
pie(Layer1);
title ('Greybos Level-1');
legend(VecLegend,'Location','southwest');

subplot(2,2,2);
b=bar(c, Layer2.real);
title ('Greybos Level-2 Real');
for i=1:Count
    text(i-0.4,Layer2.real(i),num2str(Layer2.real(i)));
end

subplot(2,2,4);
b=bar(c, Layer2.imag);
b.FaceColor = 'flat';
b.CData = [1,0.5,0];
title ('Greybos Level-2 Imag');
for i=1:Count
    text(i-0.4,Layer2.imag(i),num2str(Layer2.imag(i)));
end

TitStr = ['Mode: ',num2str(FreqSel), ' Hz'];
mtit(TitStr, 'color',[1 0 0], 'xoff', -0.3);

end
