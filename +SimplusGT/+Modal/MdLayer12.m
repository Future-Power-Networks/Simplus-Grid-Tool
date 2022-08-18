% Author(s): Yue Zhu

% Modified by Yitong:
% Debug the imag figure title.

function [Layer1, Layer2] = MdLayer12(Residue,ZmVal,N_Bus,ApparatusType,modei,ApparatusSel,FreqSel,ModeSel)

for k = 1:N_Bus
    if ApparatusType{k} <= 89  %only consider apparatus
        %Modal layer 1
        Layer1All(k) = sqrt( Residue(k).dd*conj(Residue(k).dd) + Residue(k).dq*conj(Residue(k).dq)...
            +Residue(k).qd*conj(Residue(k).qd) +Residue(k).qq*conj(Residue(k).qq) )...
            * sqrt( ZmVal(k).dd*conj(ZmVal(k).dd) + ZmVal(k).dq*conj(ZmVal(k).dq)...
            + ZmVal(k).qd*conj(ZmVal(k).qd) + ZmVal(k).qq*conj(ZmVal(k).qq) );
        %Modal layer 2
        Layer2All(k) = -1 * ( Residue(k).dd*ZmVal(k).dd + Residue(k).qd*ZmVal(k).dq ...
                    + Residue(k).dq*ZmVal(k).qd + Residue(k).qq*ZmVal(k).qq ) ;       
   end
end

%%
%%diagrams drawing
close(figure(2010+modei));
h=figure(2010+modei);
set(h,'position',[254.6,151.4,976.8,556]);
Count=0;
Layer2RealSum=0;
Layer2ImagSum=0;

for k = 1:N_Bus
   if (ismember(k,ApparatusSel)) %if selected 
       Count = Count + 1;
       Layer2RealSum = Layer2RealSum + abs(real(Layer2All(k)));
       Layer2ImagSum = Layer2ImagSum + abs(imag(Layer2All(k)));
   end
end

Count=0;
for k = 1:N_Bus
   if (ismember(k,ApparatusSel)) %if selected 
       Count = Count + 1;
       Layer1(Count) = Layer1All(k);
       Layer2.real(Count) = real(Layer2All(k));
       Layer2.imag(Count) = imag(Layer2All(k));
       Layer2.real_pu(Count) = real(Layer2All(k))/Layer2RealSum;
       Layer2.imag_pu(Count) = imag(Layer2All(k))/Layer2ImagSum;
       VecLegend{Count} = ['Apparatus',num2str(k)];
       c(Count) = categorical({['Apparatus',num2str(k)]});
   end
end
clear title
subplot(2,2,[1,3]);
pie(Layer1);
title ('Impedance Participation Level-1');
legend(VecLegend,'Location',[0.0486,0.121,0.1038,0.2437]);

subplot(2,2,2);
b=bar(c, Layer2.real_pu);
set(gca,'YLim',[-1,1]);
set(gca,'YTick',-1:0.5:1);
set(gca,'XTickLabel',[]);
title ('Impedance Participation Level-2 Real (Normalized to 1)');
% for i=1:Count
%     text(i-0.4,Layer2.real(i),num2str(Layer2.real(i)));
% end

subplot(2,2,4);
b=bar(c, Layer2.imag_pu);
set(gca,'YLim',[-1,1])
set(gca,'YTick',-1:0.5:1);
b.FaceColor = 'flat';
b.CData = [1,0.5,0];
title ('Impedance Participation Level-2 Imag (Normalized to 1)');
% for i=1:Count
%     text(i-0.4,Layer2.imag(i),num2str(Layer2.imag(i)));
% end

TitStr = ['Mode: ',num2str(ModeSel,'%.2f'), ' Hz'];
SimplusGT.mtit(TitStr, 'color',[1 0 0], 'xoff', -0.3);

end
