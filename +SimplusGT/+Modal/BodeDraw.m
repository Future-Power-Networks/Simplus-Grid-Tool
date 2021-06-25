function BodeDraw(ApparatusSel, AxisSel, GminSS, ApparatusType, ApparatusBus, N_Bus, ApparatusInputStr, ApparatusOutputStr)
%this function draws the bodeplot of the node admittance you selected
%Author: Yue Zhu.
AxisNum=length(AxisSel);

for i=1:AxisNum
    if AxisSel(i)==1 %dd axis;        
        pinBias=0;
        poutBias=0;
        title='Admittance Bode Diagram: d-d axis';
    elseif AxisSel(i)==2 %dq axis;
        pinBias=1;
        poutBias=0;
        title='Admittance Bode Diagram: d-q axis';
    elseif AxisSel(i)==3 %qd axis;
        pinBias=0;
        poutBias=1;
        title='Admittance Bode Diagram: q-d axis';
    elseif AxisSel(i)==4 %qq axis;
        pinBias=1;
        poutBias=1;
        title='Admittance Bode Diagram: q-q axis';
    else
    end  
close(figure(2000+AxisSel(i)));
figure(2000+AxisSel(i));
pin = 1;   %pointer to input
pout = 1;  %pointer to output
CountLegend = 0;
VecLegend = {};
P=bodeoptions;
P.Grid='on';
P.XLim={[0.1 1e4]};
P.FreqUnits='Hz';
%P.PhaseWrapping='on';
%P.PhaseWrappingBranch=-90;
for k = 1:N_Bus
    [~,AppIndex] = SimplusGT.CellFind(ApparatusBus,k);
    if ApparatusType{AppIndex} <= 89 %apparatuses
        if (ismember(k,ApparatusSel)) %if selected
            bode(GminSS(pout+poutBias, pin+pinBias),P)
            CountLegend = CountLegend + 1;
            VecLegend{CountLegend} = ['Node',num2str(k)];
            hold on;
        end
        pin = pin + length(ApparatusInputStr{k});
        pout = pout + length(ApparatusOutputStr{k});
    end
end
legend(VecLegend);
SimplusGT.mtit(title);
end

end