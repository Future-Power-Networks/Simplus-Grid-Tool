function BodeDraw(DeviceSel, AxisSel, GminSS, DeviceType, N_Bus, DeviceInputStr, DeviceOutputStr)
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
% P.PhaseWrapping='on';
% P.PhaseWrappingBranch=-90;
for k = 1:N_Bus
        if DeviceType{k} <= 89 %devices
            if (ismember(k,DeviceSel)) %if selected
                bode(GminSS(pout+poutBias, pin+pinBias),P)
                CountLegend = CountLegend + 1;
                VecLegend{CountLegend} = ['Node',num2str(k)];
                hold on;
            end
        end
        pin = pin + length(DeviceInputStr{k});
        pout = pout + length(DeviceOutputStr{k});
end
legend(VecLegend);
SimplexPS.mtit(title);
end

end