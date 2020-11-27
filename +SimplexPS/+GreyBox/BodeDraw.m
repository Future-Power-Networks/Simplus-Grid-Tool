function BodeDraw(DeviceSel, AxisSel, GminSS, DeviceType, N_Bus)
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
close(figure(10+i));
figure(10+AxisSel(i));
pin = 1;   %pointer to input
pout = 1;  %pointer to output
CountLegend = 0;
VecLegend = {};
P=bodeoptions;
P.Grid='on';
P.XLim={[0.1 1e4]};
P.FreqUnits='Hz';
for k = 1:N_Bus
        if DeviceType{k} <= 89 %devices
            if (ismember(k,DeviceSel)) %if selected
                bode(GminSS(pout+poutBias, pin+pinBias),P)
                CountLegend = CountLegend + 1;
                VecLegend{CountLegend} = ['Device',num2str(k)];
                hold on;
            end
            pin = pin + 4;    %4 inputs and 5 outputs.
            pout = pout + 5;
        else % floating bus, infinite bus...
            pin = pin + 2;
            pout = pout + 2;
        end
end
legend(VecLegend);
mtit(title);
end

end