ZminSS = SimplusGT.WholeSysZ_cal(GmObj,YbusObj,Port_i, Port_v);
%load('ZminSS_tuned','ZminSS');
figure(1);
clf;
VecLegend = {};
P=bodeoptions;
P.Grid='on';
P.XLim={[0.1 1e4]};
P.FreqUnits='Hz';
P.PhaseWrapping='on';
P.PhaseWrappingBranch=-180;
bus_sel=[1,2,3,6,8,11,12,13];

CountLegend=1;
for i=1:length(bus_sel)
    bode(ZminSS(bus_sel(i)*2-1,bus_sel(i)*2-1),P);
    k=num2str(ApparatusBus{bus_sel(i)});
    VecLegend{i} = ['$\hat{Z}_{',num2str(k),' ',num2str(k),'}$'];
    hold on;
end

legend(VecLegend,'interpreter','latex');
% 
% for i=1:length(bus_sel)
%     legend(VecLegend{7},'interpreter','latex');
% end

%save('ZminSS_tuned','ZminSS')

figure(2)
clf
Z1=load('ZminSS_detuned','ZminSS');
Z2=load('ZminSS_tuned','ZminSS');
Zdetuned=Z1.ZminSS;
Ztuned=Z2.ZminSS;

%bode(Zdetuned(21,21),P);
%hold on;
bode(Zdetuned(23,23),P);
hold on;
%bode(Zdetuned(25,25),P);
%hold on;
%bode(Ztuned(21,21),P);
%hold on;
bode(Ztuned(23,23),P);
hold on;
%bode(Ztuned(25,25),P);
%hold on;
title('Wholesystem Impedance model node-12 dd bode before and after tuning')

%ZminSS2=load('ZminSS_tuned','ZminSS');

% %P.PhaseWrapping='on';
% %P.PhaseWrappingBranch=-90;
% for k = 1:N_Apparatus
%     %[~,AppIndex] = SimplusGT.CellFind(ApparatusBus,k);
%     if ApparatusType{k} <= 89 %apparatuses
%         if (ismember(k,ApparatusSel)) %if selected
%             bode(GminSS(pout+poutBias, pin+pinBias),P)
%             CountLegend = CountLegend + 1;
%             VecLegend{CountLegend} = ['Node',num2str(ApparatusBus{k})];
%             hold on;
%         end
%         pin = pin + length(ApparatusInputStr{k});
%         pout = pout + length(ApparatusOutputStr{k});
%     end
% end
% legend(VecLegend);
% SimplusGT.mtit(title);
% end
% 
% end