f_low=0.5; % analysed frequency point

%% Calculation of Thevenin impedance
for i=1:length(ApparatusType)
    %Ym_matrix(i*2-1:i*2, i*2-1:i*2) = GmDSS_Cell{i}(1:2,1:2);
    if ApparatusType{i}==100 || ApparatusType{i}==90 % floating bus: machine impedance is set as 0
        Zm0(i*2-1:i*2, i*2-1:i*2) = 0;
        Zm0_svd(i,1:2) = 0;
    else
        Zm0(i*2-1:i*2, i*2-1:i*2) = inv(evalfr(GmDSS_Cell{i}(1:2,1:2), 1i*2*pi*f_low));
        Zm0_svd(i,1:2) = svd(Zm0(i*2-1:i*2, i*2-1:i*2)).';
    end
end

%Ybus0 = evalfr(YbusDSS, 1i*2*pi*f_low);
%Ysys_p = (eye(N_Bus)+Ybus0_p*Zm0_mat_p) \ Ybus0_p ;
%Ysys0 = inv(eye(2*N_Bus)+Ybus0*Zm0) * Ybus0 ;
Ysys0 = evalfr(YsysDSS, 1i*2*pi*f_low);

clear LIMR;
j=1;

for i=1:N_Bus
    if ApparatusType{i}<30 || ApparatusType{i}>=40
        if ApparatusType{i}==100 || ApparatusType{i}==90% floating bus or infinite bus
            %Zth{i} = 'NaN';%inv(Ysys0(i*2-1:i*2, i*2-1:i*2));
            %Zth_(i) = 1/(Ysys_p(i,i));
            LIMR(j).device = i;
            LIMR(j).value = 1;
            j=j+1;
        else
            LIMR(j).Zm_norm = norm(Zm0(i*2-1:i*2, i*2-1:i*2),"fro");
            LIMR(j).Zth_norm = norm( inv(Ysys0(i*2-1:i*2, i*2-1:i*2)) - Zm0(i*2-1:i*2, i*2-1:i*2), "fro");
            LIMR(j).device = i;
            LIMR(j).value = (LIMR(j).Zm_norm-LIMR(j).Zth_norm)/LIMR(j).Zm_norm;
            j=j+1;
        end
    end
    %Zth_svd(i,:) = svd(Zth{i}).';
    %SCC(i) = 1/Zth_svd(i,2);
end
%Main_K_Plot(LIMR,2);
%title('Large-Signal System Strength Heatmap');