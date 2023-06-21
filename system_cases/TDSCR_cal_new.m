
f_low=0.1; % analysed frequency point

%% Traditional SCR: from Thevenin impedance
% impedance of apparatus: in positive sequency, with frequency coupling
% effect (FCE) removed
for i=1:length(ApparatusType)
    %Ym_matrix(i*2-1:i*2, i*2-1:i*2) = GmDSS_Cell{i}(1:2,1:2);
    if ApparatusType{i}>=90 % floating bus: impedance is set as 1e7, a large value
        Zm0(i*2-1:i*2, i*2-1:i*2) = 1e7;
        Zm0_svd(i,1:2) = 1e7;
    else
        Zm0(i*2-1:i*2, i*2-1:i*2) = inv(evalfr(GmDSS_Cell{i}(1:2,1:2), 1i*2*pi*f_low));
        Zm0_svd(i,1:2) = svd(Zm0(i*2-1:i*2, i*2-1:i*2)).';
    end
end
% Ybus to positive sequence value, with FCE parts removed
Ybus0 = evalfr(YbusDSS, 1i*2*pi*f_low);
%Ysys_p = (eye(N_Bus)+Ybus0_p*Zm0_mat_p) \ Ybus0_p ;
Ysys0 = inv(eye(2*N_Bus)+Ybus0*Zm0) * Ybus0 ;

for i=1:N_Bus
    if ApparatusType{i}>=90 % floating bus
        Zth{i} = inv(Ysys0(i*2-1:i*2, i*2-1:i*2));
        %Zth_(i) = 1/(Ysys_p(i,i));
    else
        Zth{i} = inv(Ysys0(i*2-1:i*2, i*2-1:i*2)) - Zm0(i*2-1:i*2, i*2-1:i*2);
    end
    Zth_svd(i,:) = svd(Zth{i}).';
    SCC(i) = 1/Zth_svd(i,2)
end
%Zth_
Zth_svd
SCC
%Zth_p

%% Type depedent
% switch in&out for nodal admittance matrix
[HN, OrderOld2New, ApparatusSourceType] = HybridMatrix(Ybus0, ListBus, ApparatusType, N_Bus);
% switch impedance - admittance for machine impedance
for i=1:N_Bus
    if ApparatusSourceType(i) == 1 % voltage type
        HM(i*2-1:i*2,i*2-1:i*2) = Zm0(i*2-1:i*2,i*2-1:i*2);
    elseif ApparatusSourceType(i) == 2 % current type
        HM(i*2-1:i*2,i*2-1:i*2) = inv(Zm0(i*2-1:i*2,i*2-1:i*2));
    elseif ApparatusSourceType(i) == 3 % floating bus: zero admittance
        HM(i*2-1:i*2,i*2-1:i*2) = 0;
    end
end
Hsys = inv(eye(2*N_Bus)+HN*HM) * HN;

for i=1:N_Bus
    if ApparatusSourceType(i) == 1 % voltage type
        TDZth{i} = inv(Hsys(i*2-1:i*2,i*2-1:i*2))-Zm0(i*2-1:i*2, i*2-1:i*2);
    elseif ApparatusSourceType(i) == 2 % current type
        TDZth{i} = inv( inv(Hsys(i*2-1:i*2,i*2-1:i*2)) - inv(Zm0(i*2-1:i*2, i*2-1:i*2)) )
    elseif ApparatusSourceType(i) == 3 % floating bus: zero admittance
        TDZth{i} = zeros(2);
    end
    TDZth_svd(i,:) = svd(TDZth{i}).';
    TDSCC(i) = 1/TDZth_svd(i,2);
end
TDZth_svd
TDSCC
% for i=1:N_Bus
%     if ApparatusSourceType(i) == 1 % voltage type
%         SCC_H(i) = 1/(1/(HN(i,i)) - Zm0_mat_p(i,i));
%     elseif ApparatusSourceType(i) == 2 % current type
%         SCC_H(i) = 1/HN(i,i) - 1/Zm0_mat_p(i,i); %1/(HN(i,i)-Zm0_mat_p(i,i));
%     elseif ApparatusSourceType(i) == 3 % floating bus
%         SCC_H(i) = 'NAN';
%     end
% end
% %HN
% SCC_H
% 
% Zm0_mat_p
% HN
% Prat_h = Prat(OrderOld2New); % rated power in new order

%Ysys0 = evalfr(YsysDSS, 1i*2*pi*f_low);
% for i=1:N_Bus
%     if ApparatusType{i}>=90 % floating bus
%         Zth{i} = inv(Ysys0(i*2-1:i*2, i*2-1:i*2));
%     else
%         Zth{i} = inv(Ysys0(i*2-1:i*2, i*2-1:i*2)) - Zm0(i*2-1:i*2, i*2-1:i*2);
%     end
%     Zth_p(i,1:2) = (svd(Zth{i}));
% end

% Consider all the machines are working at its rated power

%YN=evalfr(YbusDSS,0); % nodal admittance matrix
%ZN=inv(YN); % nodal impedance (diagonal elements equal Thevenin impedance)
%ObjYAss=SimplusGT.ObjDss2Ss(ObjYmDSS);
%[~,YA_s]= ObjYAss.GetSS(ObjYASS);
%YA=evalfr(tf(YA_s),0); % apparatus admittance matrix

% % Rated MVA calculation: % original load - load from PF.
% PG=ListBus_(:,7) - ListPowerFlow(:,2);
% QG=ListBus_(:,8) - ListPowerFlow(:,3);
% Prat = sqrt(PG.^2+QG.^2);
% 
% SCR=zeros(1,N_Bus);
% for i=1:N_Bus
%     if ApparatusType{i}~=100 %not a load bus
%         SCR(i)=1/ norm(ZN(i*2-1:i*2, i*2-1:i*2)) / Prat(i) ;
%         %SCR(i)=1/ norm(ZN(i*2-1:i*2, i*2-1:i*2)) / Prat(i) ;
%     else
%         SCR(i)=NaN;
%     end
% end
% SCR
% 
% %% TD-SCR from H inverse --- new
% [HN, OrderOld2New, ApparatusSourceType] = HybridMatrix(YN, ListBus, ApparatusType, N_Bus);
% Prat_h = Prat(OrderOld2New); % rated power in new order
% GN=inv(HN);
% TDscr_x=zeros(1,N_Bus);
% TDescr_x=zeros(1,N_Bus);
% TDIF=zeros(N_Bus,N_Bus);
% TDIF_sum=zeros(1,N_Bus);
% for k=1:N_Bus
%     TDIF(k,k)=1;
%     if ApparatusSourceType(k)==1 % if bus-k connects a voltage source
%         for i=1:N_Bus
%                 if ApparatusSourceType(i)==1 && i~=k % if bus-i connects a voltage source
%                     TDIF(k,i) = norm( GN(2*i-1:2*i,2*k-1:2*k) / GN(2*k-1:2*k,2*k-1:2*k) ) ;               
%                 elseif ApparatusSourceType(i)==2 && i~=k % if bus-i connects a current source
%                     TDIF(k,i) = norm( GN(2*k-1:2*k,2*i-1:2*i) );
%                 end
%                 TDIF_sum(k) = TDIF_sum(k) + TDIF(k,i)*Prat_h(i);
%         end
%         TDscr_x(k)=norm(HN(k*2-1:k*2,k*2-1:k*2))/Prat_h(k);
%         %TDscr_x(k)=1/norm(GN(k*2-1:k*2,k*2-1:k*2))/Prat_h(k);
%         TDescr_x(k)=1/norm(GN(k*2-1:k*2,k*2-1:k*2))/(Prat_h(k)+TDIF_sum(k));
%         
%         
%     elseif ApparatusSourceType(k)==2 % if bus-k connects a current source       
%         for i=1:N_Bus
%                 if ApparatusSourceType(i)==1 && i~=k % if bus-i connects a voltage source
%                     TDIF(k,i) = norm( HN(2*k-1:2*k,2*i-1:2*i) );          
%                 elseif ApparatusSourceType(i)==2 && i~=k % if bus-i connects a current source
%                     TDIF(k,i) = norm( HN(2*i-1:2*i,2*k-1:2*k) /  HN(2*k-1:2*k,2*k-1:2*k) ) ; 
%                 end
%                 TDIF_sum(k) = TDIF_sum(k) + TDIF(k,i)*Prat_h(i);
%         end      
%         TDscr_x(k)=1/norm(HN(k*2-1:k*2,k*2-1:k*2))/Prat_h(k);
%         TDescr_x(k)=1/norm(HN(k*2-1:k*2,k*2-1:k*2))/(Prat_h(k)+TDIF_sum(k));
%         
%         
%     elseif ApparatusSourceType(k)==3 % load bus
%         TDscr_x(k)=NaN; % leave it for now
%         TDescr_x(k)=NaN;
%     end
% end
% 
% TDscr_x(OrderOld2New)=TDscr_x
% TDescr_x(OrderOld2New)=TDescr_x
% 
% %% TD-ESCR
% Result_scr=round([SCR;TDscr_x;TDescr_x;reshape(Prat,1,N_Bus)],3);
