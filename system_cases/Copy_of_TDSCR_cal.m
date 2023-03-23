%% Traditional SCR: from Thevenin impedance
% Consider all the machines are working at its rated power

YN=evalfr(YbusDSS,0); % nodal admittance matrix
ZN=inv(YN); % nodal impedance (diagonal elements equal Thevenin impedance)
%ObjYAss=SimplusGT.ObjDss2Ss(ObjYmDSS);
%[~,YA_s]= ObjYAss.GetSS(ObjYASS);
%YA=evalfr(tf(YA_s),0); % apparatus admittance matrix

% Rated MVA calculation: % original load - load from PF.
PG=ListBus_(:,7) - ListPowerFlow(:,2);
QG=ListBus_(:,8) - ListPowerFlow(:,3);
Prat = sqrt(PG.^2+QG.^2);

SCR=zeros(1,N_Bus);
for i=1:N_Bus
    if ApparatusType{i}~=100 %not a load bus
        SCR(i)=1/ norm(ZN(i*2-1:i*2, i*2-1:i*2)) / Prat(i) ;
    else
        SCR(i)=NaN;
    end
end
SCR

%% TD-SCR from H inverse --- new
[HN, OrderOld2New, ApparatusSourceType] = HybridMatrix(YN, ListBus, ApparatusType, N_Bus);
Prat_h = Prat(OrderOld2New); % rated power in new order
GN=inv(HN);
TDscr_x=zeros(1,N_Bus);
TDescr_x=zeros(1,N_Bus);
TDIF=zeros(N_Bus,N_Bus);
TDIF_sum=zeros(1,N_Bus);
for k=1:N_Bus
    TDIF(k,k)=1;
    if ApparatusSourceType(k)==1 % if bus-k connects a voltage source
        for i=1:N_Bus
                if ApparatusSourceType(i)==1 && i~=k % if bus-i connects a voltage source
                    TDIF(k,i) = norm( GN(2*k-1:2*k,2*i-1:2*i) / GN(2*k-1:2*k,2*k-1:2*k) ) ;               
                elseif ApparatusSourceType(i)==2 && i~=k % if bus-i connects a current source
                    TDIF(k,i) = norm( GN(2*i-1:2*i,2*k-1:2*k) );
                end
                TDIF_sum(k) = TDIF_sum(k) + TDIF(k,i)*Prat_h(i);
        end
        TDscr_x(k)=1/norm(GN(k*2-1:k*2,k*2-1:k*2))/Prat_h(k);
        TDescr_x(k)=1/norm(GN(k*2-1:k*2,k*2-1:k*2))/(Prat_h(k)+TDIF_sum(k));
        
        
    elseif ApparatusSourceType(k)==2 % if bus-k connects a current source       
        for i=1:N_Bus
                if ApparatusSourceType(i)==1 && i~=k % if bus-i connects a voltage source
                    TDIF(k,i) = norm( HN(2*i-1:2*i,2*k-1:2*k) );          
                elseif ApparatusSourceType(i)==2 && i~=k % if bus-i connects a current source
                    TDIF(k,i) = norm( HN(2*k-1:2*k,2*i-1:2*i) / HN(2*k-1:2*k,2*k-1:2*k) ) ; 
                end
                TDIF_sum(k) = TDIF_sum(k) + TDIF(k,i)*Prat_h(i);
        end      
        TDscr_x(k)=1/norm(HN(k*2-1:k*2,k*2-1:k*2))/Prat_h(k);
        TDescr_x(k)=1/norm(HN(k*2-1:k*2,k*2-1:k*2))/(Prat_h(k)+TDIF_sum(k));
        
        
    elseif ApparatusSourceType(k)==3 % load bus
        TDscr_x(k)=NaN; % leave it for now
        TDescr_x(k)=NaN;
    end
end

TDscr_x(OrderOld2New)=TDscr_x
TDescr_x(OrderOld2New)=TDescr_x

%% TD-ESCR
Result_scr=round([SCR;TDscr_x;TDescr_x;reshape(Prat,1,N_Bus)],3);
