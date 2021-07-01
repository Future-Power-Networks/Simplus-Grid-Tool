% Calculate the admittance sensitivity matrix and nodal admittance matrix at the k-th labmda, 
% Final results will be numerical matrices.
% Author: Yue Zhu
function [SensMat,Ybus_val,Ynodal_val,Yre_val,ZminSS,Mode_Hz] = SensitivityCal(Ek)

GmObj_Cell=evalin('base', 'GmObj_Cell');
YbusObj=evalin('base', 'YbusObj');
N_Apparatus=evalin('base', 'N_Apparatus');
N_Bus=evalin('base', 'N_Bus');
GmDSS_Cell=evalin('base', 'GmDSS_Cell');
ApparatusBus=evalin('base', 'ApparatusBus');


[ZsysObj,ZsysDSS] = SimplusGT.WholeSysZ_cal(GmObj_Cell,YbusObj,N_Apparatus);
ZminSS = SimplusGT.dss2ss(ZsysDSS);

[~,ZmInStr,ZmOutStr] = ZsysObj.GetString(ZsysObj);%get name string
Port_i_in = [];
Port_v_out = [];
for i = 1:N_Apparatus
    [~,in1] = SimplusGT.CellFind(ZmInStr,['i_d',num2str(i)]);
    [~,in2] = SimplusGT.CellFind(ZmInStr, ['i',num2str(i)]);
    
    [~,out1] = SimplusGT.CellFind(ZmOutStr,['v_d',num2str(i)]);
    [~,out2] = SimplusGT.CellFind(ZmOutStr,['v',num2str(i)]);
    if ~isempty(in1) %ac apparatus
        Port_i_in = [Port_i_in,in1,in1+1];
        Port_v_out = [Port_v_out,out1,out1+1];
    elseif ~isempty(in2) % dc: not activated at this moment.
        Port_i_in = [Port_i_in,in2];
        Port_v_out = [Port_v_out,out2];
    else
        error(['Error']);
    end
end

A=ZminSS.A;
B=ZminSS.B;
C=ZminSS.C;
[Phi,D]=eig(A);
Psi=inv(Phi); 
Mode=diag(D);
Mode_Hz=Mode/2/pi;

%% Trim B and C, keep only i as input and v as output
for i = 1:N_Bus*2
    Btrim(:,i)=B(:,Port_i_in(i));
    Ctrim(i,:)=C(Port_v_out(i),:);
end

%% Admittance Sensitivity Matrix = -1* residue matrix
SensMat_exp = -1*Ctrim*Phi(:,Ek)*Psi(Ek,:)*Btrim; % Residue matrix in expansion form
% pack up into d-q format
for i = 1:N_Bus
    for j = 1:N_Bus
        SensMat(i,j).dd = SensMat_exp(2*i-1,2*j-1);
        SensMat(i,j).dq = SensMat_exp(2*i-1,2*j);
        SensMat(i,j).qd = SensMat_exp(2*i,2*j-1);
        SensMat(i,j).qq = SensMat_exp(2*i,2*j);
    end
end

%% Node and branch Admittance values
lambda = Mode(Ek);
[~,YbusDSS]=YbusObj.GetDSS(YbusObj);
Ybus_val_exp = evalfr(YbusDSS, lambda); % Nodal admittance matrix Value (only passive part)

Ynodal_val_exp = Ybus_val_exp; % Nodal admittance matrix (include apparatus admittance)

for i=1:N_Apparatus
    bus_i = ApparatusBus{i};
    Ynodal_val_exp(2*bus_i-1,2*bus_i-1) = Ybus_val_exp(2*bus_i-1,2*bus_i-1)+ evalfr(GmDSS_Cell{i}(1,1),lambda); %dd
    Ynodal_val_exp(2*bus_i-1,2*bus_i  ) = Ybus_val_exp(2*bus_i-1,2*bus_i)  + evalfr(GmDSS_Cell{i}(1,2),lambda); %dq
    Ynodal_val_exp(2*bus_i  ,2*bus_i-1) = Ybus_val_exp(2*bus_i,  2*bus_i-1)+ evalfr(GmDSS_Cell{i}(2,1),lambda); %qd
    Ynodal_val_exp(2*bus_i  ,2*bus_i  ) = Ybus_val_exp(2*bus_i,  2*bus_i)  + evalfr(GmDSS_Cell{i}(2,2),lambda); %qq    
end

for i = 1:N_Bus
    for j = 1:N_Bus
        Ybus_val(i,j).dd = Ybus_val_exp(2*i-1,2*j-1);
        Ybus_val(i,j).dq = Ybus_val_exp(2*i-1,2*j);
        Ybus_val(i,j).qd = Ybus_val_exp(2*i,2*j-1);
        Ybus_val(i,j).qq = Ybus_val_exp(2*i,2*j);
        
        Ynodal_val(i,j).dd = Ynodal_val_exp(2*i-1,2*j-1);
        Ynodal_val(i,j).dq = Ynodal_val_exp(2*i-1,2*j);
        Ynodal_val(i,j).qd = Ynodal_val_exp(2*i,2*j-1);
        Ynodal_val(i,j).qq = Ynodal_val_exp(2*i,2*j);        
    end
end

% Rearranged network admittance matrix
for i = 1:N_Bus
    for j=1:N_Bus
        if i == j % node element : sum of a row / column
            Yre_val(i,i).dd = 0;%Ynodal_val(i,i).dd;
            Yre_val(i,i).dq = 0;%Ynodal_val(i,i).dq;
            Yre_val(i,i).qd = 0;%Ynodal_val(i,i).qd;
            Yre_val(i,i).qq = 0;%Ynodal_val(i,i).qq;
            for k = 1:N_Bus
               Yre_val(i,i).dd = Yre_val(i,i).dd+Ynodal_val(i,k).dd;
               Yre_val(i,i).dq = Yre_val(i,i).dq+Ynodal_val(i,k).dq;
               Yre_val(i,i).qd = Yre_val(i,i).qd+Ynodal_val(i,k).qd;
               Yre_val(i,i).qq = Yre_val(i,i).qq+Ynodal_val(i,k).qq;
            end
        else % branch element : inverse of the Ynodal
            Yre_val(i,j).dd = -1*Ynodal_val(i,j).dd;
            Yre_val(i,j).dq = -1*Ynodal_val(i,j).dq;
            Yre_val(i,j).qd = -1*Ynodal_val(i,j).qd;
            Yre_val(i,j).qq = -1*Ynodal_val(i,j).qq;
        end
    end
end



end