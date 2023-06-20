% input: 
% Ynodal0: numerical matrix of the nodal admittance matrix.
% ListBus: from Power flow
% ApparatusType, Numb bus...
% output:
% HNM: hybrid nodal matrix
% OrderOld2New: re-ordered buses.

% this is for positive sequence calculation

function [HNM, OrderOld2New, ApparatusSourceType] = HybridMatrix_ps(Ynodal0, ListBus, ApparatusType, NumBus)
[OrderOld2New, ApparatusSourceType]=ReorderBus(ListBus,ApparatusType,NumBus);
OrderOld2New_p = zeros(1,length(OrderOld2New));
for i=1:length(OrderOld2New)
    OrderOld2New_p(i) = OrderOld2New(i);
    OrderOld2New_p(i) = OrderOld2New(i);
end
Ynodal0_re = Ynodal0(OrderOld2New_p, OrderOld2New_p); % re-ordered Ynodal: Vsource, Isource, Floating.

% devide the new Ynodal into four matrix blocks.
vs_num=length(find(ApparatusSourceType==1));
Yp11 = Ynodal0_re(1:vs_num, 1:vs_num);
Yp12 = Ynodal0_re(1:vs_num, (vs_num+1):end);
Yp21 = Ynodal0_re((vs_num+1):end, 1:vs_num);
Yp22 = Ynodal0_re((vs_num+1):end, (vs_num+1):end);

% get the hybrid matrix. Attention: the order of the buses is changed now
HNM = [Yp11 - Yp12*inv(Yp22)*Yp21, Yp12*inv(Yp22);...
            -inv(Yp22)*Yp21,            inv(Yp22)];

end