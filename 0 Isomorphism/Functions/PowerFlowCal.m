% ### Power flow analysis
if Flag_PowerFlowAlgorithm == 1
    [PowerFlow,~,~,~,~,~,~,~] = SimplexPS.PowerFlow.PowerFlowGS(ListBus,ListLine,Wbase);
    % [PowerFlow,~,~,~,~,~,~,~] = PowerFlowGS(ListBus,ListLine,Wbase);
    % This power flow algorithm is wrong but I do not find the weird part.                                      % ???
    % This power flow also does not match PowerFlowNR method.
elseif Flag_PowerFlowAlgorithm == 2
    [PowerFlow] = PowerFlowNR(ListBus,ListLine,Wbase);
end
% Form of PowerFlow{i}: P, Q, V, xi, w

% Move load flow (PLi and QLi) to bus admittance matrix
[ListBus,ListLineNew,PowerFlowNew] = SimplexPS.PowerFlow.Load2SelfBranch(ListBus,ListLine,PowerFlow);

% For printting later
ListPowerFlow = SimplexPS.PowerFlow.Rearrange(PowerFlow);
ListPowerFlowNew = SimplexPS.PowerFlow.Rearrange(PowerFlowNew);

% Update V and I
[V,I] = PowerFlowUpdateVI(PowerFlowNew);
% Notes:
% The codes in this part are borrowed from the SimplexPS toolbox. The V and
% I are updated based on the new power flow.

% ### For test
VoltageTheta = ListPowerFlowNew(:,5);
Max_VoltageThetaDiff = CalDiffMax(VoltageTheta);
Max_VoltageThetaDiff = Max_VoltageThetaDiff/pi*180