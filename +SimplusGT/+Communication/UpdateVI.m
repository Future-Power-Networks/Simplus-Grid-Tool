% This function updates V and I based on the new power flow data

% Author(s): Yitong Li

function [V, I] = UpdateVI(PowerFlow)

N_Bus = length(PowerFlow);

for i = 1:N_Bus

P(i) = PowerFlow{i}(1);
Q(i) = PowerFlow{i}(2);
Vm(i) = PowerFlow{i}(3);
Ang(i) = PowerFlow{i}(4);

V(i,1) = SimplusGT.pol2rect(Vm(i),Ang(i));
S(i) = P(i)+1i*Q(i);        % S in Load convention
I(i,1) = conj(S(i)/V(i));   % I in load convention
I(i,1) = -I(i,1);           % I in source convention

end