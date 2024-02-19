% input: 1) whole-system admittance model 'Gsys' in state-space form
%        3) the selected mode (a single value)
% output: participation factors of the states to the selected mode, a vector.

function StatePF=StatePFCal(GsysSs, ModeSel)
    A=GsysSs.A;
    [Phi,~]=eig(A);
    Psi=inv(Phi); % left eigenvector
    for k=1:length(A)
        StatePF(k)=Psi(ModeSel,k) * Phi (k,ModeSel);
    end
end