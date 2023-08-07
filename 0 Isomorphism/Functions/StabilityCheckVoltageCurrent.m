% Set the default values
StableVoltageNode = 1;
StableCurrentNode = 1;

% Check voltage nodes
if ~isempty(zeta_m_V)
    if min(zeta_m_V)>=sigma_V_max
        StableVoltageNode = 1;
    else
        StableVoltageNode = 0;
    end
end

% Check current nodes
if ~isempty(zeta_m_I)
    if min(zeta_m_I)>=sigma_I_max
        StableCurrentNode = 1;
    else
        StableCurrentNode = 0;
    end
end

% Output
if StableVoltageNode==1 && StableCurrentNode==1
    fprintf('Stable!\n')
else
 	if StableVoltageNode==0
        fprintf('Unstable voltage node!\n')
    end
    if StableCurrentNode==0
        fprintf('Unstable current node!\n')
    end
end