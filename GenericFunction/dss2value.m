% This function calculates the gain value of a descriptor state space
% syste at a given frequency point jw.

% Author(s): Yitong Li

%% Notes
%
% The system can be either in state space form or descriptor state space
% form; The system can be either in real state space form or complex state
% space form.

%%
function Gvalue = dss2value(G,jw);

% Get matrices
A = G.A; B = G.B; C = G.C; D = G.D; E = G.E;

% Check if descriptor state space form
if ( isempty(E) && (~isempty(A)) )
    E = eye(length(A));
end
    
% Calculate the gain value
Gvalue = C*inv(jw*E - A)*B + D;

end