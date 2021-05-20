% This function gets the dimensions of state, input, output of a descriptor
% state space system or a state space system.

% Author(s): Yitong Li

%%  Notes:
%
% Descriptor state space system:
% E*dx/dt = A*x + B*u
% y       = C*x + D*u
%
% Normally, length(u) = column(D) = column(B), and length(y) = row(D) =
% row(C). However, this relation is not valid theoratically for static
% system with A=B=C=E=[] and D~=[], even though matlab can still calculate
% the right results sometimes because the empty matrix may also have
% non-zero dimension. Hence it is always better to use size(D) to get the
% dimensions of input and output vectors.
%
% Normally, length(x) = row(A) = column(A) = row(B) = column(C). However,
% this relation is not valid for a system without input and output but only
% state theoratically, i.e., B=C=D=[] but A~=[] and E~=[]. Hence it is
% always better to use length(A) to get the dimenstion of state vector.

%%
function [lx,lu,ly] = DssGetDim(Gdss)

if ~SimplusGT.is_dss(Gdss)
    error(['Error: The system is not in dss form.']);
end

lx = length(Gdss.A);
[ly,lu] = size(Gdss.D);

end