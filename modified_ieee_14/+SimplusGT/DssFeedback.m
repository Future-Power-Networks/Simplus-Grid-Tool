% This function calculates the feedback of two descritor systems

% Author(s): Yitong Li

%% Notes
%
% Feedback block diagram
%                      +------+
%          v --------->|      |---------> z
%                      |  G1  |
%          u --->O---->|      |-----+---> y
%                |     +------+     |
%                |                  |
%                |     +------+     |
%                +-----|  G2  |<----+
%                      +------+  
%
%
% The vector FEEDIN contains indices into the input vector of G1 and
% specifies which inputs u are involved in the feedback loop. Similarly,
% FEEDOUT specifies which outputs y of G1 are used for feedback. If SIGN=1
% then positive feedback is used. If SIGN=-1 or SIGN is omitted, then
% negative feedback is used.
%   
% In all cases, the resulting model G has the same inputs and outputs as G1
% (with their order preserved), i.e.,
% G1:
% E1*d_x1/dt = A1*x1 + B1*u1
% y1         = C1*x1 + D1*u1
% G2:
% E2*d_x2/dt = A2*x2 + B2*u2
% y2         = C2*x2 + D2*u2
% -> feedback system G: 
% State equation
% [E1 0 ]*[d_x1]*dt = [A]*[x1] + [B]*[u1]
% [0  E2] [d_x2]          [x2]
% Output equation
% [y1]              = [D]*[x1] + [D]*[u1]
%                         [x2]
% In other words, the feedback operation does not influence matrix E.

%% Function

function G = DssFeedback(G1,G2,feedin,feedout,sign)

% Check the number of input arguments
narginchk(2,5);

% Check if descriptor system
if ((isempty(G1.E) && ~isempty(G1.A)) || (isempty(G2.E) && ~isempty(G2.A)))
    error(['Error: G1 or G2 is not in descriptor state space form']);
end

% Set default values
if nargin == 2
    [ly1,lu1] = size(G1.D);
    feedin = [1:lu1]; 
    feedout = [1:ly1];
end
if nargin == 3
    error(['Error: feedout is missing'])
end
if nargin <= 5
    sign = -1;      % Default, negative feedback
end

% Feedback
G = feedback(G1,G2,feedin,feedout,sign);
G.E = blkdiag(G1.E,G2.E);

end