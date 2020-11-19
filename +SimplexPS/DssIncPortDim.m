% This function increases the dimensions of input and/or output of a
% descriptor (implicit) state space system.

% Author(s): Yitong Li

%%
% Example 1:
% Increase the tail of u
% State equation
% x = A*x + [B B']*[u ]
%                  [u']
% Output equation
% y = C*x + [D D']*[u ]
%                  [u']
% with B'=0 and D'=0

% Example 2:
% Increase the tail of y
% State equation
% x = A*x + B*u
% Output equation
% [y ] = [C ]*x + [D ]*u
% [y']   [C']     [D']
% with C'=0 and D'=0

%%
function G_Inc = DssIncPortDim(G,IncDim_u,IncDim_y,IncHeadTail)

    A = G.A; B = G.B; C = G.C; D = G.D; E = G.E;
    
    % Check if the system is in descriptor state space form
    if ismpty(E) && ~isempty(A)
        errror(['Error: the system is not in descriptor state space form']);
    end
    
    % Add Head, i.e., u_new = [u_add;u], y_new = [y_add;y]
    if IncHeadTail == 1
     	% Increase input dimension
        if IncDim_u > 0
            if ~isempty(B)  % Check for static system
                [r_B,~] = size(B);
                B = [zeros(r_B,IncDim_u),B];
            end
            [r_D,~] = size(D);
            D = [zeros(r_D,IncDim_u),D];
        end
        
        % Increase output dimension
        if IncDim_y > 0
            if ~isempty(C)
                [~,c_C] = size(C);
                C = [zeros(IncDim_y,c_C);C];
            end
            [~,c_D] = size(D);
            D = [zeros(IncDim_y,c_D);D];
        end
        
    % Add Tail, i.e., u_new = [u;u_add], y_new = [y_add;y]
    elseif IncHeadTail == 2
        % Increase input dimension
        if IncDim_u > 0
            if ~isempty(B)  % Check for static system
                [r_B,~] = size(B);
                B = [B,zeros(r_B,IncDim_u)];
            end
            [r_D,~] = size(D);
            D = [D,zeros(r_D,IncDim_u)];
        end
        
        % Increase output dimension
        if IncDim_y > 0
            if ~isempty(C)
                [~,c_C] = size(C);
                C = [C;zeros(IncDim_y,c_C)];
            end
            [~,c_D] = size(D);
            D = [D;zeros(IncDim_y,c_D)];
        end
    end
    
    G_Inc = dss(A,B,C,D,E);
    
end