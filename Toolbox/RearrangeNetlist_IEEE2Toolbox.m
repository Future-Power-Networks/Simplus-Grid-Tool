% This function transfers netlise of lines from IEEE form to toolbox form
% if required.

% Author(s): Yitong Li

function [UpdateLine] = RearrangeNetlist_IEEE2Toolbox(ListLine,ListLineIEEE)

if ListLineIEEE(1,1) == 1
    % "ListLineIEEE" is enabled, which over-writes "ListLine"
    fprintf('Warning: Line data in IEEE form is enabled, which over-writes the original line data')
    ListLineIEEE = ListLineIEEE(4:end,:);
    
    % Organize mutual branch data
    FB_m = ListLineIEEE(:,1);
    TB_m = ListLineIEEE(:,2);
    R_m = ListLineIEEE(:,3);
    X_m = ListLineIEEE(:,4);
    B_pi = ListLineIEEE(:,5);
    G_pi = ListLineIEEE(:,6);
    T_m = ListLineIEEE(:,7);
    
    % Number of bus
    N_Br_m = length(FB_m);              % Number of mutual branch
    N_Bus_ = max(max(FB_m),max(TB_m));
    N_Br_s = N_Bus_;                    % Number of self branch
    
	% Check self branch error
    for i = 1:N_Br_m
        if FB_m(i) == TB_m(i)
            error(['Error: Branch ' num2str(FB_m(i)) ' is a self branch. "Netlist IEEE" should not have self branch.'])
        end
    end
    
    % Get the self branch data
    FB_s = transpose([1:N_Br_s]);
    TB_s = FB_s;
    
    % Get B, G for mutual branch
    B_m = zeros(N_Br_m,1);
    G_m = inf(N_Br_m,1);
    
    % Initialize R, X, B, G for self branch
    R_s = zeros(N_Br_s,1);
    X_s = R_s;
    B_s = R_s;
    G_s = R_s;
    T_s = ones(N_Br_s,1);
    
    for i = 1:N_Br_m
        B_s(FB_m(i)) = B_s(FB_m(i))+B_pi(i);
        B_s(TB_m(i)) = B_s(TB_m(i))+B_pi(i);
     	G_s(FB_m(i)) = G_s(FB_m(i))+G_pi(i);
        G_s(TB_m(i)) = G_s(TB_m(i))+G_pi(i);
    end
    
    % Get ListLine
    FB = [FB_m;FB_s];
    TB = [TB_m;TB_s];
    R = [R_m;R_s];
    X = [X_m;X_s];
    B = [B_m;B_s];
    G = [G_m;G_s];
    T = [T_m;T_s];
    
    % Use IEEE form to over-write toolbox form 
    UpdateLine = [FB,TB,R,X,B,G,T];
else
    % Use the toolbox form
    UpdateLine = ListLine;
end

end