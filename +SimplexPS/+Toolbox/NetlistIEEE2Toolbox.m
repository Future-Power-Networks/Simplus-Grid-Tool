% This function transfers netlise of lines from IEEE form to toolbox form
% if required.

% Author(s): Yitong Li

function [UpdateLine] = NetlistIEEE2Toolbox(ListLine,ListLineIEEE)

if ListLineIEEE(1,1) == 1
    % "ListLineIEEE" is enabled, which converts "ListLineIEEE" to
    % "ListLine",
    fprintf('Warning: Line data in IEEE form is enabled, which over-writes the original line data.\n')
    
    ListLineIEEE = ListLineIEEE(4:end,:);
    
 	% Replace NaN by inf
    netlist_line_NaN = isnan(ListLineIEEE);
    [r,c] = find(netlist_line_NaN == 1);  	% Find the index of "inf"
    ListLineIEEE(r,c) = inf;
    
    % Organize original branch data
    FB_orig = ListLineIEEE(:,1);
    TB_orig = ListLineIEEE(:,2);
    R_orig = ListLineIEEE(:,3);
    X_orig = ListLineIEEE(:,4);
    B_pi = ListLineIEEE(:,5);
    G_pi = ListLineIEEE(:,6);
    T_orig = ListLineIEEE(:,7);
    
    % Number of bus
    N_Br_orig = length(FB_orig);              % Number of mutual branch
    N_Bus_ = max(max(FB_orig),max(TB_orig));
    N_Br_self = N_Bus_;                    % Number of self branch
    
	% Check self branch error
    for i = 1:N_Br_orig
        if FB_orig(i) == TB_orig(i)
            if ~( isinf(R_orig(i)) && isinf(X_orig(i)) )
                error(['Error: Branch ' num2str(i) ' is a self branch. Its R and X should be inf.']);
            elseif T_orig(i)~= 1
                error(['Error: Branch ' num2str(i) ' is a self branch. Its turn ratio should be 1.']);
            end
        end
        % For adding shunt B and G, the self-branch for IEEE form might
        % appear in its original data. But this self-branch will also be
        % treated as a pi-circuit mutual-branch next to calculate the new
        % self branch. The original self-branch will be deleted later,
        % which avoids the mutiple appearance of this self branch.
    end
    
  	% Set B, G for original branch
    B_orig = zeros(N_Br_orig,1);
    G_orig = inf(N_Br_orig,1);
    
    % Initialize the index for self-branch
    FB_self = transpose([1:N_Br_self]);
    TB_self = FB_self;
    
    % Initialize R, X, B, G by 0 for self branch
    R_self = zeros(N_Br_self,1);
    X_self = R_self;
    B_self = R_self;
    G_self = R_self;
    
    % Initialize turn ratio by 1
    T_self = ones(N_Br_self,1);
    
    % Move G+jB from pi-circuit to self-branch
    for i = 1:N_Br_orig
        B_self(FB_orig(i)) = B_self(FB_orig(i))+B_pi(i)/2;
        B_self(TB_orig(i)) = B_self(TB_orig(i))+B_pi(i)/2;
        
     	G_self(FB_orig(i)) = G_self(FB_orig(i))+G_pi(i)/2;
        G_self(TB_orig(i)) = G_self(TB_orig(i))+G_pi(i)/2;
    end
    
    % Get ListLine
    FB = [FB_orig;
          FB_self];
    TB = [TB_orig;
          TB_self];
    R = [R_orig;
         R_self];
    X = [X_orig;
         X_self];
    B = [B_orig;
         B_self];
    G = [G_orig;
         G_self];
    T = [T_orig;
         T_self];
    
    % Use IEEE form to over-write toolbox form 
    UpdateLine = [FB,TB,R,X,B,G,T];
    
    % Delete the open circuit branch
    N_Branch = N_Br_orig + N_Br_self;
    CountIndexDelete = 0;
    for i = 1:N_Branch
        % Find the open-circuit branch and delete it
        if ( isinf(R(i)) || isinf(X(i)) || ((B(i)==0)&&(G(i)==0)) )
            % The branch is open-circuit, jump this branch, i.e., delete
            % it.
            CountIndexDelete = CountIndexDelete + 1;
        else
            % The branch is NOT open-circuit, save this branch.
            UpdateLine_AfterDelete(i-CountIndexDelete,:) = UpdateLine(i,:);
        end
    end
    UpdateLine = UpdateLine_AfterDelete;
    
else
    % Use the toolbox line form
    UpdateLine = ListLine;
end

end