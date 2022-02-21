% This function re-arranges the netlist data of network lines.

% Author(s): Yitong Li

function [UpdateLine,N_Branch,N_Bus] = RearrangeListLine(UserData,ListBus)

%% Load data
ListLine     = xlsread(UserData,'NetworkLine');
ListLineIEEE = xlsread(UserData,'NetworkLine_IEEE');

%% IEEE Form to general form
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
        % which avoids the mutiple appearance of this self branch. This
        % requires to set the Rorig = Xorig = inf, which is different from
        % the normal network line form.
    end
    
  	% Set B, G for original branch
    B_orig = zeros(N_Br_orig,1);
    G_orig = B_orig;
    
    % Initialize the index for self-branch: normal order
    FB_self = transpose([1:N_Br_self]);
    TB_self = FB_self;
    
    % Initialize R, X, B, G by 0 for self branch
    R_self = inf(N_Br_self,1);
    X_self = R_self;
    B_self = zeros(N_Br_self,1);
    G_self = B_self;
    
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
        if ( isinf(R(i)) || isinf(X(i)) ) && (B(i)==0) && (G(i)==0)
            % The branch is open-circuit, jump this branch, i.e., delete
            % it.
            CountIndexDelete = CountIndexDelete + 1;
        else
            % The branch is NOT open-circuit, save this branch.
            UpdateLine_AfterDelete(i-CountIndexDelete,:) = UpdateLine(i,:);
        end
    end
    UpdateLine = UpdateLine_AfterDelete;
    
    % Use the IEEE form to overwrite the standard form.
    ListLine = UpdateLine;
end

%% Re-arrange netlist line
% Organize data
[N_Branch,ColumnMax_Line] = size(ListLine); 

% Replace NaN by inf
netlist_line_NaN = isnan(ListLine);
[r,c] = find(netlist_line_NaN == 1);  	% Find the index of "inf"
ListLine(r,c) = inf;

% Get data
FB  = ListLine(:,1);   % From bus
TB  = ListLine(:,2);   % To bus
Rbr = ListLine(:,3);
Xbr = ListLine(:,4);
Bbr = ListLine(:,5);
Gbr = ListLine(:,6);
Tbr = ListLine(:,7);

% Check data overflow
if (ColumnMax_Line>7)
    error(['Error: Line data overflow.']); 
end

% Check number of bus
N_Bus = max(max(FB), max(TB) );  

% Check short-circuit and open-circuit
for i = 1:N_Branch
    if ( isinf(Rbr(i)) || isinf(Xbr(i)) ) && ((Bbr(i)==0) && (Gbr(i)==0)) 
    	error(['Error: Branch' num2str(FB(i)) num2str(TB(i)) ' is open circuit']);
    end
    if ( (Rbr(i)==0) && (Xbr(i)==0) ) || isinf(Bbr(i)) || isinf(Gbr(i)) 
        error(['Error: Branch' num2str(FB(i)) num2str(TB(i)) ' is short circuit']);
    end
    if (Rbr(i)<0) || (Xbr(i)<0) || (Bbr(i)<0) || (Gbr(i)<0)
        error(['Error: Negative line paramters']);
    end
    if Tbr(i) <= 0
        error(['Error: Turns ratio can not be less than or equal to 0.']);
    end
end

% Check self-branch
for i = 1:N_Branch
    if FB(i) == TB(i)
        if Rbr(i)~=0 && (~isinf(Rbr(i))) 
            error(['Error: The self branch has to be a LCG parallel branch with R = 0 or a CG parallel branch with R = inf.'])
        end
    end
end

% Add inf inductive load to ListLine for future use
% ListLine = [ListLine,inf([N_Branch,1],'double')]; % Set all XL to inf defaultly
% ListLine = [ListLine,zeros([N_Branch,1])];

% Add area type into ListLine
for i = 1:N_Branch
    i1 = find(ListBus(:,1) == FB(i));
    i2 = find(ListBus(:,1) == TB(i));
    if ListBus(i1,12)~=ListBus(i2,12)
        error(['Error: The branch ' num2str(i1) '-' num2str(i2) ' is a hybrid branch, which is connected to ac at one side, and dc at the other side.']);
    else
        ListLine(i,8) = ListBus(i1,12);
    end
end

% Re-arrange the data
% Ensure "From Bus" <= "To Bus" is always valid
for i = 1:N_Branch
    if FB(i) > TB(i)
        % Switch "From" and "To"
        [TB(i),FB(i)] = deal(FB(i),TB(i));
        [ListLine(i,2),ListLine(i,1)] = deal(ListLine(i,1),ListLine(i,2));
        
        % Switch the positions of line impedance and transformer
        Rbr(i) = Rbr(i)*Tbr(i)^2;
        ListLine(i,3) = ListLine(i,3)*ListLine(i,7)^2;
        Xbr(i) = Xbr(i)*Tbr(i)^2;
        ListLine(i,4) = ListLine(i,4)*ListLine(i,7)^2;
        
        % Change the turns ratio
        Tbr(i) = 1/Tbr(i);
        ListLine(i,7) = 1/ListLine(i,7);
    end
end

% Re-order the branch sequence
ListLine = sortrows(ListLine,2);
ListLine = sortrows(ListLine,1);

% Output
UpdateLine = ListLine;

% The form of ListLine 
% 1          2        3   4    5    6   7   8
% From bus | To bus | R | wL | wC | G | T | Area type

end