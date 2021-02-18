% This function addes branches into simulink model.

% Author(s): Yitong Li

function [FullName_Branch,Name_Branch,Shift_ToBus] = ...
    SimAddBranch(Name_Model,Size_Branch,Shift_Branch,Pos_Bus,ListLine)

% Organize data
fb = ListLine(:,1); % From bus
tb = ListLine(:,2); % To bus
Rbr  = ListLine(:,3);
Xbr  = ListLine(:,4);
Bbr  = ListLine(:,5);
Gbr  = ListLine(:,6);
Tbr  = ListLine(:,7);
N_Branch = length(fb);

Count_ToBus = zeros(max(tb),1);

% Check if load data is combined into "ListLine"
[~,cmax_ListLine] = size(ListLine);
if cmax_ListLine>7
    XL = ListLine(:,8);
    Flag_LoadCombination = 1;
else
    XL = inf(size(ListLine(:,1)));
    Flag_LoadCombination = 0;
end

% Add branches
for i = 1:N_Branch
    
    % Initialize the postion of branch
    Pos_Branch{i} = [Pos_Bus{tb(i)}(1),Pos_Bus{fb(i)}(2)];
    Name_Branch{i} = ['Branch' num2str(fb(i)) '-' num2str(tb(i))];
    
    % Avoid the name duplication when multiple branches are connected in parallel between same two buses.
    for j = 1:(length(Name_Branch)-1)
        if fb(i)==fb(j) && tb(i)==tb(j)
            Name_Branch{i} = [Name_Branch{i},'_'];
        end
    end
    
    % Get the full name of the branch
    FullName_Branch{i} = [Name_Model '/' Name_Branch{i}];
    
    % ### Add self branch and load
    if fb(i) == tb(i)
        % Add block
        add_block(['powerlib/Elements/Three-Phase Parallel RLC Branch'],FullName_Branch{i});
        Pos_Branch{i} = Pos_Branch{i} + Shift_Branch;
        set_param(FullName_Branch{i},'position',[Pos_Branch{i},Pos_Branch{i}+Size_Branch]);
        set_param(FullName_Branch{i},'Orientation','down');
        set_param(FullName_Branch{i},'Measurements','None');
        
     	% Connect the floating terminals of self-branch to Y configuration
        add_line(Name_Model,...
            {[Name_Branch{i} '/Rconn2'],[Name_Branch{i} '/Rconn2']},...
            {[Name_Branch{i} '/Rconn1'],[Name_Branch{i} '/Rconn3']});  
        
        % Set branch type
        if (Flag_LoadCombination == 1) && (~isinf(XL(i)))
         	if (Gbr(i)==0) && (Bbr(i)==0)
                set_param(FullName_Branch{i},'BranchType','L');
            elseif Gbr(i)==0
                set_param(FullName_Branch{i},'BranchType','LC');
            elseif Bbr(i)==0
                set_param(FullName_Branch{i},'BranchType','RL');
            else
                set_param(FullName_Branch{i},'BranchType','RLC');
            end
            set_param(FullName_Branch{i},'Inductance',['(' num2str(XL(i)) ')*Zbase/Wbase']);
        else
            % Assume the self-branch is pure RC
           	if ~((Rbr(i)==0) && (Xbr(i)==0))
                error(['Error: the self branch contains L or R']);
            end
            if (Gbr(i)==0) && (Bbr(i)==0)
                error(['Error: open circuit']);
            elseif Gbr(i)==0        % Pure capacitance
                set_param(FullName_Branch{i},'BranchType','C');
            elseif Bbr(i)==0        % Pure resistance
                set_param(FullName_Branch{i},'BranchType','R');
            else                    % RC branch
                set_param(FullName_Branch{i},'BranchType','RC');
            end
        end
        
    	% Set customer data
      	set_param(FullName_Branch{i},'Capacitance',['(' num2str(Bbr(i)) ')*Ybase/Wbase']);
      	set_param(FullName_Branch{i},'Resistance',['(' num2str(1/Gbr(i)) ')*Zbase']);

        
    % ### Add mutual branch
    else
        % Check if transformer is added
        if Tbr(i) == 1
            Count_Trans = 0;    % No transformer
        else
            Count_Trans = 1;    % With transferomer
        end
        
        % Add mutual branch
        add_block(['powerlib/Elements/Three-Phase Series RLC Branch'],FullName_Branch{i});
        Count_ToBus(tb(i)) = Count_ToBus(tb(i)) + 1;
        Shift_ToBus{i} = Count_ToBus(tb(i));
        Pos_Branch{i} = Pos_Branch{i} + [Shift_Branch(1)*Shift_ToBus{i},Shift_Branch(2)*(Count_Trans+1)];
        set_param(FullName_Branch{i},'position',[Pos_Branch{i},Pos_Branch{i}+Size_Branch]);
        set_param(FullName_Branch{i},'Orientation','down');
        set_param(FullName_Branch{i},'Measurements','None');
        
        % Assume the mutual-branch is pure RL
        if ~(isinf(Gbr(i)) || isinf(Bbr(i)))
            error('Error: the mutual branch contains B or G');      
        end
        if (Rbr(i)==0) && (Xbr(i)==0)
            error('Error: short circuit')
        elseif Rbr(i)==0      % Pure inductance
            set_param(FullName_Branch{i},'BranchType','L');
        elseif Xbr(i)==0      % Pure resistance
            set_param(FullName_Branch{i},'BranchType','R');
        else                % RL branch
            set_param(FullName_Branch{i},'BranchType','RL');    
        end
        
        % Set customer data
        set_param(FullName_Branch{i},'Resistance',['(' num2str(Rbr(i)) ')*Zbase']);
     	set_param(FullName_Branch{i},'Inductance',['(' num2str(Xbr(i)) ')*Zbase/Wbase'])
    end
    
end

end