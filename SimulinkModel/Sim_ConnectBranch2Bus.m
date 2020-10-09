% This function connects branches to buses.

% Author(s): Yitong Li

function Sim_ConnectBranch2Bus(Name_Model,Name_Bus,Name_Branch,ListLine)

% Organize data
fb = ListLine(:,1);     % From bus
tb = ListLine(:,2);     % To bus
N_Branch = length(fb);  % Number fo branches

% Connect mutual-branch to the output port of bus
for i = 1:N_Branch
    if fb(i) ~= tb(i)
        
        From = fb(i);
        To = tb(i);
        
        % Connect branch and "from bus"
        add_line(Name_Model,...
            {[Name_Bus{From} '/Rconn1'],[Name_Bus{From} '/Rconn2'],[Name_Bus{From} '/Rconn3']},...
            {[Name_Branch{i} '/Lconn1'],[Name_Branch{i} '/Lconn2'],[Name_Branch{i} '/Lconn3']},...
            'autorouting','smart');
        
        % Connect branch and "to bus"
     	add_line(Name_Model,...
            {[Name_Bus{To} '/Rconn1'],[Name_Bus{To} '/Rconn2'],[Name_Bus{To} '/Rconn3']},...
            {[Name_Branch{i} '/Rconn1'],[Name_Branch{i} '/Rconn2'],[Name_Branch{i} '/Rconn3']},...
            'autorouting','smart');
    end
end

% Connect self-branch to the output port of bus
% Remarks: Doing this procedure seperately with the mutual-branch procedure
% because this will lead to a better the autorouting quality of simulink.
for i = 1:N_Branch
  	if fb(i) == tb(i)
        
        From = fb(i);     % Index number of the bus
        
        % Connect self-branch to bus
        add_line(Name_Model,...
            {[Name_Bus{From} '/Rconn1'],[Name_Bus{From} '/Rconn2'],[Name_Bus{From} '/Rconn3']},...
            {[Name_Branch{i} '/Lconn1'],[Name_Branch{i} '/Lconn2'],[Name_Branch{i} '/Lconn3']},...
            'autorouting','smart');
    end
end

end