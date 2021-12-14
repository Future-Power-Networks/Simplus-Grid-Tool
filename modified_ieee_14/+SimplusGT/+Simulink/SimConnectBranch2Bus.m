% This function connects branches to buses.

% Author(s): Yitong Li

function SimConnectBranch2Bus(Name_Model,Name_Bus,Name_Branch,Name_Trans,ListLine)

% Organize data
fb = ListLine(:,1);     % From bus
tb = ListLine(:,2);     % To bus
Tbr = ListLine(:,7);    % Transformer tap ratio
AreaType = ListLine(:,8);

N_Branch = length(fb);  % Number fo branches

% Connect mutual-branch to the output port of bus
for i = 1:N_Branch
    if fb(i) ~= tb(i)
        
        From = fb(i);
        To = tb(i);
        
        % ### Ac branch
        if AreaType(i) == 1
            % Connect branch to "to bus"
            add_line(Name_Model,...
                {[Name_Bus{To} '/Rconn1'],[Name_Bus{To} '/Rconn2'],[Name_Bus{To} '/Rconn3']},...
                {[Name_Branch{i} '/Rconn1'],[Name_Branch{i} '/Rconn2'],[Name_Branch{i} '/Rconn3']},...
                'autorouting','smart');
            
            if Tbr(i) == 1
                % Connect branch to "from bus"
                add_line(Name_Model,...
                    {[Name_Bus{From} '/Rconn1'],[Name_Bus{From} '/Rconn2'],[Name_Bus{From} '/Rconn3']},...
                    {[Name_Branch{i} '/Lconn1'],[Name_Branch{i} '/Lconn2'],[Name_Branch{i} '/Lconn3']},...
                    'autorouting','smart');
            else
                % Connect branch to transformer
                add_line(Name_Model,...
                    {[Name_Trans{i} '/Rconn1'],[Name_Trans{i} '/Rconn2'],[Name_Trans{i} '/Rconn3']},...
                    {[Name_Branch{i} '/Lconn1'],[Name_Branch{i} '/Lconn2'],[Name_Branch{i} '/Lconn3']},...
                    'autorouting','smart');
                % Connect transformer to "from bus"
                add_line(Name_Model,...
                    {[Name_Bus{From} '/Rconn1'],[Name_Bus{From} '/Rconn2'],[Name_Bus{From} '/Rconn3']},...
                    {[Name_Trans{i} '/Lconn1'],[Name_Trans{i} '/Lconn2'],[Name_Trans{i} '/Lconn3']},...
                    'autorouting','smart');
            end
            
        % ### Dc branch
        elseif AreaType(i) == 2
            % Connect branch to "to bus"
            add_line(Name_Model,...
                {[Name_Bus{To} '/Rconn1'],[Name_Bus{To} '/Rconn2']},...
                {[Name_Branch{i} '/Rconn1'],[Name_Branch{i} '/Rconn2']},...
                'autorouting','smart');
            % Connect branch to "from bus"
            add_line(Name_Model,...
                    {[Name_Bus{From} '/Rconn1'],[Name_Bus{From} '/Rconn2']},...
                    {[Name_Branch{i} '/Lconn1'],[Name_Branch{i} '/Lconn2']},...
                    'autorouting','smart');
                
        % ### Error
        else
            error(['Error.']);
                
        end

    end
end

% Connect self-branch to the output port of bus
% Remarks: Doing this procedure seperately with the mutual-branch procedure
% because this will lead to a better the autorouting quality of simulink.
for i = 1:N_Branch
  	if fb(i) == tb(i)
        
        From = fb(i);     % Index number of the bus
        
        % ### Ac branch
        if AreaType(i) == 1
            % Connect self-branch to bus
            add_line(Name_Model,...
                {[Name_Bus{From} '/Rconn1'],[Name_Bus{From} '/Rconn2'],[Name_Bus{From} '/Rconn3']},...
                {[Name_Branch{i} '/Lconn1'],[Name_Branch{i} '/Lconn2'],[Name_Branch{i} '/Lconn3']},...
                'autorouting','smart');
        
        % ### Dc branch
        elseif AreaType(i) == 2
            % Connect self-branch to bus
            add_line(Name_Model,...
                {[Name_Bus{From} '/Rconn1'],[Name_Bus{From} '/Rconn2']},...
                {[Name_Branch{i} '/Rconn1'],[Name_Branch{i} '/Lconn1']},...
                'autorouting','smart');
        
        end
        
    end
end

end