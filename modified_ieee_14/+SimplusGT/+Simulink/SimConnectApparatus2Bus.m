% This function connects apparatus to bus

% Author(s): Yitong Li

function SimConnectApparatus2Bus(Name_Model,Name_Bus,Name_Apparatus,ApparatusBus,ApparatusType)

% Organize data
N_Apparatus = length(ApparatusType);

% Add block
for i = 1:N_Apparatus
    
    % Get the bus index of the apparatus
    Bus = ApparatusBus{i};
    
    % If the apparatus is NOT a "floating bus"
    if ApparatusType{i}<=90
        % For ac apparatuses
        add_line(Name_Model,...
            {[Name_Apparatus{i} '/Lconn1'],[Name_Apparatus{i} '/Lconn2'],[Name_Apparatus{i} '/Lconn3']},...
            {[Name_Bus{Bus} '/Lconn1'],[Name_Bus{Bus} '/Lconn2'],[Name_Bus{Bus} '/Lconn3']},...
            'autorouting','smart');
    elseif 1000<=ApparatusType{i} && ApparatusType{i}<=1090  
        % For dc apparatuses
    	add_line(Name_Model,...
            {[Name_Apparatus{i} '/Lconn1'],[Name_Apparatus{i} '/Lconn2']},...
            {[Name_Bus{Bus} '/Lconn1'],[Name_Bus{Bus} '/Lconn2']},...
            'autorouting','smart');
    elseif 2000<=ApparatusType{i} && ApparatusType{i}<=2090
        % For interlink apparatuses
        % Connect to ac bus
        add_line(Name_Model,...
            {[Name_Apparatus{i} '/Lconn1'],[Name_Apparatus{i} '/Lconn2'],[Name_Apparatus{i} '/Lconn3']},...
            {[Name_Bus{Bus(1)} '/Lconn1'],[Name_Bus{Bus(1)} '/Lconn2'],[Name_Bus{Bus(1)} '/Lconn3']},...
            'autorouting','smart');
        % Connect to dc bus
     	add_line(Name_Model,...
            {[Name_Apparatus{i} '/Lconn5'],[Name_Apparatus{i} '/Lconn6']},...
            {[Name_Bus{Bus(2)} '/Lconn1'],[Name_Bus{Bus(2)} '/Lconn2']},...
            'autorouting','smart');
    end
    
end

end