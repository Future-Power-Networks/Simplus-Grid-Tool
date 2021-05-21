% This function connects device to bus

% Author(s): Yitong Li

function SimConnectDevice2Bus(Name_Model,Name_Bus,Name_Device,DeviceBus,DeviceType)

% Organize data
N_Device = length(DeviceType);

% Add block
for i = 1:N_Device
    
    % Get the bus index of the device
    Bus = DeviceBus{i};
    
    % If the device is NOT a "floating bus"
    if DeviceType{i}<=90
        % For ac apparatuses
        add_line(Name_Model,...
            {[Name_Device{i} '/Lconn1'],[Name_Device{i} '/Lconn2'],[Name_Device{i} '/Lconn3']},...
            {[Name_Bus{Bus} '/Lconn1'],[Name_Bus{Bus} '/Lconn2'],[Name_Bus{Bus} '/Lconn3']},...
            'autorouting','smart');
    elseif 1000<=DeviceType{i} && DeviceType{i}<=1090  
        % For dc apparatuses
    	add_line(Name_Model,...
            {[Name_Device{i} '/Lconn1'],[Name_Device{i} '/Lconn2']},...
            {[Name_Bus{Bus} '/Lconn1'],[Name_Bus{Bus} '/Lconn2']},...
            'autorouting','smart');
    elseif 2000<=DeviceType{i} && DeviceType{i}<=2090
        % For interlink apparatuses
        % Connect to ac bus
        add_line(Name_Model,...
            {[Name_Device{i} '/Lconn1'],[Name_Device{i} '/Lconn2'],[Name_Device{i} '/Lconn3']},...
            {[Name_Bus{Bus(1)} '/Lconn1'],[Name_Bus{Bus(1)} '/Lconn2'],[Name_Bus{Bus(1)} '/Lconn3']},...
            'autorouting','smart');
        % Connect to dc bus
     	add_line(Name_Model,...
            {[Name_Device{i} '/Lconn5'],[Name_Device{i} '/Lconn6']},...
            {[Name_Bus{Bus(2)} '/Lconn1'],[Name_Bus{Bus(2)} '/Lconn2']},...
            'autorouting','smart');
    end
    
end

end